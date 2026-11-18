"""
Allele normalisation and equivalence helpers shared by the NCBI integration, the
legacy MyVariant seeding path and the ClinVar release scan.

Everything here is pure (no I/O) so the single definition of "is this the same
allele?" can be unit-tested and reused. Coordinates are explicit about their base:
`pos0` arguments are 0-based (SPDI), `pos` arguments are 1-based (VCF).
"""

from __future__ import annotations

import re
from typing import List, Optional, Sequence, Tuple

# Shown when the sample has a variant call at the locus, but not the allele of the
# ClinVar / NCBI record being displayed. The wording is part of the UI contract.
DIFFERENT_ALT_STATUS = "Different ALT Allele (Not the NCBI Variant)"
# The record's own alleles are unusable (e.g. a placeholder 'N'), so the call cannot
# be attributed to it. Never counted as carried.
UNVERIFIED_ALLELE_STATUS = "Allele Not Verified"

_ACGT = re.compile(r"[ACGT]+")


def is_plain_allele(seq: Optional[str]) -> bool:
    """True for a non-empty A/C/G/T string (no '.', '*', '<SYMBOLIC>', 'N')."""
    return bool(seq) and _ACGT.fullmatch(seq.upper()) is not None


def minimal_alleles(pos0: int, deleted: str, inserted: str) -> Tuple[int, str, str]:
    """Trim common suffix then common prefix (left-most minimal form). Returns (pos0, del, ins)."""
    d, i = deleted, inserted
    while d and i and d[-1] == i[-1]:
        d, i = d[:-1], i[:-1]
    shift = 0
    while d and i and d[0] == i[0]:
        d, i = d[1:], i[1:]
        shift += 1
    return (pos0 + shift, d, i)


def normalize_vcf_allele(pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
    """
    VCF (1-based POS, REF, ALT) -> minimal (pos, deleted, inserted), either side possibly
    empty. Two VCF records that differ only in their anchor base(s) normalise to the same
    triple, which makes it a safe exact-match key. It does NOT shift an event through a
    repeat; use `alleles_equivalent` for that.
    """
    p0, d, i = minimal_alleles(pos - 1, ref.upper(), alt.upper())
    return (p0 + 1, d, i)


def alleles_equivalent(
    spdi_pos0: int, spdi_del: str, spdi_ins: str, vcf_pos: int, vcf_ref: str, vcf_alt: str
) -> Optional[bool]:
    """
    Decide whether an SPDI allele and a VCF allele describe the same haplotype, without
    needing the reference genome: the reference over the union window is reconstructed
    from the SPDI deleted sequence and the VCF REF. Returns None if the window has gaps
    (cannot decide) or the two reference sequences disagree.
    """
    if not vcf_ref or vcf_alt in (None, "", ".", "*") or vcf_alt.startswith("<"):
        return False
    v0 = vcf_pos - 1
    vref, valt = vcf_ref.upper(), vcf_alt.upper()
    lo = min(spdi_pos0, v0)
    hi = max(spdi_pos0 + len(spdi_del), v0 + len(vref))
    if hi - lo > 5000:
        return None
    ref: List[Optional[str]] = [None] * (hi - lo)
    for start, seq in ((spdi_pos0, spdi_del), (v0, vref)):
        for k, base in enumerate(seq):
            idx = start - lo + k
            if ref[idx] is not None and ref[idx] != base:
                return None
            ref[idx] = base
    if any(b is None for b in ref):
        return None
    ref_s = "".join(ref)  # type: ignore[arg-type]

    def apply(start: int, length: int, ins: str) -> str:
        o = start - lo
        return ref_s[:o] + ins + ref_s[o + length:]

    return apply(v0, len(vref), valt) == apply(spdi_pos0, len(spdi_del), spdi_ins)


def vcf_alleles_equivalent(
    pos_a: int, ref_a: str, alt_a: str, pos_b: int, ref_b: str, alt_b: str
) -> Optional[bool]:
    """`alleles_equivalent` for two VCF-style alleles (both 1-based). True / False / None."""
    return alleles_equivalent(pos_a - 1, ref_a.upper(), alt_a.upper(), pos_b, ref_b, alt_b)


def spans_touch(pos_a: int, ref_a: str, pos_b: int, ref_b: str) -> bool:
    """True when two VCF REF spans overlap or are directly adjacent (a precondition for
    `alleles_equivalent` to be able to reconstruct a contiguous reference window)."""
    a_lo, a_hi = pos_a, pos_a + len(ref_a)
    b_lo, b_hi = pos_b, pos_b + len(ref_b)
    return a_lo <= b_hi and b_lo <= a_hi


def find_alt_index(
    row_pos: int, row_ref: str, row_alts: Sequence[str], cv_pos: int, cv_ref: str, cv_alt: str
) -> Tuple[Optional[int], str]:
    """
    Which ALT of a VCF row (1-based index) is the ClinVar allele (cv_pos, cv_ref, cv_alt)?

    Returns (index, how):
      (k, "exact")       identical after normalisation
      (k, "equivalent")  same haplotype under a different indel representation
      (None, "none")     the row's ALTs are decidedly a different allele
      (None, "undecidable")  overlapping indel representations that cannot be compared
                              without the reference sequence, or the ClinVar allele is unusable
    """
    if not (is_plain_allele(cv_ref) and is_plain_allele(cv_alt)):
        return (None, "undecidable")
    want = normalize_vcf_allele(cv_pos, cv_ref, cv_alt)
    undecided = False
    for k, a in enumerate(row_alts, start=1):
        if not is_plain_allele(a) or not is_plain_allele(row_ref):
            continue
        if normalize_vcf_allele(row_pos, row_ref, a) == want:
            return (k, "exact")
    for k, a in enumerate(row_alts, start=1):
        if not is_plain_allele(a) or not is_plain_allele(row_ref):
            continue
        if len(row_ref) == len(a) and len(cv_ref) == len(cv_alt):
            continue  # substitutions have a single representation: exact match or different
        if not spans_touch(row_pos, row_ref, cv_pos, cv_ref):
            continue
        eq = vcf_alleles_equivalent(cv_pos, cv_ref, cv_alt, row_pos, row_ref, a)
        if eq:
            return (k, "equivalent")
        if eq is None:
            undecided = True
    return (None, "undecidable" if undecided else "none")
