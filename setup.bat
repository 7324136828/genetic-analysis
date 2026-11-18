@echo off
rem GenomeScope setup (Windows): creates .venv, installs backend + frontend dependencies.
setlocal
cd /d "%~dp0"

set "PY="
where py >nul 2>nul && set "PY=py -3"
if not defined PY (
    where python >nul 2>nul && set "PY=python"
)
if not defined PY (
    echo [SETUP.BAT] Python 3.10+ was not found on PATH. Install it from https://www.python.org/downloads/
    pause
    exit /b 1
)

echo [SETUP.BAT] Dispatching setup.py with %PY% ...
%PY% setup.py %*
if errorlevel 1 (
    echo [SETUP.BAT] Setup failed.
    pause
    exit /b 1
)
echo [SETUP.BAT] Setup complete. Start the app with run.bat
endlocal
