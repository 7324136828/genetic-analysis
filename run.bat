@echo off
rem Start GenomeScope: backend (FastAPI) + frontend (Vite) together.
rem   run.bat              development servers
rem   run.bat --build      production build served by vite preview
rem   run.bat --no-browser do not open the browser
rem Ports: first free port of each list in port.json. Ctrl+C stops both servers.
setlocal
cd /d "%~dp0"

if not exist ".venv\Scripts\python.exe" (
    echo [RUN.BAT] Virtual environment not found. Running setup first ...
    call setup.bat
    if errorlevel 1 exit /b 1
)
if not exist "frontend\node_modules" (
    echo [RUN.BAT] Frontend dependencies missing. Running setup first ...
    call setup.bat --skip-backend
    if errorlevel 1 exit /b 1
)

".venv\Scripts\python.exe" launch.py %*
if errorlevel 1 pause
endlocal
