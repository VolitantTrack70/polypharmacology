@echo off
REM Start the Rust API.
REM
REM cargo and protoc are prepended to PATH because a freshly-installed
REM toolchain is not always visible to spawned processes until re-login.
REM
REM Must run from services/api: main.rs loads ../../.env relative to cwd.

setlocal
set "PATH=%USERPROFILE%\.cargo\bin;%ProgramFiles%\protoc\bin;%PATH%"
cd /d "%~dp0..\services\api" || exit /b 1

where cargo >nul 2>nul
if errorlevel 1 (
    echo [api] cargo not found. Install with: winget install Rustlang.Rustup
    exit /b 1
)

cargo run %*
