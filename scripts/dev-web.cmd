@echo off
REM Start the SvelteKit dev server.
REM
REM Two Windows-specific workarounds are baked in here:
REM   1. Node is prepended to PATH -- spawned processes don't always inherit
REM      the machine PATH updated by the installer until a full re-login.
REM   2. We cd into services/web rather than passing a root argument, because
REM      SvelteKit overrides Vite's `root` with process.cwd().

setlocal
set "PATH=%ProgramFiles%\nodejs;%PATH%"
cd /d "%~dp0..\services\web" || exit /b 1
node node_modules\vite\bin\vite.js dev --port 5173 --strictPort %*
