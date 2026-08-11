@echo off
REM Start the gRPC cheminformatics worker.
REM
REM Requires a built fingerprint index. If it is missing, run:
REM     services\ingest\.venv\Scripts\chemmed-ingest index
REM
REM The generated protobuf stubs are git-ignored build output. If the import
REM fails, regenerate them:
REM     services\chemworker\.venv\Scripts\python services\chemworker\gen_proto.py

setlocal
cd /d "%~dp0..\services\chemworker" || exit /b 1

if not exist ".venv\Scripts\python.exe" (
    echo [chemworker] venv missing. Run:
    echo     uv venv --python 3.12 .venv ^&^& uv pip install --python .venv -e ".[dev]"
    exit /b 1
)

if not exist "src\chemworker\chemworker_pb2_grpc.py" (
    echo [chemworker] protobuf stubs missing, generating...
    .venv\Scripts\python.exe gen_proto.py || exit /b 1
)

.venv\Scripts\python.exe -m chemworker.server %*
