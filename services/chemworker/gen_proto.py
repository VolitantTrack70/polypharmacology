"""Generate the protobuf/gRPC stubs, then fix grpcio's broken import.

Run with the chemworker venv:

    .venv/Scripts/python gen_proto.py

WHY THIS SCRIPT EXISTS
----------------------
`grpc_tools.protoc` emits `import chemworker_pb2` into the generated
`_pb2_grpc.py` -- a flat, top-level import. That only resolves if the output
directory happens to be on sys.path, which it is not when the stubs live inside
a package. Importing `chemworker.server` then dies with ModuleNotFoundError.

This is a long-standing grpcio codegen issue (protocolbuffers/protobuf#1491).
The fix is to rewrite the import to be package-relative after generation.
The stubs are build output and git-ignored, so this must be re-runnable.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent
PROTO_DIR = HERE / "proto"
OUT_DIR = HERE / "src" / "chemworker"
PROTO = PROTO_DIR / "chemworker.proto"


def generate() -> None:
    cmd = [
        sys.executable,
        "-m",
        "grpc_tools.protoc",
        f"-I{PROTO_DIR}",
        f"--python_out={OUT_DIR}",
        f"--grpc_python_out={OUT_DIR}",
        f"--pyi_out={OUT_DIR}",
        str(PROTO),
    ]
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


def fix_imports() -> None:
    """Rewrite `import X_pb2 as Y` to `from chemworker import X_pb2 as Y`."""
    target = OUT_DIR / "chemworker_pb2_grpc.py"
    source = target.read_text(encoding="utf-8")

    patched, n = re.subn(
        r"^import (\w+_pb2) as (\w+)$",
        r"from chemworker import \1 as \2",
        source,
        flags=re.MULTILINE,
    )

    if n == 0 and "from chemworker import" not in source:
        raise SystemExit(
            f"Could not find the flat import in {target.name}. grpcio's codegen "
            f"may have changed -- inspect the file and update this script."
        )

    target.write_text(patched, encoding="utf-8")
    print(f"patched {n} import(s) in {target.name}")


if __name__ == "__main__":
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    generate()
    fix_imports()
    print("done")
