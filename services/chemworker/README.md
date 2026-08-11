# chemmed-chemworker

gRPC service owning RDKit and the in-memory fingerprint index. The Rust API
delegates all structure handling here.

## Why this is a separate service

RDKit has no mature Rust equivalent, so structure parsing, standardisation, and
fingerprinting have to happen in Python. Keeping them behind a gRPC boundary
means the API layer never has to embed a Python runtime, and the ~600 MB
fingerprint index is loaded once per worker rather than once per request.

## Setup

```bash
uv venv --python 3.12 .venv && uv pip install --python .venv -e ".[dev]"
```

Generate the protobuf stubs (required before first run, and after any change to
`proto/chemworker.proto`):

```bash
.venv/Scripts/python -m grpc_tools.protoc -I proto --python_out=src/chemworker --grpc_python_out=src/chemworker proto/chemworker.proto
```

The generated `chemworker_pb2*.py` files are git-ignored — they are build output.

## Run

Requires a built fingerprint index (`chemmed-ingest index`):

```bash
.venv/Scripts/chemmed-chemworker
```

Listens on `[::]:50051` by default; override with `CHEMWORKER_BIND`.

## Contract notes

- **A bad SMILES is not a gRPC error.** `Standardize` returns `ok=false` with a
  populated `error` field, because unparseable user input is a client-side
  condition the UI shows inline. `SimilaritySearch` does abort with
  `INVALID_ARGUMENT`, since there is no partial result to return.
- **`fp_signature` in `Status` must match the ingest config.** A mismatch means
  the index was built with different fingerprint parameters and every score
  would be silently wrong. The index refuses to load in that case.
