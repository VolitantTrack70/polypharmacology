"""gRPC cheminformatics service.

Holds the fingerprint index in memory and answers similarity queries. This is
the only process that touches RDKit at request time.

Generate the protobuf stubs before first run:

    python -m grpc_tools.protoc -I proto \
        --python_out=src/chemworker --grpc_python_out=src/chemworker \
        proto/chemworker.proto
"""

from __future__ import annotations

import logging
import os
import signal
import time
from concurrent import futures
from pathlib import Path

import grpc
import psycopg
from grpc_health.v1 import health, health_pb2, health_pb2_grpc

from chemmed_ingest.chem.fingerprint import (
    StandardizationError,
    fingerprint_bytes,
    standardize,
)
from chemmed_ingest.chem.similarity import FingerprintIndex
from chemmed_ingest.config import DATABASE_URL, FP, PATHS

from chemworker import chemworker_pb2 as pb
from chemworker import chemworker_pb2_grpc as pb_grpc

log = logging.getLogger("chemworker")

MAX_LIMIT = 5000


class ChemWorkerServicer(pb_grpc.ChemWorkerServicer):
    def __init__(self, index: FingerprintIndex) -> None:
        self._index = index

    # -- Standardize ---------------------------------------------------------

    def Standardize(self, request, context):  # noqa: N802 (gRPC naming)
        try:
            std = standardize(request.smiles)
        except StandardizationError as exc:
            # A bad SMILES is a client error, not a server failure, so it comes
            # back as a populated `error` field rather than a gRPC status.
            return pb.StandardizeResponse(ok=False, error=str(exc))

        return pb.StandardizeResponse(
            ok=True,
            canonical_smiles=std.canonical_smiles,
            inchikey=std.inchikey,
            parent_inchikey=std.parent_inchikey,
        )

    # -- SimilaritySearch ----------------------------------------------------

    def SimilaritySearch(self, request, context):  # noqa: N802
        started = time.perf_counter()

        which = request.WhichOneof("query")
        if which == "smiles":
            try:
                std = standardize(request.smiles)
            except StandardizationError as exc:
                context.abort(grpc.StatusCode.INVALID_ARGUMENT, str(exc))
                return None
            query_fp = fingerprint_bytes(std.mol)
            canonical = std.canonical_smiles
        elif which == "chembl_id":
            row = self._lookup_fingerprint(request.chembl_id)
            if row is None:
                context.abort(
                    grpc.StatusCode.NOT_FOUND,
                    f"no fingerprint stored for {request.chembl_id}",
                )
                return None
            query_fp, canonical = row
        else:
            context.abort(grpc.StatusCode.INVALID_ARGUMENT, "query must set smiles or chembl_id")
            return None

        threshold = request.threshold or 0.40
        if not 0.0 < threshold <= 1.0:
            context.abort(grpc.StatusCode.INVALID_ARGUMENT, f"bad threshold {threshold}")
            return None

        limit = min(request.limit or 250, MAX_LIMIT)
        hits = self._index.search(query_fp, threshold=threshold, limit=limit)

        return pb.SimilaritySearchResponse(
            hits=[pb.SimilarityHit(chembl_id=h.chembl_id, tanimoto=h.tanimoto) for h in hits],
            canonical_smiles=canonical,
            searched=len(self._index),
            considered=len(hits),
            elapsed_ms=(time.perf_counter() - started) * 1000.0,
        )

    # -- Status --------------------------------------------------------------

    def Status(self, request, context):  # noqa: N802
        return pb.StatusResponse(
            index_loaded=True,
            compound_count=len(self._index),
            fp_signature=FP.signature,
            version="0.1.0",
        )

    # -- helpers -------------------------------------------------------------

    @staticmethod
    def _lookup_fingerprint(chembl_id: str) -> tuple[bytes, str] | None:
        with psycopg.connect(DATABASE_URL) as conn, conn.cursor() as cur:
            cur.execute(
                """
                SELECT f.fp, c.canonical_smiles
                FROM chem.compound_fingerprint f
                JOIN chem.compound c USING (chembl_id)
                WHERE f.chembl_id = %s
                """,
                (chembl_id,),
            )
            row = cur.fetchone()
        if row is None:
            return None
        return bytes(row[0]), row[1] or ""


def _load_index() -> FingerprintIndex:
    path = PATHS.fpsim_index.with_suffix(".npz")
    if not path.exists():
        raise SystemExit(
            f"fingerprint index not found at {path}.\n"
            f"Build it with: chemmed-ingest index"
        )
    log.info("loading fingerprint index from %s", path)
    started = time.perf_counter()
    index = FingerprintIndex.load(Path(path))
    log.info("loaded %s in %.1fs", index, time.perf_counter() - started)
    return index


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-7s %(name)s :: %(message)s",
    )

    index = _load_index()
    addr = os.environ.get("CHEMWORKER_BIND", "[::]:50051")

    # The search is NumPy-bound and releases the GIL during popcount, so a
    # modest thread pool does give real concurrency here.
    server = grpc.server(
        futures.ThreadPoolExecutor(max_workers=8),
        options=[("grpc.max_send_message_length", 32 * 1024 * 1024)],
    )
    pb_grpc.add_ChemWorkerServicer_to_server(ChemWorkerServicer(index), server)

    health_servicer = health.HealthServicer()
    health_pb2_grpc.add_HealthServicer_to_server(health_servicer, server)
    health_servicer.set("", health_pb2.HealthCheckResponse.SERVING)

    server.add_insecure_port(addr)
    server.start()
    log.info("chemworker listening on %s", addr)

    def _shutdown(*_: object) -> None:
        log.info("shutting down")
        server.stop(grace=5)

    signal.signal(signal.SIGINT, _shutdown)
    signal.signal(signal.SIGTERM, _shutdown)
    server.wait_for_termination()


if __name__ == "__main__":
    main()
