"""Deterministic archives: the same members give the same bytes (sorted names, no timestamps or owners)."""

import gzip
import hashlib
import io
import tarfile
from pathlib import Path


def write_tar(path: Path, members: dict[str, bytes]) -> None:
    """``.tar.gz`` of ``members`` (name → bytes), sorted, mtime 0, mode 644, no owner."""
    with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz, \
            tarfile.open(fileobj=gz, mode="w", format=tarfile.PAX_FORMAT) as tf:
        for name in sorted(members):
            data = members[name]
            info = tarfile.TarInfo(name)
            info.size, info.mtime, info.mode = len(data), 0, 0o644
            tf.addfile(info, io.BytesIO(data))


def read_tar(path: Path, suffix: str) -> dict[str, bytes]:
    """Member name → bytes of the regular files ending in ``suffix``."""
    with tarfile.open(path) as tf:
        return {m.name: tf.extractfile(m).read() for m in tf  # type: ignore[union-attr]
                if m.isfile() and m.name.endswith(suffix)}


def sha256(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()
