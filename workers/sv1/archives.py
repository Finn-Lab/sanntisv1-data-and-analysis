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


def limit_pair(gbk: Path, json_: Path, out: Path, n: int) -> tuple[Path, Path]:
    """A trial subset of a MIBiG-format archive pair (LIMIT runs): ``n`` BGCs evenly spaced over the sorted
    accessions, so both kingdoms and the whole accession range stay represented. Written to ``out`` under the
    same file names."""
    g, j = read_tar(gbk, ".gbk"), read_tar(json_, ".json")
    accs = sorted(Path(m).stem for m in g)
    step = max(len(accs) // max(n, 1), 1)
    keep = set(accs[::step][:n])
    out.mkdir(parents=True, exist_ok=True)
    paths = (out / Path(gbk).name, out / Path(json_).name)
    write_tar(paths[0], {m: b for m, b in g.items() if Path(m).stem in keep})
    write_tar(paths[1], {m: b for m, b in j.items() if Path(m).stem in keep})
    return paths


if __name__ == "__main__":  # python -m sv1.archives limit <gbk> <json> <out dir> <n>  → prints the two paths
    import sys

    if sys.argv[1:2] != ["limit"] or len(sys.argv) != 6:
        sys.exit("usage: python -m sv1.archives limit <gbk.tar.gz> <json.tar.gz> <out dir> <n>")
    print("\n".join(str(p) for p in limit_pair(Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]),
                                               int(sys.argv[5]))))
