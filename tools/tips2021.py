# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Build absrbrs/data/tips2021.bin from the TIPS-2021 tables (Gamache et al. 2021).

Source: Zenodo record 10.5281/zenodo.4708099 (MIT licence, Copyright (C) 2021 Robert Gamache),
one `<molecule>_<isotopologue>.QTpy` file per isotopologue: a pickled dict of strings, temperature
(K, 1 K steps from 1 K) to Q. Each file's MD5 is checked against Zenodo's record, its opcodes are
inspected (no globals or reductions allowed), and it is loaded with an unpickler that refuses
all globals. Q is stored as f32 (≤ 6e-8 relative rounding, at or below the tables' 6-9
significant digits). Three tables have a non-positive Q at 1 K (31_2 and 31_3: negative; 57_1:
zero), where ln Q is undefined; leading non-positive values are dropped, so those tables start
at 2 K. Run with `uv run tools/tips2021.py`; the outputs are committed.
"""

from __future__ import annotations

import hashlib
import io
import json
import math
import pathlib
import pickle
import pickletools
import struct
import urllib.request

RECORD = "4708099"
API = f"https://zenodo.org/api/records/{RECORD}"
UA = {"User-Agent": "absrbrs-tools"}
DATA = pathlib.Path(__file__).resolve().parents[1] / "absrbrs" / "data"
MAGIC = b"ABSRBRS-TIPS2021"
FORBIDDEN = {
    "GLOBAL",
    "STACK_GLOBAL",
    "REDUCE",
    "INST",
    "OBJ",
    "NEWOBJ",
    "NEWOBJ_EX",
    "BUILD",
    "EXT1",
    "EXT2",
    "EXT4",
}


class Refuse(pickle.Unpickler):
    def find_class(self, module, name):
        raise pickle.UnpicklingError(f"global {module}.{name} refused")


def fetch(url: str) -> bytes:
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA)) as r:
        return r.read()


def load_table(data: bytes) -> list[tuple[float, float]]:
    ops = {op.name for op, _, _ in pickletools.genops(data)}
    if ops & FORBIDDEN:
        raise ValueError(f"unexpected pickle opcodes {sorted(ops & FORBIDDEN)}")
    table = Refuse(io.BytesIO(data)).load()
    if not isinstance(table, dict) or not all(
        isinstance(k, str) and isinstance(v, str) for k, v in table.items()
    ):
        raise ValueError("expected a dict of str to str")
    return sorted((float(k), float(v)) for k, v in table.items())


def main() -> None:
    record = json.loads(fetch(API))
    files = {
        f["key"]: f
        for f in record["files"]
        if f["key"].endswith(".QTpy") and f["key"] != "QT_HITRAN.QTpy"
    }
    entries, values, sources = [], [], {}
    for key in sorted(files, key=lambda k: tuple(int(p) for p in k[:-5].split("_"))):
        meta = files[key]
        data = fetch(meta["links"]["self"])
        md5 = hashlib.md5(data).hexdigest()  # Zenodo publishes MD5
        if f"md5:{md5}" != meta["checksum"]:
            raise ValueError(f"{key}: MD5 {md5} does not match Zenodo's {meta['checksum']}")
        rows = load_table(data)
        temps = [t for t, _ in rows]
        if temps != [1.0 + i for i in range(len(temps))]:
            raise ValueError(f"{key}: temperatures are not 1, 2, …, N K")
        while rows and rows[0][1] <= 0.0:
            rows.pop(0)
        if len(rows) < 2 or not all(0.0 < q < math.inf for _, q in rows):
            raise ValueError(f"{key}: Q is not finite and positive after the leading values")
        mol, iso = (int(p) for p in key[:-5].split("_"))
        entries.append((mol, iso, len(rows), rows[0][0], 1.0))
        values.extend(q for _, q in rows)
        sources[key] = {
            "md5": md5,
            "sha256": hashlib.sha256(data).hexdigest(),
            "points": len(rows),
            "t_first": rows[0][0],
        }
    blob = bytearray(MAGIC)
    blob += struct.pack("<I", len(entries))
    for mol, iso, n, t0, dt in entries:
        blob += struct.pack("<HHIff", mol, iso, n, t0, dt)
    blob += struct.pack(f"<{len(values)}f", *values)
    DATA.mkdir(parents=True, exist_ok=True)
    (DATA / "tips2021.bin").write_bytes(bytes(blob))
    manifest = {
        "source": f"https://doi.org/10.5281/zenodo.{RECORD}",
        "licence": "MIT, Copyright (C) 2021 Robert Gamache"
        " (see absrbrs/LICENSES/TIPS-2021-MIT.txt)",
        "isotopologues": len(entries),
        "points": len(values),
        "sha256": hashlib.sha256(bytes(blob)).hexdigest(),
        "files": sources,
    }
    (DATA / "tips2021.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"{len(entries)} isotopologues, {len(values)} points, {len(blob)} bytes")


if __name__ == "__main__":
    main()
