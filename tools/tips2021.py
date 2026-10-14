# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Build absrbrs/data/tips2021.bin from the TIPS-2021 tables (Gamache et al. 2021).

Source: Zenodo record 10.5281/zenodo.4708099 (MIT licence, Copyright (C) 2021 Robert Gamache),
one `<molecule>_<isotopologue>.QTpy` file per isotopologue: a pickled dict of strings, temperature
(K, 1 K steps from 1 K) to Q. Each file's MD5 is checked against Zenodo's record, its opcodes are
checked against an allowlist (a plain dict of strings: no globals, reductions or object
construction), and it is loaded with an unpickler that refuses all globals. Q is stored as f32
(≤ 6e-8 relative rounding, at or below the tables' 6-9 significant digits).

Low temperatures: three tables have a non-positive Q at 1 K (31_2 and 31_3: negative; 57_1:
zero), where ln Q is undefined; leading non-positive values are dropped. The H2S tables 31_2 and
31_3 are also unphysical up to about 30 K, so they start at 30 K (`MIN_TEMPERATURE`). 57_1
starts at 2 K, but its low-temperature values carry few digits (Q(2 K) = 0.000416, about 1e-3
relative rounding).

Cross-check: the Fortran DATA tables of `BD_TIPS_2021.FOR` (in `TIPS_2021.zip`, same record,
MD5-checked) give Q to 5 significant figures on 1-20 K by 1 K, then 22-9000 K by 2 K. Every stored
value at a shared node must agree with them to within half a unit in the 5th significant figure
(plus the f32 rounding of the stored value). The result is recorded in the manifest.

Run with `uv run tools/tips2021.py`; the outputs are committed.
"""

from __future__ import annotations

import hashlib
import http.client
import io
import json
import math
import pathlib
import pickle
import pickletools
import re
import struct
import time
import urllib.error
import urllib.request
import zipfile

RECORD = "4708099"
API = f"https://zenodo.org/api/records/{RECORD}"
UA = {"User-Agent": "absrbrs-tools"}
DATA = pathlib.Path(__file__).resolve().parents[1] / "absrbrs" / "data"
MAGIC = b"ABSRBRS-TIPS2021"
# The opcodes of a pickled dict of str to str (protocols 2-4); anything else is refused.
ALLOWED = {
    "PROTO",
    "FRAME",
    "EMPTY_DICT",
    "MARK",
    "BINUNICODE",
    "SHORT_BINUNICODE",
    "BINPUT",
    "LONG_BINPUT",
    "MEMOIZE",
    "SETITEM",
    "SETITEMS",
    "STOP",
}
# Tables whose lowest temperatures are unphysical start here (K). For H2S the ratio
# Q(31_2)/Q(31_1) should be about 1 and Q(31_3)/Q(31_1) about 4 (the ratio of the nuclear-spin
# degeneracies), but they are 0.44 and 1.77 at 2 K, 1.43 and 5.72 at 6 K and 1.03 and 4.11 at
# 20 K; from 30 K they are within 0.6% of it.
MIN_TEMPERATURE = {"31_2": 30.0, "31_3": 30.0}
FORTRAN_ZIP = "TIPS_2021.zip"
FORTRAN_FILE = "BD_TIPS_2021.FOR"
# The temperatures of the Fortran QOFT tables (TDAT in BD_TIPS_2021.FOR).
FORTRAN_T = [float(t) for t in range(1, 21)] + [float(t) for t in range(22, 9001, 2)]


class Refuse(pickle.Unpickler):
    def find_class(self, module, name):
        raise pickle.UnpicklingError(f"global {module}.{name} refused")


def fetch(url: str, attempts: int = 4) -> bytes:
    """Downloads `url`, retrying (after 5, 10, 20 s) on server errors and network failures,
    including truncated transfers."""
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code < 500 or attempt == attempts - 1:
                raise
        except (http.client.HTTPException, OSError):  # incl. URLError and timeouts
            if attempt == attempts - 1:
                raise
        time.sleep(5 * 2**attempt)
    raise AssertionError("unreachable")


def fetch_checked(meta: dict) -> tuple[bytes, str]:
    """Downloads a file of the record, checks its MD5 against Zenodo's, returns both."""
    data = fetch(meta["links"]["self"])
    md5 = hashlib.md5(data).hexdigest()  # Zenodo publishes MD5
    if f"md5:{md5}" != meta["checksum"]:
        raise ValueError(f"{meta['key']}: MD5 {md5} does not match Zenodo's {meta['checksum']}")
    return data, md5


def load_table(data: bytes) -> list[tuple[float, float]]:
    ops = {op.name for op, _, _ in pickletools.genops(data)}
    if ops - ALLOWED:
        raise ValueError(f"unexpected pickle opcodes {sorted(ops - ALLOWED)}")
    table = Refuse(io.BytesIO(data)).load()
    if not isinstance(table, dict) or not all(
        isinstance(k, str) and isinstance(v, str) for k, v in table.items()
    ):
        raise ValueError("expected a dict of str to str")
    return sorted((float(k), float(v)) for k, v in table.items())


def fortran_tables(zip_data: bytes) -> dict[tuple[int, int], dict[int, str]]:
    """The QOFT DATA blocks of BD_TIPS_2021.FOR: (molecule, isotopologue) -> {node index: text}.

    Molecules are numbered by the order of the `SUBROUTINE QT_` blocks (QT_H2O is 1, ...).
    """
    with zipfile.ZipFile(io.BytesIO(zip_data)) as z:
        names = [
            n
            for n in z.namelist()
            if n.rsplit("/", 1)[-1] == FORTRAN_FILE and not n.startswith("__MACOSX")
        ]
        if len(names) != 1:
            raise ValueError(f"{FORTRAN_ZIP}: expected one {FORTRAN_FILE}, found {names}")
        lines = z.read(names[0]).decode("latin-1").splitlines()
    tables: dict[tuple[int, int], dict[int, str]] = {}
    molecule, i = 0, 0
    while i < len(lines):
        line = lines[i]
        if re.match(r"\s+SUBROUTINE QT_", line):
            molecule += 1
        m = re.match(r"\s+DATA \(QOFT\(\s*(\d+),J\)\s*,\s*J=\s*(\d+),\s*(\d+)\)/(.*)", line)
        if m:
            iso, first, last, body = int(m[1]), int(m[2]), int(m[3]), m[4]
            while "/" not in body:  # fixed-form continuation lines
                i += 1
                body += " " + lines[i][6:]
            values = re.findall(r"[-+]?\d\.\d+E[-+]\d+", body[: body.index("/")])
            if len(values) != last - first + 1:
                raise ValueError(f"{FORTRAN_FILE}: QOFT({iso}) block {first}-{last} is malformed")
            nodes = tables.setdefault((molecule, iso), {})
            nodes.update(zip(range(first, last + 1), values, strict=True))
        i += 1
    return tables


def cross_check(
    fortran: dict[tuple[int, int], dict[int, str]],
    stored: dict[tuple[int, int], tuple[float, tuple[float, ...]]],
) -> tuple[int, float]:
    """Compares the stored f32 values with the Fortran tables; returns (nodes, worst ratio)."""
    if set(fortran) != set(stored):
        raise ValueError(f"{FORTRAN_FILE}: isotopologues differ ({set(fortran) ^ set(stored)})")
    compared, worst = 0, 0.0
    for key, nodes in sorted(fortran.items()):
        t_first, q = stored[key]
        for j, text in sorted(nodes.items()):
            t = FORTRAN_T[j - 1]
            if t < t_first:  # dropped from the stored table
                continue
            k = int(t - t_first)
            value = float(text)
            if k >= len(q):  # past the end of the .QTpy table: the Fortran pads with zeros
                if value != 0.0:
                    raise ValueError(f"{key} at {t} K: Fortran {text} past the end of the table")
                continue
            # Half a unit in the 5th significant figure of d.dddd E±xx, with a little slack for
            # exact decimal ties, plus the f32 rounding of the stored value.
            half_unit = 0.5 * 10.0 ** (int(text.split("E")[1]) - 4)
            bound = half_unit * (1 + 1e-9) + abs(q[k]) * 2.0**-24
            ratio = abs(q[k] - value) / bound
            if not ratio <= 1.0:
                raise ValueError(f"{key} at {t} K: stored {q[k]} vs Fortran {text}")
            compared += 1
            worst = max(worst, ratio)
    return compared, worst


def main() -> None:
    record = json.loads(fetch(API))
    by_key = {f["key"]: f for f in record["files"]}
    files = {k: f for k, f in by_key.items() if k.endswith(".QTpy") and k != "QT_HITRAN.QTpy"}
    entries, values, sources = [], [], {}
    for key in sorted(files, key=lambda k: tuple(int(p) for p in k[:-5].split("_"))):
        data, md5 = fetch_checked(files[key])
        rows = load_table(data)
        temps = [t for t, _ in rows]
        if temps != [1.0 + i for i in range(len(temps))]:
            raise ValueError(f"{key}: temperatures are not 1, 2, …, N K")
        name = key[:-5]
        rows = [(t, q) for t, q in rows if t >= MIN_TEMPERATURE.get(name, 0.0)]
        while rows and rows[0][1] <= 0.0:
            rows.pop(0)
        if len(rows) < 2 or not all(0.0 < q < math.inf for _, q in rows):
            raise ValueError(f"{key}: Q is not finite and positive after the leading values")
        mol, iso = (int(p) for p in name.split("_"))
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
    packed = struct.pack(f"<{len(values)}f", *values)
    blob += packed

    # Cross-check the stored (f32) values against the Fortran tables.
    zip_data, zip_md5 = fetch_checked(by_key[FORTRAN_ZIP])
    as_f32 = struct.unpack(f"<{len(values)}f", packed)
    stored, offset = {}, 0
    for mol, iso, n, t0, _ in entries:
        stored[(mol, iso)] = (t0, as_f32[offset : offset + n])
        offset += n
    compared, worst = cross_check(fortran_tables(zip_data), stored)

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
        "fortran_check": {
            "file": f"{FORTRAN_ZIP} ({FORTRAN_FILE})",
            "md5": zip_md5,
            "sha256": hashlib.sha256(zip_data).hexdigest(),
            "bound": "half a unit in the 5th significant figure of the Fortran value"
            " plus the f32 rounding of the stored value",
            "nodes_compared": compared,
            "worst_ratio_to_bound": round(worst, 6),
        },
    }
    (DATA / "tips2021.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"{len(entries)} isotopologues, {len(values)} points, {len(blob)} bytes")
    print(f"Fortran cross-check: {compared} nodes, worst ratio to bound {worst:.6f}")


if __name__ == "__main__":
    main()
