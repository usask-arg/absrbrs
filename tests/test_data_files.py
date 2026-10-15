"""The embedded TIPS-2021 file matches its manifest."""

import hashlib
import json
import struct
from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "absrbrs" / "data"


def test_tips2021_blob_matches_its_manifest():
    blob = (DATA / "tips2021.bin").read_bytes()
    manifest = json.loads((DATA / "tips2021.json").read_text("utf-8"))
    assert hashlib.sha256(blob).hexdigest() == manifest["sha256"]
    assert blob[:16] == b"ABSRBRS-TIPS2021"
    (count,) = struct.unpack_from("<I", blob, 16)
    assert count == manifest["isotopologues"] == 181
    points = sum(struct.unpack_from("<HHIff", blob, 20 + 16 * k)[2] for k in range(count))
    assert points == manifest["points"]
    assert len(blob) == 20 + 16 * count + 4 * points
