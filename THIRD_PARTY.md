# Third-party material

Material that absrbrs uses, ports or redistributes, with its licence and what that allows.
Checked as part of phase 0 (spec §13); update this file before committing any new third-party
code or data.

Status values: **cleared** (use as described, no further action), **conditions** (use as
described, provided the obligations are met), **blocked** (do not commit; the fallback applies).
"Redistribution" means shipping the item, or a port of it, in the MIT-licensed absrbrs crate,
wheel or repository, for use by anyone, including commercial users. A university's own
operational use counts as commercial use for licence purposes.

Every row is pending confirmation by a named maintainer. When someone confirms a row, replace
"pending maintainer confirmation" with their name and the date.

| Item | Use in absrbrs | Licence / terms (SPDX) | Classification | Redistribution and use by others | Status | Obligations | Determined by | Checked |
|---|---|---|---|---|---|---|---|---|
| Faddeeva Package, S. G. Johnson (`Faddeeva.cc`, revision of 12 May 2015) | Rust port (Accurate tier, phase 1a) | `MIT`, Copyright (c) 2012 Massachusetts Institute of Technology | permissive | Permitted | **conditions** | Keep the copyright notice and the full MIT permission notice in the header of the ported Rust file, and list it here. Port from `Faddeeva.cc` only. ACM TOMS Algorithm 916 and the Poppe–Wijers code are under ACM's restrictive terms (spec §6.3), and Johnson re-implemented them from the papers alone; do not consult or copy them. | pending maintainer confirmation | 2026-10-14 |
| HAPI (HITRAN Application Programming Interface) | Optional, user-installed dependency of the Python helper that builds partition tables; cross-code comparisons in development | `MIT`, Copyright 2018 HITRAN team. GitHub reports `NOASSERTION` because of an extra title line; the text is the MIT licence. | permissive | Permitted | **cleared** | None while HAPI is imported, not copied. No HAPI code is in the repository: the prototype's "hum1wei24 (HAPI)" baseline is assembled from absrbrs's own methods. If HAPI code is ever copied, keep its copyright and MIT notice with the copy and change this row to **conditions**. | pending maintainer confirmation | 2026-10-14 |
| TIPS-2021 tables, Gamache et al. (2021) (Zenodo record 10.5281/zenodo.4708099) | Built-in partition functions (phase 1b) | `MIT`: the licence of the Zenodo record, which holds the tables and code. `TIPS_2021.py` names the holder: Copyright (c) 2021 Robert Gamache. | permissive | Permitted | **conditions** | Ship the copyright notice and the MIT permission notice with the embedded tables, in the crate and in the wheel (a notice file beside the data, or the header of the generated source). Cite Gamache et al. (2021), doi:10.1016/j.jqsrt.2021.107713. Build the tables with a script from the Zenodo files, and record the record DOI and the file hashes. | pending maintainer confirmation | 2026-10-14 |
| HITRAN line data (`.par` records) | Small snippets as loader golden tests (phase 1b) | `NONE`. No licence or terms of use stated. HITRAN is free of charge and requires citation; redistribution is not addressed. | no-license | Not permitted (no licence means all rights reserved) | **blocked** | Do not commit HITRAN records. Use synthetic `.par` records in the loader golden tests. Users supply their own HITRAN files and cite HITRAN themselves. To unblock: an explicit licence or written permission from the HITRAN team. | pending maintainer confirmation | 2026-10-14 |
| sasktran2 8-pole J-pole constants | Prototype baseline only (`prototypes/faddeeva/faddeeva_prototype/candidates.py`); not in the crate or the wheel | `MIT`, Copyright (c) 2024 University of Saskatchewan Atmospheric Research Group | permissive | Permitted | **conditions** | Keep the provenance comment in `candidates.py`. If the constants ever ship in the crate or the wheel, include sasktran2's copyright notice and MIT permission notice with them. | pending maintainer confirmation | 2026-10-14 |
| Weideman (1994), "Computation of the complex error function", SIAM J. Numer. Anal. 31(5), 1497–1518, doi:10.1137/0731077 | Algorithm in the prototype and in phase 1a (Fast tier) | Not applicable. This is a published method, implemented independently; no code or text is copied, and the coefficients are computed by absrbrs's own code. | not applicable (method, not a copied work) | Permitted (nothing is copied) | **cleared** | Cite the paper (scholarly credit, not a licence term). Any code printed in the paper is under the publisher's copyright: do not copy it. | pending maintainer confirmation | 2026-10-14 |
| python-flint, with FLINT (and the GMP and MPFR libraries bundled in its wheels) | Development tool only: generates the reference data in `reference/`. Not distributed. | python-flint: `MIT`, Copyright (c) 2014-2016 Fredrik Johansson. FLINT: `LGPL-3.0-or-later`. GMP and MPFR: LGPL family (not checked individually). | permissive (python-flint); copyleft (FLINT, GMP, MPFR) | Not redistributed; nothing to permit | **cleared** | None, provided absrbrs never ships or links python-flint, FLINT, GMP or MPFR. The reference data is absrbrs's own output, not a copy of these libraries. | pending maintainer confirmation | 2026-10-14 |

## Sources checked

All checked on 2026-10-14.

- Faddeeva Package: <http://ab-initio.mit.edu/Faddeeva> (served `Faddeeva.c`, which includes
  `Faddeeva.cc`); <http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package> (the
  documentation page, with a licence section); and the licence header of
  <http://ab-initio.mit.edu/Faddeeva.cc> (latest revision 12 May 2015; SHA-256
  `46dc0d4e7d81ddcce5a95953b97d31696ffb424228e751aac10dfd9aeb3430a0`). The page says the package
  is distributed "under the 'MIT License' (also called the Expat License)". The host refused
  HTTPS on this date, so the pages were fetched over plain HTTP.
- HAPI: <https://github.com/hitranonline/hapi> and
  <https://github.com/hitranonline/hapi/blob/master/LICENSE.txt>. "MIT License … Copyright
  2018 HITRAN team". The source, `hapi/hapi.py`, embeds the TIPS-2021 tables and credits Gamache
  et al. (2021) and the Zenodo record below.
- TIPS-2021: <https://zenodo.org/records/4708099> (the distribution that HAPI names; the record's
  licence is "mit-license", and it is the only version of concept DOI 10.5281/zenodo.4708098),
  its `read_me.txt` and the header of `TIPS_2021.py`, which reads "COPYRIGHT (C) [2021]
  [ROBERT GAMACHE]". Also checked <https://doi.org/10.1016/j.jqsrt.2021.107713> (via Crossref
  metadata): a subscription Elsevier article. The data-availability statement could not be read,
  because ScienceDirect refused automated access. The licence above therefore rests on the
  Zenodo record, not on the paper.
- HITRAN: <https://hitran.org/>, <https://hitran.org/faq/>, <https://hitran.org/citepolicy/>
  and <https://hitran.org/register>. None of these pages states a licence, terms of use or
  redistribution terms. The FAQ says "The HITRAN database and its associated compilations are
  free" (free of charge, which is not a licence).
- sasktran2: `~/dev/sasktran2/license.md` (a local file; sasktran2 commit 70d6bb8a). "MIT
  License … Copyright (c) 2024 University of Saskatchewan Atmospheric Research Group". The
  constants are in `rust/sasktran2-rs/src/math/errorfunctions/optimized.rs`. sasktran2 does not
  say where the constants originally come from; this was not determined.
- Weideman (1994): Crossref metadata for doi:10.1137/0731077 confirms the citation. The
  implementation is `weideman_coefficients` in `candidates.py`, which computes the coefficients
  by FFT.
- python-flint and FLINT: <https://github.com/flintlib/python-flint> (`LICENSE`: "MIT License",
  Copyright (c) 2014-2016 Fredrik Johansson) and <https://github.com/flintlib/flint> (README:
  "LGPL (GNU Lesser General Public License) version 3 or later"). python-flint's README still
  says "LGPL v2.1+" for FLINT; FLINT 3's own README takes precedence. For a development tool
  this makes no difference.

## Related, not used

- TIPS-2025 (Gamache et al. 2025, doi:10.1016/j.jqsrt.2025.109568; Zenodo record
  10.5281/zenodo.17191976, for HITRAN2024) is licensed `CC-BY-4.0` according to its Zenodo
  record. HAPI 1.3 embeds it. absrbrs does not use it; check it properly before adopting it.
