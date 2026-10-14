# Third-party material

Material that absrbrs uses, ports or redistributes, with its licence and what that allows.
Checked as part of phase 0 (spec §13); update this file before committing any new third-party
code or data.

Status values: **cleared** (use as described, no further action), **conditions** (use as
described, provided the obligations are met), **blocked** (do not commit; the fallback applies).
"Redistribution" means shipping the item, or a port of it, in the MIT-licensed absrbrs crate,
wheel or repository, for use by anyone, including commercial users. A university's own
operational use counts as commercial use for licence purposes.

Each row records the maintainer who confirmed it, and the date. A new or changed row is pending
until a maintainer confirms it.

| Item | Use in absrbrs | Licence / terms (SPDX) | Classification | Redistribution and use by others | Status | Obligations | Determined by | Checked |
|---|---|---|---|---|---|---|---|---|
| Faddeeva Package, S. G. Johnson (`Faddeeva.cc`, revision of 12 May 2015) | Rust port (Accurate tier, phase 1a) | `MIT`, Copyright (c) 2012 Massachusetts Institute of Technology | permissive | Permitted | **conditions** | Keep the copyright notice and the full MIT permission notice in the header of the ported Rust file, and list it here. Wheels contain no Rust source, and `pyproject.toml` ships only `LICENSE`. So in the same commit as the port, add a notice file inside the crate package (e.g. `absrbrs/LICENSES/Faddeeva-MIT.txt`) and list it in `pyproject.toml` `license-files`. Port from `Faddeeva.cc` only. The authors' code for ACM TOMS Algorithm 916 (Zaghloul and Ali) and Algorithm 680 (Poppe and Wijers) is under restrictive ACM copyright terms, as Johnson's `Faddeeva.cc` header and the Faddeeva Package page both say. Johnson re-implemented both algorithms from the papers alone; do not consult or copy the ACM code. Ported files: `absrbrs/src/faddeeva/accurate/{erfcx,w_im,w}.rs` (each with the notice in its header); notice file `absrbrs/LICENSES/Faddeeva-MIT.txt`, listed in `pyproject.toml` `license-files`. | Daniel Zawada, 2026-10-14 | 2026-10-14 |
| HAPI (HITRAN Application Programming Interface) | Optional, user-installed dependency of the Python helper that builds partition tables, if that helper is built; cross-code comparisons in development | `MIT`, Copyright 2018 HITRAN team. GitHub reports `NOASSERTION`, probably because of an extra title line; the text is the MIT licence. | permissive | Permitted | **cleared** | None while HAPI is imported, not copied. No HAPI code is in the repository: the prototype's "hum1wei24 (HAPI)" baseline is assembled from absrbrs's own methods. If HAPI code is ever copied, keep its copyright and MIT notice with the copy and change this row to **conditions**. | Daniel Zawada, 2026-10-14 | 2026-10-14 |
| TIPS-2021 tables, Gamache et al. (2021) (Zenodo record 10.5281/zenodo.4708099) | Built-in partition functions (phase 1b) | `MIT`: the licence of the Zenodo record, which holds the tables and code. The paper (section 4) also says the codes all use the MIT licence. `TIPS_2021.py` and `BD_TIPS_2021.FOR` name the holder: Copyright (c) 2021 Robert Gamache. | permissive | Permitted | **conditions** | Ship the copyright notice and the MIT permission notice with the embedded tables. In the same commit as the data, add a notice file inside the crate package (e.g. `absrbrs/LICENSES/TIPS-2021-MIT.txt`) and list it in `pyproject.toml` `license-files`, so the notice reaches the crate and the wheel. Cite Gamache et al. (2021), doi:10.1016/j.jqsrt.2021.107713. Build the tables offline, with a script, from the `.QTpy` files. They carry 6–8 significant figures, against 5 in the Fortran tables. They are third-party pickles: inspect each with `pickletools`, then load it with a restricted `pickle.Unpickler` whose `find_class` refuses all globals. Never unpickle at runtime. Use the Fortran DATA tables in `TIPS_2021.zip` (`BD_TIPS_2021.FOR`) as an independent cross-check, which must agree to their 5 significant figures. Record the record DOI, Zenodo's published MD5 and a locally computed SHA-256 for each file used. | Daniel Zawada, 2026-10-14 | 2026-10-14 |
| HITRAN line data (`.par` records) | Small snippets as loader golden tests (phase 1b) | `NONE`. No licence or terms of use stated. HITRAN is free of charge and requires citation; redistribution is not addressed. | no-license | Not permitted (no licence means all rights reserved) | **blocked** | Blocked means: do not commit, ship or publish HITRAN records, including in CI artifacts and caches. It does not cover users, or the spec §8.5 benchmark jobs, downloading HITRAN data at runtime. Use synthetic `.par` records in the loader golden tests. Users supply their own HITRAN files and cite HITRAN themselves. To unblock: an explicit licence or written permission from the HITRAN team. | Daniel Zawada, 2026-10-14 | 2026-10-14 |
| sasktran2 8-pole J-pole constants | Prototype baseline only (`prototypes/faddeeva/faddeeva_prototype/candidates.py`), which is in the repository but not in the crate or the wheel | `MIT`, Copyright (c) 2024 University of Saskatchewan Atmospheric Research Group | permissive | Permitted | **conditions** | The comment above `_JPOLE_B` in `candidates.py` keeps the provenance and the copyright notice (sasktran2 path, copyright holder, MIT licence); keep it with the constants. If the constants ever ship in the crate or the wheel, include sasktran2's copyright notice and MIT permission notice with them. | Daniel Zawada, 2026-10-14 | 2026-10-14 |
| Weideman (1994), "Computation of the complex error function", SIAM J. Numer. Anal. 31(5), 1497–1518, doi:10.1137/0731077 | Algorithm in the prototype and in phase 1a (Fast tier) | Not applicable. `weideman_coefficients` follows the published algorithm, whose few-line listing in the paper is the direct expression of the formula; low risk. | not applicable (published method) | Permitted (low risk) | **cleared** | Cite the paper (scholarly credit, not a licence term). Phase 1a derives the Rust from the paper's equations, not from the listing. | Daniel Zawada, 2026-10-14 | 2026-10-14 |
| python-flint, with FLINT (and the GMP and MPFR libraries bundled in its wheels) | Development tool only: generates the reference data in `reference/`. Not distributed. | python-flint: `MIT`, Copyright (c) 2014-2016 Fredrik Johansson. FLINT: `LGPL-3.0-or-later`. GMP and MPFR: LGPL family (not checked individually). | permissive (python-flint); copyleft (FLINT, GMP, MPFR) | Not redistributed; nothing to permit | **cleared** | None, provided absrbrs never ships or links python-flint, FLINT, GMP or MPFR. The reference data is absrbrs's own output, not a copy of these libraries. | Daniel Zawada, 2026-10-14 | 2026-10-14 |

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
  its `read_me.txt`, the header of `TIPS_2021.py`, which reads "COPYRIGHT (C) [2021]
  [ROBERT GAMACHE]", and `TIPS_2021.zip` (Zenodo MD5 `b8e5baaf59d17d6713b74deb6fba7ae5`, which
  matched the download), whose plain-text `BD_TIPS_2021.FOR` holds the tables under the same
  MIT header, to 5 significant figures (H2O isotopologue 1 at 10 K: `1.3284E+00`). `1_1.QTpy`
  (Zenodo MD5 `6c7100b120a61bc0889179737735ea03`, matched) was inspected with `pickletools`
  without being executed. It is a plain dict of str to str values, with 6–8 significant
  figures (`'1.32837'` at 10 K), and contains no `GLOBAL` or `REDUCE` opcodes. The paper, doi:10.1016/j.jqsrt.2021.107713, is a subscription Elsevier article, and
  ScienceDirect refused automated access. Its copy in the University of Namur repository
  (<https://pure.unamur.be/ws/files/63761906/1_s2.0_S0022407321002065_main.pdf>, found through
  OpenAlex) says in section 4: "The codes all use the MIT license (see codes for details)". It
  names Zenodo 10.5281/zenodo.4708099 as one download location.
- HITRAN: <https://hitran.org/>, <https://hitran.org/faq/>, <https://hitran.org/citepolicy/>
  and <https://hitran.org/register>. None of these pages states a licence, terms of use or
  redistribution terms. The FAQ says "The HITRAN database and its associated compilations are
  free" (free of charge, which is not a licence).
- sasktran2: `~/dev/sasktran2/license.md` (a local file; sasktran2 commit 70d6bb8a). "MIT
  License … Copyright (c) 2024 University of Saskatchewan Atmospheric Research Group". The
  constants are in `rust/sasktran2-rs/src/math/errorfunctions/optimized.rs`, and the repository
  copy carries the notice in a comment in `candidates.py`. sasktran2 does not
  say where the constants originally come from; this was not determined.
- Weideman (1994): Crossref metadata for doi:10.1137/0731077 confirms the citation. The
  implementation is `weideman_coefficients` in `candidates.py`, which computes the coefficients
  by FFT. It follows the paper's few-line listing step by step, because that listing is the
  formula written out.
- python-flint and FLINT: <https://github.com/flintlib/python-flint> (`LICENSE`: "MIT License",
  Copyright (c) 2014-2016 Fredrik Johansson) and <https://github.com/flintlib/flint> (README:
  "LGPL (GNU Lesser General Public License) version 3 or later"). python-flint's README still
  says "LGPL v2.1+" for FLINT; FLINT 3's own README takes precedence. For a development tool
  this makes no difference.

## Related, not used

- TIPS-2025 (Gamache et al. 2025, doi:10.1016/j.jqsrt.2025.109568; Zenodo record
  10.5281/zenodo.17191976, for HITRAN2024) is licensed `CC-BY-4.0` according to its Zenodo
  record. HAPI 1.3 embeds it, and its `partitionSum` uses it by default. absrbrs does not use
  it; check it properly before adopting it.
