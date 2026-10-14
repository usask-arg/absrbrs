| Algorithm | Max value error | Max derivative error | Holdout value error | Holdout derivative error | Expected cost |
|---|---|---|---|---|---|
| jpole8 (sasktran2) | 8.4e+294 | 3.3e-06 | 3.9e+297 | 3.3e-06 | 240 |
| hum1wei24 (HAPI) | 2.6e+86 | 1.7e-04 | 2.4e+86 | 1.9e-04 | 66 |
| Weideman 32 | 1.3e+290 | 9.8e-12 | 5.7e+292 | 9.8e-12 | 568 |
| Weideman 40 | 2.5e+287 | 1.4e-14 | 1.1e+290 | 1.4e-14 | 696 |
| **Selected composite** | 4.8e-07 | 1.8e-06 | 4.8e-07 | 1.8e-06 | 61 |

Regions of the selected composite (first match wins):

| Region | Method | Max value error | Max derivative error | Workload share |
|---|---|---|---|---|
| exp + CF | gauss_cf 5 | 1.4e-07 | 1.4e-07 | 0.28% |
| CF |z| in [8, 10) | cf 4 | 5.1e-08 | 5.1e-08 | 0.05% |
| CF |z| in [10, 12) | cf 3 | 1.5e-07 | 1.5e-07 | 0.05% |
| CF |z| in [12, 15) | cf 3 | 3.4e-08 | 3.4e-08 | 0.11% |
| CF |z| in [15, 20) | cf 2 | 4.8e-07 | 4.8e-07 | 0.15% |
| CF |z| in [20, 30) | cf 2 | 8.3e-08 | 8.3e-08 | 0.33% |
| CF |z| in [30, 50) | cf 2 | 7.3e-09 | 7.3e-09 | 0.85% |
| CF |z| in [50, 100) | cf 1 | 4.0e-07 | 4.0e-07 | 2.11% |
| CF |z| in [100, 300) | cf 1 | 2.5e-08 | 2.5e-08 | 8.43% |
| CF |z| in [300, 1000) | cf 1 | 3.1e-10 | 3.1e-10 | 29.46% |
| CF |z| in [1000, 3000) | cf 1 | 2.5e-12 | 2.5e-12 | 44.64% |
| CF |z| in [3000, inf) | cf 0 | 1.7e-07 | 1.7e-07 | 13.27% |
| Taylor | taylor 3 | 8.7e-08 | 1.8e-06 | 0.08% |
| weideman 24 | weideman 24 | 2.4e-07 | 9.6e-09 | 0.19% |
