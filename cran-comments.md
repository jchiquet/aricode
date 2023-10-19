This submission follows the reporting of a bug in the calculation of the AMI.

## aricode 1.0.3

- minor fix in package documentation due to evolution of roxygen2 7.0.0 <https://github.com/r-lib/roxygen2/issues/1491>.

## Tested environments

* tested locally on Ubuntu Linux 22.04.1 LTS, R-release, GCC

* tested remotely with github-action

- Linux ubuntu 22.04, R-release (github-action)
- Linux ubuntu 22.04, R-oldrel (github-action)
- Linux ubuntu 22.04, R-devel (github-action)
- Windows Server 2022, R-release, 64 bit
- macOS 12, R-release (github action)
- Linux ubuntu 22.04, R-release, gcc + unit test with sanitizers (github-action)
- Linux ubuntu 22.04, R-release, clang + unit test with sanitizers (github-action)

* tested remotely with win-builder (R version 4.3.1, R unstable, R version 4.2.3)

all status OK

## Local R CMD check results

── R CMD check results ──────────────────────────────────────────────── aricode 1.0.3 ────
Duration: 43.1s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

R CMD check succeeded
