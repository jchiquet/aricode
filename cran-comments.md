This submission follows the reporting of a bug in the calculation of the AMI.

## aricode 1.0.2

- fix serious bug in computation of AMI, thanks to Florent Bersani

## Tested environments

* tested locally on Ubuntu Linux 22.04.1 LTS, R-release, GCC

* tested remotely with github-action

- Linux ubuntu 20.04, R-release (github-action)
- Linux ubuntu 20.04, R-oldrel (github-action)
- Linux ubuntu 20.04, R-devel (github-action)
- Windows Server 2022, R-release, 64 bit
- macOS Big Sur 11, R-release (github action)
- Linux ubuntu 20.04, R-release, gcc + unit test with sanitizers (github-action)
- Linux ubuntu 20.04, R-release, clang + unit test with sanitizers (github-action)

* tested remotely with win-builder (R version 4.2.1, R dev, R version 4.1.3)

all status OK 

## Local R CMD check results

── R CMD check results ─────────────────────────────────────────────────────────────────────────────── aricode 1.0.2 ────
Duration: 22.8s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

R CMD check succeeded

