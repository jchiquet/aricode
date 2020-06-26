
## Package update

This package uses a bucket-sorting algorithm to avoid useless storage and computations in sparse contingency tables when comparing two classifications with many entries. It is then used to compute fastly adjusted Rand index and other clustering comparison measures. In this new version we

  -	added the Modified adjusted Rand Index
  -	added the Chi-Square statistics
  - fixed bugs

Thank you for your time and work on CRAN.

## Tested environments

- macOS 10.13.6 High Sierra, R-release, CRAN's setup (R hub)
- macOS Catalina 10.15, R-release (github action)
- macOS Catalina 10.15, R-devel (github action)
- Linux ubuntu 16.04, R-release (github-action)
- Linux ubuntu 18.04 R-release, (local)
- Linux Debian GCC  R-devel, (R hub)
- Windows Server 2019, R-release (github action)
- Windows Server 2008 R2 SP1, R-release  (R hub)
- Windows, R-oldrelease (winbuilder)
- Windows, R-release (winbuilder)
- Windows, R-devel  (winbuilder)

all status OK 

## R CMD check results

── R CMD check results ────────────────────────────────────── aricode 1.0.0 ────
Duration: 26.5s

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

R CMD check succeeded
