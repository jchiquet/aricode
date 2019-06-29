
## Package update

This package uses a bucket-sorting algorithm to avoid useless storage and computations in sparse contingency tables when comparing two classifications with many entries. It is then used to compute fastly adjusted Rand index and other clustering comparison measures.

I recently added a new measure of clustering comparison (the adjusted mutual information - AMI)

Thank you for your time and work on CRAN.

## Tested environments

- local ubuntu 18.04 install, R 3.6.0
- ubuntu 14.04 (on travis-ci), R 3.6.0, oldrelease and R devel
- Mac OS X (on travis-ci), R 3.6.0 and oldrelease
- win-builder, R 3.6.0 and R devel
- Debian (R-hub), R devel

## R CMD check results

── R CMD check results ────────────────────────────────────── aricode 0.1.2 ────
Duration: 1m 2.8s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

R CMD check succeeded
