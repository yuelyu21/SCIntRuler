## Resubmission

This is a resubmission following CRAN's manual review of SCIntRuler 0.99.7.
In response to the review, this update:

* ensures that the package name in the Description field is written in single
  quotes; and
* replaces commented-out code in documentation examples with executable toy
  examples that are run during package checks.

SCIntRuler 0.99.6 was archived on 2025-07-23 after its vignette failed to
rebuild with a newer cowplot release. Version 0.99.7 corrected the legend
extraction that caused the archived check error.

The update also:

* replaces the defunct SeuratObject 5 `slot` interface with version-compatible
  count-layer access; and
* fixes and tests the orientation used when `firstn` selects nearest-neighbor
  distances, and ensures `FindNNDist()` processes every sampled cell; and
* replaces the deprecated `batchelor::cosineNorm()` dependency path with the
  equivalent direct column-wise L2 normalization.

## Test environments

* local macOS arm64, R 4.3.2, Seurat 4.4.0, cowplot 1.1.2
* local macOS arm64, R 4.3.2, Seurat 5.2.1, SeuratObject 5.1.0,
  cowplot 1.2.0

## R CMD check results

0 errors | 0 warnings | 1 note

The note records that this is a new submission of a package archived on
2025-07-23 because issues were not corrected in time. Those archived check
issues are addressed in this release.
