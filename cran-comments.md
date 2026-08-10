## Resubmission

SCIntRuler 0.99.6 was archived on 2025-07-23 after its vignette failed to
rebuild with a newer cowplot release. This update corrects the legend
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

The note was environmental: `unable to verify current time`.
