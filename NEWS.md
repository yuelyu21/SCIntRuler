# SCIntRuler 0.99.7

- Fixed vignette rendering with current `cowplot` releases.
- Updated count-layer access for SeuratObject 5 while retaining compatibility
  with SeuratObject 4.
- Corrected `PermTest()` and `FindCell()` so `firstn` selects nearest-neighbor
  rows, as defined in the published method, rather than sampled-cell columns.
- Corrected `FindNNDist()` to calculate distances for every sampled cell.
- Added validation and regression tests for `firstn`.
