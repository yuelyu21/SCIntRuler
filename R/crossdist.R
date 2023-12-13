#' Calculate Cross-Distance Matrix
#'
#' This function computes the pairwise Euclidean distance between rows of two matrices.
#' @name crossdist
#' @param m1 A numeric matrix.
#' @param m2 Another numeric matrix.
#' @return A numeric matrix where the entry at `[i, j]` is the Euclidean distance between row i of `m1` and row j of `m2`.
#' @export
#' @examples
#' \dontrun{
#' mat1 <- matrix(1:4, ncol = 2)
#' mat2 <- matrix(5:8, ncol = 2)
#' crossdist(mat1, mat2)
#'}
Rcpp::cppFunction('NumericMatrix crossdist(NumericMatrix m1, NumericMatrix m2) {
    int nrow1 = m1.nrow();
    int nrow2 = m2.nrow();
    int ncol = m1.ncol();

    if (ncol != m2.ncol()) {
      throw std::runtime_error("Incompatible number of dimensions");
    }

    NumericMatrix out(nrow1, nrow2);

    for (int r1 = 0; r1 < nrow1; r1++) {
      for (int r2 = 0; r2 < nrow2; r2++) {
        double total = 0;
        for (int c12 = 0; c12 < ncol; c12++) {
          total += pow(m1(r1, c12) - m2(r2, c12), 2);
        }
        out(r1, r2) = sqrt(total);
      }
    }

    return out;
  }')
