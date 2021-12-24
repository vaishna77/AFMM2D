#include "kernel.hpp"
// #include "../include/kernel.hpp"
Eigen::VectorXd kernel::getRow(const int j, std::vector<int> col_indices) {
  int n_cols = col_indices.size();
  Eigen::VectorXd row(n_cols);
  #pragma omp parallel for
  for(int k = 0; k < n_cols; k++) {
      row(k) = this->getMatrixEntry(j, col_indices[k]);
  }
  return row;
}

Eigen::VectorXd kernel::getCol(const int k, std::vector<int> row_indices) {
  int n_rows = row_indices.size();
  Eigen::VectorXd col(n_rows);
  #pragma omp parallel for
  for (int j=0; j<n_rows; ++j) {
    col(j) = this->getMatrixEntry(row_indices[j], k);
  }
  return col;
}

Eigen::MatrixXd kernel::getMatrix(int row_start_index, int col_start_index, int row_end_index, int col_end_index) {
  Eigen::MatrixXd mat(row_end_index-row_start_index, col_end_index-col_start_index);
  #pragma omp parallel for
  for (int j=row_start_index; j < row_end_index; ++j) {
      #pragma omp parallel for
      for (int k=col_start_index; k < col_end_index; ++k) {
          mat(j,k) = this->getMatrixEntry(j, k);
      }
  }
  return mat;
}

Eigen::MatrixXd kernel::getMatrix(std::vector<int> row_indices, std::vector<int> col_indices) {
  int n_rows = row_indices.size();
  int n_cols = col_indices.size();
  Eigen::MatrixXd mat(n_rows, n_cols);
  #pragma omp parallel for
  for (int j=0; j < n_rows; ++j) {
      #pragma omp parallel for
      for (int k=0; k < n_cols; ++k) {
          mat(j,k) = this->getMatrixEntry(row_indices[j], col_indices[k]);
      }
  }
  return mat;
}

double userkernel::defineVector(const pts2D r) {
  double q = r.x; //user defined
  return q;
}
