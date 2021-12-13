#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cstdlib>
#define EIGEN_DONT_PARALLELIZE

const double PI	=	3.1415926535897932384;
#include <map>
struct pts2D {
	double x,y;
};

class kernel {
public:
  double a;
  std::vector<pts2D> particles_X;
	std::vector<pts2D> particles_Y;

	kernel(std::vector<pts2D>& particles_X, std::vector<pts2D>& particles_Y) {
			this->particles_X = particles_X;
			this->particles_Y = particles_Y;
	}

	virtual double getMatrixEntry(const unsigned i, const unsigned j) {
		std::cout << "virtual getInteraction" << std::endl;
		return 0.0;
	}

	Eigen::VectorXd getRow(const int j, std::vector<int> col_indices) {
		int n_cols = col_indices.size();
		Eigen::VectorXd row(n_cols);
    #pragma omp parallel for
    for(int k = 0; k < n_cols; k++) {
        row(k) = this->getMatrixEntry(j, col_indices[k]);
    }
    return row;
  }

  Eigen::VectorXd getCol(const int k, std::vector<int> row_indices) {
		int n_rows = row_indices.size();
    Eigen::VectorXd col(n_rows);
    #pragma omp parallel for
    for (int j=0; j<n_rows; ++j) {
			col(j) = this->getMatrixEntry(row_indices[j], k);
    }
    return col;
  }

  Eigen::MatrixXd getMatrix(std::vector<int> row_indices, std::vector<int> col_indices) {
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
  ~kernel() {};
};

class userkernel: public kernel {
public:
	double chargesFunction(const pts2D r) {
		double q = r.x; //user defined
		return q;
	};
	// #ifdef ONEOVERR
	// userkernel(std::vector<pts2D>& particles_X, std::vector<pts2D>& particles_Y): kernel(particles_X, particles_Y) {
	// };
	// double getMatrixEntry(const unsigned i, const unsigned j) {
	// 	pts2D r1 = particles_X[i];
	// 	pts2D r2 = particles_X[j];
	// 	double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
	// 	double R	=	sqrt(R2);
	// 	if (R < a) {
	// 		return R/a;
	// 	}
	// 	else {
	// 		return a/R;
	// 	}
	// }
	// #elif LOGR
	userkernel(std::vector<pts2D> particles_X, std::vector<pts2D> particles_Y): kernel(particles_X, particles_Y) {
	};
	double getMatrixEntry(const unsigned i, const unsigned j) {
		pts2D r1 = particles_X[i];
		pts2D r2 = particles_X[j];
		double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
		if (R2 < 1e-10) {
			return 0.0;
		}
		else if (R2 < a*a) {
			return 0.5*R2*log(R2)/a/a;
		}
		else {
			return 0.5*log(R2);
		}
	}
	// #endif
	~userkernel() {};
};
