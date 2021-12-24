#ifndef _kernel_HPP__
#define _kernel_HPP__

#include <iostream>
#include <vector>
#include <Eigen/Dense>
// struct pts2D {
// 	double x,y;
// };
#include "definitions.hpp"

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

	Eigen::VectorXd getRow(const int j, std::vector<int> col_indices);
  Eigen::VectorXd getCol(const int k, std::vector<int> row_indices);
  Eigen::MatrixXd getMatrix(int row_start_index, int row_end_index, int col_start_index, int col_end_index);
  Eigen::MatrixXd getMatrix(std::vector<int> row_indices, std::vector<int> col_indices);

  ~kernel() {};
};

class userkernel: public kernel {
public:
	double defineVector(const pts2D r);
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

	double getMatrixEntry(const unsigned i, const unsigned j);
	// #endif
	~userkernel() {};
};

#endif
