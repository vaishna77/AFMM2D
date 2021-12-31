#ifndef _AFMM2D_HPP__
#define _AFMM2D_HPP__

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "definitions.hpp"
#include "kernel.hpp"
#include "ACA.hpp"
#include "AFMM2DBox.hpp"
#include "AFMM2DTree.hpp"
#include "KDTree.hpp"

class AFMM {
public:
	int N;
	int MinParticlesInLeaf;
	int TOL_POW;
	unsigned nLevels;
	unsigned n_Dimension;  //      Dimension.
	unsigned n_Properties;  //      Number of properties.
	double* locations;   //      Stores all the locations.
	double* properties;  //      Stores all the properties.
	double* sorted_Locations;
	double* sorted_Properties;
	userkernel* mykernel;
	FMM2DTree* afmm2dtree;
  AFMM(int N, int MinParticlesInLeaf, int TOL_POW, Eigen::MatrixXd& loc);
  void assemble();
  Eigen::VectorXd computeMatVecProduct(Eigen::VectorXd inputVecUnsorted);
  void evaluateError();
  ~AFMM();
};

#endif
