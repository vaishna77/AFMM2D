#ifndef _AFMM2DBox_HPP__
#define _AFMM2DBox_HPP__

#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include "definitions.hpp"

class FMM2DBox {
public:
	int boxNumber;
	int parentNumber;
	int childrenNumbers[4];
	int neighborNumbers[8];
	int innerNumbers[16];
	int outerNumbers[24];

	FMM2DBox ();

	pts2D center;
	Eigen::VectorXd charges;

	std::map<int, Eigen::MatrixXd> M2L;
	Eigen::MatrixXd L2P;					//	Transfer from multipoles of 4 children to multipoles of parent.
	//	The following will be stored only at the leaf nodes
	std::vector<pts2D> chebNodes;
	std::vector<pts2D> particle_loc;
  std::vector<int> chargeLocations;

	Eigen::VectorXd outgoing_charges;//equivalent densities {f_{k}^{B,o}}
	Eigen::VectorXd incoming_charges;//equivalent densities {f_{k}^{B,i}}
	Eigen::VectorXd incoming_potential;//check potentials {u_{k}^{B,i}}
	std::vector<int> incoming_chargePoints;//equivalent points {y_{k}^{B,i}}
	std::vector<int> incoming_checkPoints;//check points {x_{k}^{B,i}}
	Eigen::VectorXd potential;//used only at leaf level
};

#endif
