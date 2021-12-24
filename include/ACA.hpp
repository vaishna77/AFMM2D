#ifndef _ACA_HPP__
#define _ACA_HPP__

#include <bits/stdc++.h>
#include <iostream>
#include <Eigen/Dense>

class LowRank {
public:
	userkernel* K;
	double tol_ACA;
	std::vector<int> row_indices;
	std::vector<int> col_indices;

	LowRank(userkernel* K, int tol_pow, std::vector<int>& row_indices, std::vector<int>& col_indices);
	void maxAbsVector(const Eigen::VectorXd& v, const std::set<int>& allowed_indices,
																double max, int& index
															);
	void ACA_only_nodes(std::vector<int>& row_bases, std::vector<int>& col_bases, int &computed_rank, Eigen::MatrixXd &Ac, Eigen::MatrixXd &Ar);
};

#endif
