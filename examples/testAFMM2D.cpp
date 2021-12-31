#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <map>

#include "AFMM2D.hpp"

double userkernel::getMatrixEntry(const unsigned i, const unsigned j) {
	pts2D r1 = particles_X[i];
	pts2D r2 = particles_X[j];
	double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
	double R	=	sqrt(R2);
	if (R < a) {
		return R/a;
	}
	else {
		return a/R;
	}
}

// double userkernel::getMatrixEntry(const unsigned i, const unsigned j) {
//   pts2D r1 = particles_X[i];
//   pts2D r2 = particles_X[j];
//   double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
//   if (R2 < 1e-10) {
//     return 0.0;
//   }
//   else if (R2 < a*a) {
//     return 0.5*R2*log(R2)/a/a;
//   }
//   else {
//     return 0.5*log(R2);
//   }
// }

int main(int argc, char* argv[]) {
  unsigned N;
  unsigned MinParticlesInLeaf;
  int TOL_POW;
  if(argc < 4)
  {
      std::cout << "All arguments weren't passed to executable!" << std::endl;
      std::cout << "Using Default Arguments:" << std::endl;
      // Size of the Matrix in consideration:
      N          = 6400;
      // Size of Matrices at leaf level:
      MinParticlesInLeaf          = 200;
      // Tolerance of problem
      TOL_POW  = 12;
  }

  else
  {
	N    =       atoi(argv[1]);  //      Number of particles.
	MinParticlesInLeaf   =       atoi(argv[2]); // minimum particles in each leaf of KD Tree
	TOL_POW = atoi(argv[3]);
}

////////////////////////////
Eigen::VectorXd b = Eigen::VectorXd::Random(N);
unsigned Dimension = 2;
Eigen::MatrixXd loc(N,Dimension);
for (size_t j = 0; j < N; j++) {
  for (size_t k = 0; k < Dimension; k++) {
    loc(j,k) = 2.0*double(rand())/double(RAND_MAX)-1.0;
  }
}
double start, end;
///////////////////////// AFMM /////////////////////////////
start	=	omp_get_wtime();
AFMM* afmm = new AFMM(N, MinParticlesInLeaf, TOL_POW, loc);
end	=	omp_get_wtime();
double timeCreateTreeAFMM = end-start;
std::cout << std::endl << "AFMM tree creation time: " << timeCreateTreeAFMM << std::endl;

std::cout << "========================= Problem Parameters =========================" << std::endl;
std::cout << "Matrix Size                        :" << N << std::endl;
std::cout << "Leaf Size                          :" << MinParticlesInLeaf << std::endl;
std::cout << "Tolerance                          :" << pow(10,-TOL_POW) << std::endl << std::endl;
unsigned nLevels = log(N/MinParticlesInLeaf)/log(4);
// std::cout << "nLevels: " << nLevels << std::endl;

start	=	omp_get_wtime();
afmm->assemble();
end	=	omp_get_wtime();
double timeAssemble = end-start;
std::cout << "========================= Assembly Time =========================" << std::endl;
std::cout << "Time for assemble in AFMM form    :" << timeAssemble << std::endl;


Eigen::VectorXd outputVec;
start = omp_get_wtime();
outputVec = afmm->computeMatVecProduct(b);
end   = omp_get_wtime();
double timeMatVecProduct = (end - start);

std::cout << "========================= Matrix-Vector Multiplication =========================" << std::endl;
std::cout << "Time for MatVec in HODLR form      :" << timeMatVecProduct << std::endl;
afmm->evaluateError();
}
