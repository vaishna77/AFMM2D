#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <map>

#include "kernel.hpp"
#include "ACA.hpp"
#include "AFMM2DBox.hpp"
#include "AFMM2DTree.hpp"
#include "KDTree.hpp"

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
  double L;
	AFMM(int N, int MinParticlesInLeaf, int TOL_POW, Eigen::MatrixXd& loc, double L) {
		this->N = N;
    this->L = L;
		this->MinParticlesInLeaf = MinParticlesInLeaf;
		this->TOL_POW = TOL_POW;
		nLevels      = log(N/MinParticlesInLeaf)/log(4);
		n_Dimension  =       2;  //      Dimension.
		n_Properties =       1;  //      Number of properties.
		locations    =       new double [N*n_Dimension];   //      Stores all the locations.
		properties   =       new double [N*n_Properties];  //      Stores all the properties.
		// Generates random locations and random values of property.
		// change the locations as per your need
		unsigned count_Location =       0;
		unsigned count_Property =       0;
		for (unsigned j=0; j<N; ++j) {
						for (unsigned k=0; k<n_Dimension; ++k) {
										locations[count_Location] = loc(j,k);//       =       2*double(rand())/double(RAND_MAX)-1;//2*(int(rand())%2)-1;//double(rand())/double(RAND_MAX);
										++count_Location;
						}
						for (unsigned k=0; k<n_Properties; ++k) {
										properties[count_Property]      =       j; //indices 0,1,2,3; after sorting this will be the indices of sorted array
										++count_Property;
						}
		}
		//////////////////////////////// KD Tree //////////////////////////////////////////
		std::vector<std::vector<int> > boxNumbers; // boxNumbers[nLevels] contains box numbers in N ordering. ex: [0 3 1 2]
		sorted_Locations        =       new double[N*n_Dimension];
		sorted_Properties       =       new double[N*n_Properties];
		std::vector<int> NumberOfParticlesInLeaves;// contains number of particels in each leaf in N ordering
		sort_KDTree(N, n_Dimension, locations, n_Properties, properties, MinParticlesInLeaf, nLevels, sorted_Locations, sorted_Properties, boxNumbers, NumberOfParticlesInLeaves); // Creates a KDTree given the locations. This KD Tree class generates a uniform tree - all leaves are at level nLevels. Number of particles in boxes at a given level differ by atmost 1.
		std::vector<pts2D> particles_X, particles_Y;
		mykernel		=	new userkernel(particles_X, particles_Y);
		afmm2dtree	=	new FMM2DTree(mykernel, int(N), int(nLevels), L, TOL_POW, sorted_Locations, boxNumbers[nLevels], NumberOfParticlesInLeaves);
	}

	void assemble() {
		afmm2dtree->createTree();
		afmm2dtree->assign_Tree_Interactions();
		afmm2dtree->assign_Center_Location();
		afmm2dtree->assignLeafChargeLocations();
		afmm2dtree->getNodes();
		afmm2dtree->assemble_M2L();
	}

	Eigen::VectorXd computeMatVecProduct(Eigen::VectorXd inputVecUnsorted) {
 		Eigen::VectorXd outputVec;
    Eigen::VectorXd inputVec(N); //converting array to Eigen vector
		for (size_t i = 0; i < N; i++) {
      inputVec(i) = inputVecUnsorted(int(sorted_Properties[i]));
		}
		afmm2dtree->assignLeafCharges(inputVec);
		afmm2dtree->evaluate_M2M();
		afmm2dtree->evaluate_M2L();
		afmm2dtree->evaluate_L2L();
		afmm2dtree->evaluate_NearField();
		outputVec = Eigen::VectorXd(N);
		Eigen::VectorXd outputVecSorted;
		outputVecSorted = afmm2dtree->getMatVecProductOutput();//potential sorted
		outputVec = Eigen::VectorXd(N);
		for (size_t i = 0; i < N; i++) {
			outputVec(int(sorted_Properties[i])) = outputVecSorted(i);
		}
		return outputVec;
	}

	void evaluateError() {
		srand(time(NULL));
    std::cout << std::endl << "Performing error calculation in box: " << std::endl;
    for (size_t i = 0; i < 20; i++) {
      int nBox	=	rand()%afmm2dtree->nBoxesPerLevel[nLevels];
      std::cout << "nBox: " << nBox << "	err: " << afmm2dtree->compute_error(nBox) << std::endl;
    }
	}

	~AFMM() {
		delete locations;
		delete properties;
		delete sorted_Locations;
		delete sorted_Properties;
		delete afmm2dtree;
	};
};

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
Eigen::VectorXd b(N); //vector definition
for (size_t i = 0; i < N; i++) {
  b(i) = 2*double(rand())/double(RAND_MAX)-1;
}
double L = 1.0;
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
AFMM* afmm = new AFMM(N, MinParticlesInLeaf, TOL_POW, loc, L);
end	=	omp_get_wtime();
double timeCreateTreeAFMM = end-start;
std::cout << std::endl << "AFMM tree creation time: " << timeCreateTreeAFMM << std::endl;

std::cout << "========================= Problem Parameters =========================" << std::endl;
std::cout << "Matrix Size                        :" << N << std::endl;
std::cout << "Leaf Size                          :" << MinParticlesInLeaf << std::endl;
std::cout << "Tolerance                          :" << pow(10,-TOL_POW) << std::endl << std::endl;
unsigned nLevels = log(N/MinParticlesInLeaf)/log(4);
std::cout << "nLevels: " << nLevels << std::endl;

start	=	omp_get_wtime();
afmm->assemble();
end	=	omp_get_wtime();
double timeAssemble = end-start;
std::cout << "========================= Assembly Time =========================" << std::endl;
std::cout << "Time for assemble in AFMM form    :" << timeAssemble << std::endl;


// What we are doing here is explicitly generating
// the matrix from its entries
start = omp_get_wtime();
Eigen::MatrixXd B = afmm->afmm2dtree->K->getMatrix(0, 0, N, N);
end   = omp_get_wtime();
double exact_time = (end - start);
std::cout << "Time for direct matrix generation  :" << exact_time << std::endl;
std::cout << "Magnitude of Speed-Up              :" << (exact_time / timeAssemble) << std::endl << std::endl;

Eigen::VectorXd outputVec;
start = omp_get_wtime();
outputVec = afmm->computeMatVecProduct(b);
end   = omp_get_wtime();
double timeMatVecProduct = (end - start);
// afmm->evaluateError();

std::cout << "========================= Matrix-Vector Multiplication =========================" << std::endl;
std::cout << "Time for MatVec in HODLR form      :" << timeMatVecProduct << std::endl;

start = omp_get_wtime();
Eigen::VectorXd bSorted(N);
for (size_t i = 0; i < N; i++) {
  bSorted(i) = b(int(afmm->sorted_Properties[i]));
}
Eigen::MatrixXd r_exact_Sorted = B * bSorted;
Eigen::VectorXd r_exact(N);
for (size_t i = 0; i < N; i++) {
  r_exact(int(afmm->sorted_Properties[i])) = r_exact_Sorted(i);
}

end   = omp_get_wtime();
exact_time = (end - start);
std::cout << "Time for direct MatVec             :" << exact_time << std::endl;
std::cout << "Magnitude of Speed-Up              :" << (exact_time / timeMatVecProduct) << std::endl;
// Computing the relative error in the solution obtained:
std::cout << "Error in the solution is           :" << (outputVec-r_exact).norm() / (r_exact.norm()) << std::endl << std::endl;
}
