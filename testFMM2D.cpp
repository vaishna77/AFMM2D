#include "kernel.hpp"
#include "ACA.hpp"
#include "FMM2DTree.hpp"
#include "KDTree.cpp"

int main(int argc, char* argv[]) {
	unsigned N    =       atoi(argv[1]);  //      Number of particles.
	unsigned MinParticlesInLeaf   =       atoi(argv[2]); // minimum particles in each leaf of KD Tree
	int TOL_POW = atoi(argv[3]);
	double start, end;
	unsigned nLevels = log(N/MinParticlesInLeaf)/log(4);
  unsigned n_Dimension    =       2;  //      Dimension.
	unsigned n_Properties   =       0;  //      Number of properties.
  double* locations       =       new double [N*n_Dimension];   //      Stores all the locations.
	double* properties      =       new double [N*n_Properties];  //      Stores all the properties.

	//      Generates random locations and random values of property.
  // change the locations as per your need
	unsigned count_Location =       0;
	unsigned count_Property =       0;
	for (unsigned j=0; j<N; ++j) {
					for (unsigned k=0; k<n_Dimension; ++k) {
									locations[count_Location]       =       2*double(rand())/double(RAND_MAX)-1;//2*(int(rand())%2)-1;//double(rand())/double(RAND_MAX);
									++count_Location;
					}
					for (unsigned k=0; k<n_Properties; ++k) {
									properties[count_Property]      =       double(rand())/double(RAND_MAX);
									++count_Property;
					}
	}

  //////////////////////////////// KD Tree //////////////////////////////////////////
  std::vector<std::vector<int> > boxNumbers; // boxNumbers[nLevels] contains box numbers in N ordering. ex: [0 3 1 2]
  double* sorted_Locations        =       new double[N*n_Dimension];
  double* sorted_Properties       =       new double[N*n_Properties];
  std::vector<int> NumberOfParticlesInLeaves;// contains number of particels in each leaf in N ordering
  // Creates a KDTree given the locations. This KD Tree class generates a uniform tree - all leaves are at level nLevels. Number of particles in boxes at a given level differ by atmost 1.
  sort_KDTree(N, n_Dimension, locations, n_Properties, properties, MinParticlesInLeaf, nLevels, sorted_Locations, sorted_Properties, boxNumbers, NumberOfParticlesInLeaves);

	///////////////////////////// FMM2D ////////////////////////////////////////////
	start	=	omp_get_wtime();
	std::vector<pts2D> particles_X, particles_Y;
	userkernel* mykernel		=	new userkernel(particles_X, particles_Y);
	FMM2DTree<userkernel>* A	=	new FMM2DTree<userkernel>(mykernel, int(N), int(nLevels), TOL_POW, sorted_Locations, boxNumbers[nLevels], NumberOfParticlesInLeaves);

	A->createTree();
	A->assign_Tree_Interactions();

	end		=	omp_get_wtime();
	double timeCreateTree	=	(end-start);

	std::cout << std::endl << "Number of particles is: " << A->N << std::endl;
	std::cout << std::endl << "Time taken to create the tree is: " << timeCreateTree << std::endl;
	/////////////////////////////////////////////////////////////////////////
	start	=	omp_get_wtime();

	A->assign_Center_Location();
	A->assignLeafChargeLocations();
	A->getNodes();
	A->assemble_M2L();
	// int N = A->N;
	Eigen::VectorXd b(N);
	for (size_t i = 0; i < N; i++) {
		b(i) = A->K->chargesFunction(A->gridPoints[i]);
	}
	A->assignLeafCharges(b);

	end		=	omp_get_wtime();
	double timeAssignCharges=	(end-start);
	std::cout << std::endl << "Time taken to assemble is: " << timeAssignCharges << std::endl;
	/////////////////////////////////////////////////////////////////////////
	start	=	omp_get_wtime();
	A->evaluate_M2M();
	A->evaluate_M2L();
	A->evaluate_L2L();
	A->evaluate_NearField();
	end		=	omp_get_wtime();
	double timeMatVecProduct=	(end-start);
	std::cout << std::endl << "Time taken to do Mat-Vec product is: " << timeMatVecProduct << std::endl;
	/////////////////////////////////////////////////////////////////////////
	srand(time(NULL));
	int nBox	=	rand()%A->nBoxesPerLevel[nLevels];
	std::cout << std::endl << "Performing error calculation in box: " << nBox << std::endl;
	std::cout << "err: " << A->compute_error(nBox) << std::endl;

  delete sorted_Locations;
  delete sorted_Properties;
  delete A;
}
