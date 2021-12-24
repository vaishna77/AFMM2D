#include "AFMM2DTree.hpp"

FMM2DTree::FMM2DTree(userkernel* K, int N, int nLevels, double L, int TOL_POW, double* locations, std::vector<int>& boxNumbers, std::vector<int>& NumberOfParticlesInLeaves) {
  this->K					=	K;
  this->nLevels		=	nLevels;
  this->L					=	L;
  this->locations = locations;
  this->boxNumbers= boxNumbers;
  this->NumberOfParticlesInLeaves = NumberOfParticlesInLeaves;
  // this->nParticlesInLeafAlong1D = nParticlesInLeafAlong1D;
  // this->nParticlesInLeaf = nParticlesInLeafAlong1D*nParticlesInLeafAlong1D;
  this->TOL_POW = TOL_POW;
  nBoxesPerLevel.push_back(1);
  boxRadius.push_back(L);
  for (int k=1; k<=nLevels; ++k) {
    nBoxesPerLevel.push_back(4*nBoxesPerLevel[k-1]);
    boxRadius.push_back(0.5*boxRadius[k-1]);
  }
  this->smallestBoxSize	=	boxRadius[nLevels];
  K->a					=	smallestBoxSize;
  // this->N					=	nParticlesInLeaf*nBoxesPerLevel[nLevels];
  this->N					= N;
}

void FMM2DTree::createTree() {
  //	First create root and add to tree
  FMM2DBox root;
  root.boxNumber		=	0;
  root.parentNumber	=	-1;
  #pragma omp parallel for
  for (int l=0; l<4; ++l) {
    root.childrenNumbers[l]	=	l;
  }
  #pragma omp parallel for
  for (int l=0; l<8; ++l) {
    root.neighborNumbers[l]	=	-1;
  }
  #pragma omp parallel for
  for (int l=0; l<16; ++l) {
    root.innerNumbers[l]	=	-1;
  }
  #pragma omp parallel for
  for (int l=0; l<24; ++l) {
    root.outerNumbers[l]	=	-1;
  }
  std::vector<FMM2DBox> rootLevel;
  rootLevel.push_back(root);
  tree.push_back(rootLevel);

  for (int j=1; j<=nLevels; ++j) {
    std::vector<FMM2DBox> level;
    for (int k=0; k<nBoxesPerLevel[j]; ++k) {
      FMM2DBox box;
      box.boxNumber		=	k;
      box.parentNumber	=	k/4;
      for (int l=0; l<4; ++l) {
        box.childrenNumbers[l]	=	4*k+l;
      }
      level.push_back(box);
    }
    tree.push_back(level);
  }
}

//	Assigns the interactions for child0 of a box
void FMM2DTree::assign_Child0_Interaction(int j, int k) {
  int nL	=	j+1;
  int nC	=	4*k;
  int nN;

  //	Assign siblings
  {
    /************************/
    /*	 _____________		*/
    /*	|	   |	  |		*/
    /*	|  N5  |  N4  |		*/
    /*	|______|______|		*/
    /*	|	   |	  |		*/
    /*	|  **  |  N3  |		*/
    /*	|______|______|		*/
    /*						*/
    /************************/
    tree[nL][nC].neighborNumbers[3]	=	nC+1;
    tree[nL][nC].neighborNumbers[4]	=	nC+2;
    tree[nL][nC].neighborNumbers[5]	=	nC+3;
  }

  //	Assign children of parent's zeroth neighbor
  {
    /****************************/
    /*				   ______	*/
    /*				  |		 |	*/
    /*				  |	 **  |	*/
    /*	 _____________|______|  */
    /*	|	   |	  |			*/
    /*	|  I15 |  N0  |			*/
    /*	|______|______|			*/
    /*	|	   |	  |			*/
    /*	|  I0  |  I1  |			*/
    /*	|______|______|			*/
    /*							*/
    /****************************/
    nN	=	tree[j][k].neighborNumbers[0];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[15]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's first neighbor
  {
    /****************************/
    /*	 ______			  		*/
    /*	|	   |	  			*/
    /*	|  **  |				*/
    /*	|______|______			*/
    /*	|	   |	  |			*/
    /*	|  N1  |  N2  |			*/
    /*	|______|______|			*/
    /*	|	   |	  |			*/
    /*	|  I2  |  I3  |			*/
    /*	|______|______|			*/
    /*							*/
    /****************************/
    nN	=	tree[j][k].neighborNumbers[1];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].neighborNumbers[1]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's second neighbor
  {
    /************************************/
    /*	 ______			  				*/
    /*	|	   |	  					*/
    /*	|  **  |						*/
    /*	|______|	   _____________	*/
    /*				  |	     |	    |	*/
    /*				  |  I5  |  O8  |	*/
    /*				  |______|______|	*/
    /*				  |	     |	    |	*/
    /*				  |  I4  |  O7  |	*/
    /*				  |______|______|	*/
    /*									*/
    /************************************/
    nN	=	tree[j][k].neighborNumbers[2];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].outerNumbers[7]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].outerNumbers[8]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's third neighbor
  {
    /************************************/
    /*				   _____________	*/
    /*				  |	     |	    |	*/
    /*				  |  I7  |  O10 |	*/
    /*	 ______		  |______|______|	*/
    /*	|	   |	  |	     |	    |	*/
    /*	|  **  |	  |  I6  |  O9  |	*/
    /*	|______|	  |______|______|	*/
    /*									*/
    /************************************/
    nN	=	tree[j][k].neighborNumbers[3];
    if (nN!=-1) {
      tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].outerNumbers[9]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].outerNumbers[10]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[3];
    }

  }

  //	Assign children of parent's fourth neighbor
  {
    /************************************/
    /*				   _____________	*/
    /*				  |	     |	    |	*/
    /*				  |  O13 |  O12 |	*/
    /*				  |______|______|	*/
    /*				  |	     |	    |	*/
    /*		    	  |  I8  |  O11 |	*/
    /*		    	  |______|______|	*/
    /*									*/
    /*									*/
    /*	 ______							*/
    /*  |      |						*/
    /*  |  **  |						*/
    /*  |______|						*/
    /*  								*/
    /************************************/
    nN	=	tree[j][k].neighborNumbers[4];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].outerNumbers[11]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].outerNumbers[12]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].outerNumbers[13]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's fifth neighbor
  {
    /************************/
    /*	 _____________		*/
    /*	|	   |	  |		*/
    /*	|  O15 |  O14 |		*/
    /*	|______|______|		*/
    /*	|	   |	  |		*/
    /*	|  I10 |  I9  |		*/
    /*	|______|______|		*/
    /*						*/
    /*						*/
    /*	 ______				*/
    /*  |	   |			*/
    /*	|  **  |			*/
    /*	|______|			*/
    /*  					*/
    /************************/
    nN	=	tree[j][k].neighborNumbers[5];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].outerNumbers[14]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].outerNumbers[15]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's sixth neighbor
  {
    /****************************/
    /*	 _____________			*/
    /*	|	   |	  |			*/
    /*	|  O17 |  O16 |			*/
    /*	|______|______|			*/
    /*	|	   |	  |			*/
    /*	|  I12 |  I11 |			*/
    /*	|______|______|			*/
    /*							*/
    /*							*/
    /*				   ______	*/
    /*  			  |		 |	*/
    /*				  |	 **  |	*/
    /*				  |______|	*/
    /*  						*/
    /****************************/
    nN	=	tree[j][k].neighborNumbers[6];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[12]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].outerNumbers[16]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].outerNumbers[17]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's seventh neighbor
  {
    /****************************/
    /*							*/
    /*	 _____________			*/
    /*	|	   |	  |			*/
    /*	|  I13 |  N6  |			*/
    /*	|______|______|______	*/
    /*  |	   |	  |		 |	*/
    /*	|  I14 |  N7  |	 **  |	*/
    /*	|______|______|______|	*/
    /*  						*/
    /****************************/
    nN	=	tree[j][k].neighborNumbers[7];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[14]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].neighborNumbers[7]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[13]	=	tree[j][nN].childrenNumbers[3];
    }
  }
}

//	Assigns the interactions for child1 of a box
void FMM2DTree::assign_Child1_Interaction(int j, int k) {
  int nL	=	j+1;
  int nC	=	4*k+1;
  int nN;

  //	Assign siblings
  {
    /************************/
    /*	 _____________		*/
    /*	|	   |	  |		*/
    /*	|  N6  |  N5  |		*/
    /*	|______|______|		*/
    /*	|	   |	  |		*/
    /*	|  N7  |  **  |		*/
    /*	|______|______|		*/
    /*						*/
    /************************/
    tree[nL][nC].neighborNumbers[7]	=	nC-1;
    tree[nL][nC].neighborNumbers[5]	=	nC+1;
    tree[nL][nC].neighborNumbers[6]	=	nC+2;
  }

  //	Assign children of parent's zeroth neighbor
  {
    /************************************/
    /*				   		  ______	*/
    /*				  	     |		|	*/
    /*				         |	**  |	*/
    /*	 _____________       |______|  	*/
    /*	|	   |	  |					*/
    /*	|  O22 |  I15 |					*/
    /*	|______|______|					*/
    /*	|	   |	  |					*/
    /*	|  O23 |  I0  |					*/
    /*	|______|______|					*/
    /*									*/
    /************************************/
    nN	=	tree[j][k].neighborNumbers[0];
    if (nN != -1) {
      tree[nL][nC].outerNumbers[23]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[15]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].outerNumbers[22]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's first neighbor
  {
    /****************************/
    /*	 		______		  	*/
    /*		   |	  |			*/
    /*	       |  **  |			*/
    /*	 ______|______|			*/
    /*	|	   |	  |			*/
    /*	|  N0  |  N1  |			*/
    /*	|______|______|			*/
    /*	|	   |	  |			*/
    /*	|  I1  |  I2  |			*/
    /*	|______|______|			*/
    /*							*/
    /****************************/
    nN	=	tree[j][k].neighborNumbers[1];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].neighborNumbers[1]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's second neighbor
  {
    /****************************/
    /*	 ______		  			*/
    /*	|	   |				*/
    /*	|  **  |	  			*/
    /*	|______|_____________	*/
    /*		   |	  |	     |	*/
    /*		   |  N2  |  I5  |	*/
    /*		   |______|______|	*/
    /*		   |	  |	     |	*/
    /*	       |  I3  |  I4  |	*/
    /*		   |______|______|	*/
    /*							*/
    /****************************/
    nN	=	tree[j][k].neighborNumbers[2];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's third neighbor
  {
    /****************************/
    /*			_____________	*/
    /*		   |	  |	     |	*/
    /*		   |  N4  |	 I7	 |  */
    /*	 ______|______|______|	*/
    /*	|	   |	  |	     |	*/
    /*	|  **  |  N3  |  I6  |  */
    /*	|______|______|______| 	*/
    /*							*/
    /****************************/
    nN	=	tree[j][k].neighborNumbers[3];
    if (nN != -1){
      tree[nL][nC].neighborNumbers[3]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].neighborNumbers[4]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's fourth neighbor
  {
    /****************************/
    /*			_____________	*/
    /*		   |	  |	     |	*/
    /*		   |  O14 |  O13 |	*/
    /*		   |______|______|	*/
    /*		   |	  |	     |	*/
    /*		   |  I9  |  I8  |	*/
    /*		   |______|______|	*/
    /*				  			*/
    /*				  			*/
    /*	 ______					*/
    /*	|	   |				*/
    /*  |  **  |				*/
    /*  |______|				*/
    /*  						*/
    /****************************/
    nN	=	tree[j][k].neighborNumbers[4];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].outerNumbers[13]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].outerNumbers[14]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's fifth neighbor
  {
    /************************/
    /*	 _____________		*/
    /*	|	   |	  |		*/
    /*	|  O16 |  O15 |		*/
    /*	|______|______|		*/
    /*	|	   |	  |		*/
    /*	|  I11 |  I10 |		*/
    /*	|______|______|		*/
    /*						*/
    /*						*/
    /*		    ______		*/
    /* 		   |	  |		*/
    /*		   |  **  |		*/
    /*		   |______|		*/
    /*  					*/
    /************************/
    nN	=	tree[j][k].neighborNumbers[5];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].outerNumbers[15]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].outerNumbers[16]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's sixth neighbor
  {
    /************************************/
    /*	 _____________					*/
    /*	|	   |	  |					*/
    /*	|  O18 |  O17 |					*/
    /*	|______|______|					*/
    /*	|	   |	  |					*/
    /*	|  O19 |  I12 |					*/
    /*	|______|______|					*/
    /*									*/
    /*									*/
    /*				   		  ______	*/
    /*  			  		 |		|	*/
    /*				  		 |	** 	|	*/
    /*				  		 |______|	*/
    /*  								*/
    /************************************/
    nN	=	tree[j][k].neighborNumbers[6];
    if (nN != -1) {
      tree[nL][nC].outerNumbers[19]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[12]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].outerNumbers[17]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].outerNumbers[18]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's seventh neighbor
  {
    /************************************/
    /*									*/
    /*	 _____________					*/
    /*	|	   |	  |					*/
    /*	|  O20 |  I13 |					*/
    /*	|______|______|		  ______	*/
    /*  |	   |	  |		 |		|	*/
    /*	|  O21 |  I14 |	 	 |	**  |	*/
    /*	|______|______|		 |______|	*/
    /*  								*/
    /************************************/
    nN	=	tree[j][k].neighborNumbers[7];
    if (nN != -1) {
      tree[nL][nC].outerNumbers[21]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[14]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[13]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].outerNumbers[20]	=	tree[j][nN].childrenNumbers[3];
    }
  }
}

//	Assigns the interactions for child2 of a box
void FMM2DTree::assign_Child2_Interaction(int j, int k) {
  int nL	=	j+1;
  int nC	=	4*k+2;
  int nN;

  //	Assign siblings
  {
    /************************/
    /*	 _____________		*/
    /*	|	   |	  |		*/
    /*	|  N7  |  **  |		*/
    /*	|______|______|		*/
    /*	|	   |	  |		*/
    /*	|  N0  |  N1  |		*/
    /*	|______|______|		*/
    /*						*/
    /************************/
    tree[nL][nC].neighborNumbers[0]	=	nC-2;
    tree[nL][nC].neighborNumbers[1]	=	nC-1;
    tree[nL][nC].neighborNumbers[7]	=	nC+1;
  }

  //	Assign children of parent's zeroth neighbor
  {
    /************************************/
    /*				   		  ______	*/
    /*				  	     |		|	*/
    /*				         |	**  |	*/
    /*				         |______|  	*/
    /*									*/
    /*									*/
    /*	 _____________					*/
    /*	|	   |	  |					*/
    /*	|  O23 |  I0  |					*/
    /*	|______|______|					*/
    /*	|	   |	  |					*/
    /*	|  O0  |  O1  |					*/
    /*	|______|______|					*/
    /*									*/
    /************************************/
    nN	=	tree[j][k].neighborNumbers[0];
    if (nN != -1) {
      tree[nL][nC].outerNumbers[0]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].outerNumbers[1]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].outerNumbers[23]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's first neighbor
  {
    /****************************/
    /*	 		______		  	*/
    /*		   |	  |			*/
    /*	       |  **  |			*/
    /*	 	   |______|			*/
    /*							*/
    /*							*/
    /*	 _____________			*/
    /*	|	   |	  |			*/
    /*	|  I1  |  I2  |			*/
    /*	|______|______|			*/
    /*	|	   |	  |			*/
    /*	|  O2  |  O3  |			*/
    /*	|______|______|			*/
    /*							*/
    /****************************/
    nN	=	tree[j][k].neighborNumbers[1];
    if (nN != -1) {
      tree[nL][nC].outerNumbers[2]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].outerNumbers[3]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's second neighbor
  {
    /****************************/
    /*	 ______		  			*/
    /*	|	   |				*/
    /*	|  **  |	  			*/
    /*	|______|				*/
    /*							*/
    /*							*/
    /*			_____________	*/
    /*		   |	  |	     |	*/
    /*		   |  I3  |  I4  |	*/
    /*		   |______|______|	*/
    /*		   |	  |	     |	*/
    /*	       |  O4  |  O5  |	*/
    /*		   |______|______|	*/
    /*							*/
    /****************************/
    nN	=	tree[j][k].neighborNumbers[2];
    if (nN != -1) {
      tree[nL][nC].outerNumbers[4]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].outerNumbers[5]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's third neighbor
  {
    /****************************/
    /*	 ____________________	*/
    /*	|	   |	  |	     |	*/
    /*	|  **  |  N3  |	 I6	 |  */
    /*	|______|______|______|	*/
    /*		   |	  |	     |	*/
    /*		   |  N2  |  I5  |  */
    /*		   |______|______| 	*/
    /*							*/
    /****************************/
    nN	=	tree[j][k].neighborNumbers[3];
    if (nN != -1) {
      tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].neighborNumbers[3]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's fourth neighbor
  {
    /****************************/
    /*			_____________	*/
    /*		   |	  |	     |	*/
    /*		   |  I9  |  I8  |	*/
    /*		   |______|______|	*/
    /*		   |	  |	     |	*/
    /*		   |  N4  |  I7  |	*/
    /*	 ______|______|______|	*/
    /*	|	   |	  			*/
    /*	|  **  |	  			*/
    /*	|______|	  			*/
    /*							*/
    /****************************/
    nN	=	tree[j][k].neighborNumbers[4];
    if (nN != -1) {
      tree[nL][nC].neighborNumbers[4]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's fifth neighbor
  {
    /************************/
    /*	 _____________		*/
    /*	|	   |	  |		*/
    /*	|  I11 |  I10 |		*/
    /*	|______|______|		*/
    /*	|	   |	  |		*/
    /*	|  N6  |  N5  |		*/
    /*	|______|______|		*/
    /*		   |	  |		*/
    /*		   |  **  |		*/
    /*		   |______|		*/
    /*  					*/
    /************************/
    nN	=	tree[j][k].neighborNumbers[5];
    if (nN != -1) {
      tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].neighborNumbers[5]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's sixth neighbor
  {
    /************************************/
    /*	 _____________					*/
    /*	|	   |	  |					*/
    /*	|  O19 |  I12 |					*/
    /*	|______|______|					*/
    /*	|	   |	  |					*/
    /*	|  O20 |  I13 |					*/
    /*	|______|______|		  ______	*/
    /*  			  		 |		|	*/
    /*				  		 |	** 	|	*/
    /*				  		 |______|	*/
    /*  								*/
    /************************************/
    nN	=	tree[j][k].neighborNumbers[6];
    if (nN != -1) {
      tree[nL][nC].outerNumbers[20]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[13]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[12]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].outerNumbers[19]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's seventh neighbor
  {
    /************************************/
    /*									*/
    /*	 _____________		  ______	*/
    /*	|	   |	  |		 |	    |	*/
    /*	|  O21 |  I14 |		 |	**	|	*/
    /*	|______|______|		 |______|	*/
    /*  |	   |	  |		 			*/
    /*	|  O22 |  I15 |	 	 			*/
    /*	|______|______|		 			*/
    /*  								*/
    /************************************/
    nN	=	tree[j][k].neighborNumbers[7];
    if (nN != -1) {
      tree[nL][nC].outerNumbers[22]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[15]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[14]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].outerNumbers[21]	=	tree[j][nN].childrenNumbers[3];
    }
  }
}

//	Assigns the interactions for child3 of a box
void FMM2DTree::assign_Child3_Interaction(int j, int k) {
  int nL	=	j+1;
  int nC	=	4*k+3;
  int nN;

  //	Assign siblings
  {
    /************************/
    /*	 _____________		*/
    /*	|	   |	  |		*/
    /*	|  **  |  N3  |		*/
    /*	|______|______|		*/
    /*	|	   |	  |		*/
    /*	|  N1  |  N2  |		*/
    /*	|______|______|		*/
    /*						*/
    /************************/
    tree[nL][nC].neighborNumbers[1]	=	nC-3;
    tree[nL][nC].neighborNumbers[2]	=	nC-2;
    tree[nL][nC].neighborNumbers[3]	=	nC-1;
  }

  //	Assign children of parent's zeroth neighbor
  {
    /****************************/
    /*				   ______	*/
    /*				  |		 |	*/
    /*				  |	 **  |	*/
    /*				  |______|  */
    /*							*/
    /*							*/
    /*	 _____________			*/
    /*	|	   |	  |			*/
    /*	|  I0  |  I1  |			*/
    /*	|______|______|			*/
    /*	|	   |	  |			*/
    /*	|  O1  |  O2  |			*/
    /*	|______|______|			*/
    /*							*/
    /****************************/
    nN	=	tree[j][k].neighborNumbers[0];
    if (nN != -1) {
      tree[nL][nC].outerNumbers[1]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].outerNumbers[2]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's first neighbor
  {
    /****************************/
    /*	 ______		  			*/
    /*	|	   |				*/
    /*	|  **  |				*/
    /*	|______|				*/
    /*							*/
    /*							*/
    /*	 _____________			*/
    /*	|	   |	  |			*/
    /*	|  I2  |  I3  |			*/
    /*	|______|______|			*/
    /*	|	   |	  |			*/
    /*	|  O3  |  O4  |			*/
    /*	|______|______|			*/
    /*							*/
    /****************************/
    nN	=	tree[j][k].neighborNumbers[1];
    if (nN != -1) {
      tree[nL][nC].outerNumbers[3]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].outerNumbers[4]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's second neighbor
  {
    /************************************/
    /*	 ______		  					*/
    /*	|	   |						*/
    /*	|  **  |	  					*/
    /*	|______|						*/
    /*									*/
    /*									*/
    /*				   _____________	*/
    /*		   		  |	     |	    |	*/
    /*		   		  |  I4  |  O7  |	*/
    /*		   		  |______|______|	*/
    /*		   		  |	  	 |	    |	*/
    /*	       		  |  O5  |  O6  |	*/
    /*		   		  |______|______|	*/
    /*									*/
    /************************************/
    nN	=	tree[j][k].neighborNumbers[2];
    if (nN != -1) {
      tree[nL][nC].outerNumbers[5]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].outerNumbers[6]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].outerNumbers[7]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's third neighbor
  {
    /************************************/
    /*	 ______		   _____________	*/
    /*	|	   |	  |	     |		|	*/
    /*	|  **  |      |	 I6	 |  O9	|	*/
    /*	|______|	  |______|______|	*/
    /*		   		  |	  	 |	    |	*/
    /*		   		  |  I5  |  O8  |  	*/
    /*		   		  |______|______| 	*/
    /*									*/
    /************************************/
    nN	=	tree[j][k].neighborNumbers[3];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].outerNumbers[8]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].outerNumbers[9]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's fourth neighbor
  {
    /************************************/
    /*				   _____________	*/
    /*		   		  |	  	 |	    |	*/
    /*		   		  |  I8  |  O11 |	*/
    /*		   		  |______|______|	*/
    /*		   		  |	  	 |	    |	*/
    /*		   		  |  I7  |  O10 |	*/
    /*	 ______	      |______|______|	*/
    /*	|	   |	  					*/
    /*	|  **  |	  					*/
    /*	|______|	  					*/
    /*									*/
    /************************************/
    nN	=	tree[j][k].neighborNumbers[4];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].outerNumbers[10]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].outerNumbers[11]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's fifth neighbor
  {
    /************************/
    /*	 _____________		*/
    /*	|	   |	  |		*/
    /*	|  I10 |  I9  |		*/
    /*	|______|______|		*/
    /*	|	   |	  |		*/
    /*	|  N5  |  N4  |		*/
    /*	|______|______|		*/
    /*	|	   |			*/
    /*	|  **  |			*/
    /*	|______|			*/
    /*  					*/
    /************************/
    nN	=	tree[j][k].neighborNumbers[5];
    if (nN != -1) {
      tree[nL][nC].neighborNumbers[5]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].neighborNumbers[4]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's sixth neighbor
  {
    /****************************/
    /*	 _____________			*/
    /*	|	   |	  |			*/
    /*	|  I12 |  I11 |			*/
    /*	|______|______|			*/
    /*	|	   |	  |			*/
    /*	|  I13 |  N6  |			*/
    /*	|______|______|______	*/
    /*  			  |		 |	*/
    /*				  |	 **  |	*/
    /*				  |______|	*/
    /*  						*/
    /****************************/
    nN	=	tree[j][k].neighborNumbers[6];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[13]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].neighborNumbers[6]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[12]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's seventh neighbor
  {
    /****************************/
    /*							*/
    /*	 ____________________	*/
    /*	|	   |	  |		 |	*/
    /*	|  I14 |  N7  |	 **	 |	*/
    /*	|______|______|______|	*/
    /*  |	   |	  |		 	*/
    /*	|  I15 |  N0  |	 	 	*/
    /*	|______|______|		 	*/
    /*  						*/
    /****************************/
    nN	=	tree[j][k].neighborNumbers[7];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[15]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].neighborNumbers[7]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[14]	=	tree[j][nN].childrenNumbers[3];
    }
  }
}


//	Assigns the interactions for the children of a box
void FMM2DTree::assign_Box_Interactions(int j, int k) {
  this->assign_Child0_Interaction(j,k);
  this->assign_Child1_Interaction(j,k);
  this->assign_Child2_Interaction(j,k);
  this->assign_Child3_Interaction(j,k);
}

//	Assigns the interactions for the children all boxes at a given level
void FMM2DTree::assign_Level_Interactions(int j) {
  #pragma omp parallel for
  for (int k=0; k<nBoxesPerLevel[j]; ++k) {
    this->assign_Box_Interactions(j,k);
  }
}

//	Assigns the interactions for the children all boxes in the tree
void FMM2DTree::assign_Tree_Interactions() {
  for (int j=0; j<nLevels; ++j) {
    this->assign_Level_Interactions(j);
  }
}

void FMM2DTree::assign_Center_Location() {
  int J;
  tree[0][0].center.x	=	0.0;
  tree[0][0].center.y	=	0.0;
  for (int j=0; j<nLevels; ++j) {
    J	=	j+1;
    double shift	=	0.5*boxRadius[j];
    #pragma omp parallel for
    for (int k=0; k<nBoxesPerLevel[j]; ++k) {
      tree[J][4*k].center.x		=	tree[j][k].center.x-shift;
      tree[J][4*k+1].center.x	=	tree[j][k].center.x+shift;
      tree[J][4*k+2].center.x	=	tree[j][k].center.x+shift;
      tree[J][4*k+3].center.x	=	tree[j][k].center.x-shift;

      tree[J][4*k].center.y		=	tree[j][k].center.y-shift;
      tree[J][4*k+1].center.y	=	tree[j][k].center.y-shift;
      tree[J][4*k+2].center.y	=	tree[j][k].center.y+shift;
      tree[J][4*k+3].center.y	=	tree[j][k].center.y+shift;
    }
  }
}

void FMM2DTree::getNodes() {
  for (int j=nLevels; j>=2; j--) {
    this->getNodes_incoming_level(j);
  }
}

void FMM2DTree::getParticlesFromChildren_incoming(int j, int k, std::vector<int>& searchNodes) {
  if (j==nLevels) {
    searchNodes.insert(searchNodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
  }
  else {
    int J = j+1;
    for (int c = 0; c < 4; c++) {
      searchNodes.insert(searchNodes.end(), tree[J][4*k+c].incoming_checkPoints.begin(), tree[J][4*k+c].incoming_checkPoints.end());
    }
  }
}

void FMM2DTree::getNodes_incoming_box(int j, int k, int& n_rows, int& n_cols, int& ComputedRank) {
  std::vector<int> boxA_Nodes;
  this->getParticlesFromChildren_incoming(j, k, boxA_Nodes);
  std::vector<int> IL_Nodes;//indices
  for(int in=0; in<16; ++in) {
    if(tree[j][k].innerNumbers[in] != -1) {
      int kIL = tree[j][k].innerNumbers[in];
      std::vector<int> chargeLocations;
      this->getParticlesFromChildren_incoming(j, kIL, chargeLocations);
      IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
    }
  }
  for(int on=0; on<24; ++on) {
    if(tree[j][k].outerNumbers[on] != -1) {
      int kIL = tree[j][k].outerNumbers[on];
      std::vector<int> chargeLocations;
      this->getParticlesFromChildren_incoming(j, kIL, chargeLocations);
      IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
    }
  }
  n_rows = boxA_Nodes.size();
  n_cols = IL_Nodes.size();
  std::vector<int> row_indices, col_indices;
  row_indices = boxA_Nodes;//object of base class G_LowRank
  col_indices = IL_Nodes;
  std::vector<int> row_bases, col_bases;
  Eigen::MatrixXd Ac, Ar;
  LowRank* LR		=	new LowRank(K, TOL_POW, row_indices, col_indices);
  LR->ACA_only_nodes(row_bases, col_bases, ComputedRank, Ac, Ar);
  int minN = n_rows;
  if (n_rows > n_cols) {
    minN = n_cols;
  }
  if(ComputedRank > 0) {
    for (int r = 0; r < row_bases.size(); r++) {
      tree[j][k].incoming_checkPoints.push_back(boxA_Nodes[row_bases[r]]);
    }
    for (int c = 0; c < col_bases.size(); c++) {
      tree[j][k].incoming_chargePoints.push_back(IL_Nodes[col_bases[c]]);
    }
    std::vector<int> row_indices_local;
    for (size_t r = 0; r < row_bases.size(); r++) {
      row_indices_local.push_back(boxA_Nodes[row_bases[r]]);
    }
    std::vector<int> col_indices_local;
    for (size_t c = 0; c < col_bases.size(); c++) {
      col_indices_local.push_back(IL_Nodes[col_bases[c]]);
    }
    Eigen::MatrixXd D = K->getMatrix(row_indices_local, col_indices_local);
    Eigen::PartialPivLU<Eigen::MatrixXd> D_T = Eigen::PartialPivLU<Eigen::MatrixXd>(D.transpose());
    tree[j][k].L2P = D_T.solve(Ac.transpose()).transpose();
  }
}

void FMM2DTree::getNodes_incoming_level(int j) {
  int ComputedRank, n_rows, n_cols;
  for (int k=0; k<nBoxesPerLevel[j]; ++k) {
    this->getNodes_incoming_box(j, k, n_rows, n_cols, ComputedRank);
  }
}

void FMM2DTree::assemble_M2L() {
  #pragma omp parallel for
  for (size_t j = 2; j <= nLevels; j++) {
    #pragma omp parallel for
    for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
      #pragma omp parallel for
      for(int in=0; in<16; ++in) {
        if(tree[j][k].innerNumbers[in] != -1) {
          int kIL = tree[j][k].innerNumbers[in];
          if (tree[j][k].M2L[kIL].size() == 0)
            tree[j][k].M2L[kIL] = K->getMatrix(tree[j][k].incoming_checkPoints, tree[j][kIL].incoming_checkPoints);
          if (tree[j][kIL].M2L[k].size() == 0)
            tree[j][kIL].M2L[k] = tree[j][k].M2L[kIL].transpose();
        }
      }
      #pragma omp parallel for
      for(int on=0; on<24; ++on) {
        if(tree[j][k].outerNumbers[on] != -1) {
          int kIL = tree[j][k].outerNumbers[on];
          if (tree[j][k].M2L[kIL].size() == 0)
            tree[j][k].M2L[kIL] = K->getMatrix(tree[j][k].incoming_checkPoints, tree[j][kIL].incoming_checkPoints);
          if (tree[j][kIL].M2L[k].size() == 0)
            tree[j][kIL].M2L[k] = tree[j][k].M2L[kIL].transpose();
        }
      }
    }
  }
}

void FMM2DTree::assignLeafChargeLocations() {
  for (size_t i = 0; i < N*2; i+=2) {
    pts2D temp;
    temp.x = locations[i];
    temp.y = locations[i+1];
    gridPoints.push_back(temp);
  }
  int startIndex = 0;
  for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
    int boxNum = boxNumbers[k];
    int NumParticles = NumberOfParticlesInLeaves[k];
    for (size_t i = 0; i < NumParticles; i++) {
      tree[nLevels][boxNum].chargeLocations.push_back(startIndex+i);
    }
    startIndex += NumParticles;
  }
  K->particles_X = gridPoints;//object of base class FMM_Matrix
  K->particles_Y = gridPoints;

}

void FMM2DTree::assignLeafCharges(Eigen::VectorXd &charges) {
  int start = 0;
  for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
    int boxNum = boxNumbers[k];
    tree[nLevels][boxNum].charges	=	charges.segment(start, tree[nLevels][boxNum].chargeLocations.size());
    start += tree[nLevels][boxNum].chargeLocations.size();
  }
}

void FMM2DTree::evaluate_M2M() {
  for (size_t j = nLevels; j >= 2; j--) {
    #pragma omp parallel for
    for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
      Eigen::VectorXd source_densities;
      if (j==nLevels) {
        source_densities = tree[j][k].charges;
      }
      else {
        int J = j+1;
        int Veclength = 0;
        for (int child = 0; child < 4; child++) {
          Veclength += tree[J][4*k+child].outgoing_charges.size();
        }
        source_densities = Eigen::VectorXd::Zero(Veclength);// = tree[j][k].multipoles//source densities
        int start = 0;
        for (int child = 0; child < 4; child++) {
          int NumElem = tree[J][4*k+child].outgoing_charges.size();
          source_densities.segment(start, NumElem) = tree[J][4*k+child].outgoing_charges;
          start += NumElem;
        }
      }
      tree[j][k].outgoing_charges = tree[j][k].L2P.transpose()*source_densities;//f^{B,o} //solve system: A\tree[j][k].outgoing_potential
    }
  }
}

void FMM2DTree::evaluate_M2L() {
  #pragma omp parallel for
  for (int j=2; j<=nLevels; ++j) {
    #pragma omp parallel for
    for (int k=0; k<nBoxesPerLevel[j]; ++k) {//BoxA
      tree[j][k].incoming_potential	=	Eigen::VectorXd::Zero(tree[j][k].incoming_checkPoints.size());
      #pragma omp parallel for
      for(int in=0; in<16; ++in) {
        if(tree[j][k].innerNumbers[in] != -1) {
          int kIL = tree[j][k].innerNumbers[in];
          tree[j][k].incoming_potential += tree[j][k].M2L[kIL]*tree[j][kIL].outgoing_charges;
        }
      }
      #pragma omp parallel for
      for(int on=0; on<24; ++on) {
        if(tree[j][k].outerNumbers[on] != -1) {
          int kIL = tree[j][k].outerNumbers[on];
          tree[j][k].incoming_potential += tree[j][k].M2L[kIL]*tree[j][kIL].outgoing_charges;
        }
      }
    }
  }
}

void FMM2DTree::evaluate_L2L() {
  for (size_t j = 2; j <= nLevels; j++) {
    #pragma omp parallel for
    for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
      if (j != nLevels) {
        Eigen::VectorXd temp = tree[j][k].L2P*tree[j][k].incoming_potential;
        int start = 0;
        for (size_t c = 0; c < 4; c++) {
          tree[j+1][4*k+c].incoming_potential += temp.segment(start, tree[j+1][4*k+c].incoming_checkPoints.size());
          start += tree[j+1][4*k+c].incoming_checkPoints.size();
        }
      }
      else {
        tree[j][k].potential = tree[j][k].L2P*tree[j][k].incoming_potential;
      }
    }
  }
}

void FMM2DTree::evaluate_NearField() { // evaluating at chargeLocations
  #pragma omp parallel for
  for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
    #pragma omp parallel for
    for (size_t n = 0; n < 8; n++) {
      int nn = tree[nLevels][k].neighborNumbers[n];
      if(nn != -1) {
        Eigen::MatrixXd R = K->getMatrix(tree[nLevels][k].chargeLocations, tree[nLevels][nn].chargeLocations);
        tree[nLevels][k].potential += R*tree[nLevels][nn].charges;
      }
    }
    Eigen::MatrixXd R = K->getMatrix(tree[nLevels][k].chargeLocations, tree[nLevels][k].chargeLocations);
    tree[nLevels][k].potential += R*tree[nLevels][k].charges; //self Interaction
  }
}

Eigen::VectorXd FMM2DTree::getMatVecProductOutput() { // evaluating at chargeLocations
  Eigen::VectorXd r(N);
  int offset = 0;
  for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
    int boxNum = boxNumbers[k];
    // std::cout << boxNum << std::endl;
    r.segment(offset,tree[nLevels][boxNum].potential.size()) = tree[nLevels][boxNum].potential;
    offset += tree[nLevels][boxNum].potential.size();
  }
  return r;
}

double FMM2DTree::compute_error(int nBox) { // evaluating at chargeLocations
  Eigen::VectorXd true_potential = Eigen::VectorXd::Zero(tree[nLevels][nBox].chargeLocations.size());
  for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
    Eigen::MatrixXd R = K->getMatrix(tree[nLevels][nBox].chargeLocations, tree[nLevels][k].chargeLocations);
    true_potential += R*tree[nLevels][k].charges;
  }
  Eigen::VectorXd errVec = true_potential - tree[nLevels][nBox].potential;
  double error = errVec.norm()/true_potential.norm();
  return error;
}
