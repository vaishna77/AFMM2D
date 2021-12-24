//
//  KDTree.cpp
//
//
//  Created by Sivaram Ambikasaran on 12/3/13.
//
//

#include "KDTree.hpp"
// #include "../include/KDTree.hpp"
#include <Eigen/Dense>

KDTree::KDTree(const unsigned n_Locations, const unsigned n_Dimension, double* locations, const unsigned n_Properties, double* properties, const unsigned MinParticlesInLeaf, const unsigned nLevels) {
        this->n_Locations       =       n_Locations;
        this->n_Dimension       =       n_Dimension;
        this->n_Properties      =       n_Properties;
        this->MinParticlesInLeaf=       MinParticlesInLeaf;
        this->nLevels           =       nLevels;

        sorted_Contents         =       Eigen::MatrixXd(n_Locations, n_Dimension+n_Properties);

        unsigned count_Location =       0;
        unsigned count_Property =       0;

        for (unsigned j=0; j<n_Locations; ++j) {
                for (unsigned k=0; k<n_Dimension; ++k) {
                        sorted_Contents(j,k)    =       locations[count_Location];
                        ++count_Location;
                }
                for (unsigned k=n_Dimension; k<n_Dimension+n_Properties; ++k) {
                        sorted_Contents(j,k)    =       properties[count_Property];
                        ++count_Property;
                }
        }
}

void KDTree::merge_Sorted_Lists(unsigned n_Left_Start, unsigned n_Left_Size, unsigned n_Right_Start, unsigned n_Right_Size, unsigned n_Index) {

        Eigen::MatrixXd temp_List(n_Left_Size+n_Right_Size, n_Dimension+n_Properties);

        unsigned j_Left =       n_Left_Start;
        unsigned j_Right=       n_Right_Start;

        unsigned j      =       0;

        while (j_Left < n_Left_Start + n_Left_Size && j_Right < n_Right_Start + n_Right_Size) {
                if (sorted_Contents(j_Left, n_Index) < sorted_Contents(j_Right, n_Index)) {
                        temp_List.row(j)=       sorted_Contents.row(j_Left);
                        ++j_Left;
                }
                else {
                        temp_List.row(j)=       sorted_Contents.row(j_Right);
                        ++j_Right;
                }
                ++j;
        }

        while (j_Left < n_Left_Start + n_Left_Size) {
                temp_List.row(j)        =       sorted_Contents.row(j_Left);
                ++j_Left;
                ++j;
        }

        while (j_Right < n_Right_Start + n_Right_Size) {
                temp_List.row(j)        =       sorted_Contents.row(j_Right);
                ++j_Right;
                ++j;
        }

        sorted_Contents.block(n_Left_Start, 0, n_Left_Size + n_Right_Size, n_Dimension+n_Properties)        =       temp_List;
}

void KDTree::merge_Sort(unsigned n_Start, unsigned n_Size, unsigned n_Index) {
        if (n_Size<=1) {
                //      Do nothing
                return;
        }
        else {
                unsigned n_Left_Start   =       n_Start;
                unsigned n_Left_Size    =       n_Size/2;
                unsigned n_Right_Start  =       n_Start+n_Left_Size;
                unsigned n_Right_Size   =       n_Size-n_Left_Size;

                merge_Sort(n_Left_Start, n_Left_Size, n_Index);
                merge_Sort(n_Right_Start, n_Right_Size, n_Index);

                merge_Sorted_Lists(n_Left_Start, n_Left_Size, n_Right_Start, n_Right_Size, n_Index);
        }
}

int KDTree::sort_KDTree(unsigned n_Start, unsigned n_Size, unsigned n_Index, unsigned local_level, std::vector<int>& NumberOfParticlesInLeaves) {
    // static int nLevels = 0;
    n_Index                 =       n_Index%n_Dimension;
    merge_Sort(n_Start, n_Size, n_Index);

    unsigned n_Left_Start   =       n_Start;
    unsigned n_Left_Size    =       n_Size/2;
    unsigned n_Right_Start  =       n_Start+n_Left_Size;
    unsigned n_Right_Size   =       n_Size-n_Left_Size;

    ++n_Index;
    n_Index                 =       n_Index%n_Dimension;

    merge_Sort(n_Left_Start, n_Left_Size, n_Index);
    unsigned n_Left_Bottom_Start   =       n_Left_Start;
    unsigned n_Left_Bottom_Size    =       n_Left_Size/2;
    unsigned n_Left_Top_Start  =       n_Left_Start+n_Left_Bottom_Size;
    unsigned n_Left_Top_Size   =       n_Left_Size-n_Left_Bottom_Size;

    merge_Sort(n_Right_Start, n_Right_Size, n_Index);
    unsigned n_Right_Bottom_Start   =       n_Right_Start;
    unsigned n_Right_Bottom_Size    =       n_Right_Size/2;
    unsigned n_Right_Top_Start  =       n_Right_Start+n_Right_Bottom_Size;
    unsigned n_Right_Top_Size   =       n_Right_Size-n_Right_Bottom_Size;

    ++n_Index;
    n_Index                 =       n_Index%n_Dimension;

    if (local_level < nLevels) {
    // if (n_Left_Bottom_Size >= MinParticlesInLeaf &&  n_Left_Top_Size >= MinParticlesInLeaf && n_Right_Bottom_Size >= MinParticlesInLeaf && n_Right_Top_Size >= MinParticlesInLeaf && local_level < nLevels) {
      local_level++;
      if(local_level==nLevels) {
        // std::cout << n_Left_Bottom_Size << std::endl;
        // std::cout << n_Left_Top_Size << std::endl;
        // std::cout << n_Right_Bottom_Size << std::endl;
        // std::cout << n_Right_Top_Size << std::endl;
        // std::cout << "--------------------------------" << std::endl;
        NumberOfParticlesInLeaves.push_back(n_Left_Bottom_Size);
        NumberOfParticlesInLeaves.push_back(n_Left_Top_Size);
        NumberOfParticlesInLeaves.push_back(n_Right_Bottom_Size);
        NumberOfParticlesInLeaves.push_back(n_Right_Top_Size);
      }
      // std::cout << "local_level: " << local_level << std::endl;
      this->sort_KDTree(n_Left_Bottom_Start, n_Left_Bottom_Size, n_Index, local_level, NumberOfParticlesInLeaves);
      this->sort_KDTree(n_Left_Top_Start, n_Left_Top_Size, n_Index, local_level, NumberOfParticlesInLeaves);
      this->sort_KDTree(n_Right_Bottom_Start, n_Right_Bottom_Size, n_Index, local_level, NumberOfParticlesInLeaves);
      this->sort_KDTree(n_Right_Top_Start, n_Right_Top_Size, n_Index, local_level, NumberOfParticlesInLeaves);
      // nLevels += 1;
    }
    return 0;
    // return nLevels;
}

int KDTree::sort_KDTree(std::vector<int>& NumberOfParticlesInLeaves) {
        if (n_Locations<=1) {
                //      Do nothing
                return 0;
        }
        else {
                //      Number of point on the left cluster
                return this->sort_KDTree(0, n_Locations, 0, 0, NumberOfParticlesInLeaves);
        }
}

void KDTree::get_Location_Properties(double* locations, double* properties) {
        unsigned count_Location         =       0;
        unsigned count_Property         =       0;
        for (unsigned j=0; j<n_Locations; ++j) {
                for (unsigned k=0; k<n_Dimension; ++k) {
                        locations[count_Location]       =       sorted_Contents(j,k);
                        ++count_Location;
                }
                for (unsigned k=n_Dimension; k<n_Dimension+n_Properties; ++k) {
                        properties[count_Property]      =       sorted_Contents(j,k);
                        ++count_Property;
                }
        }
}

void KDTree::get_Location_Properties(const unsigned n_Index, double* location, double* properties) {
        for (unsigned k=0; k<n_Dimension; ++k) {
                location[k]     =       sorted_Contents(n_Index,k);
        }
        for (unsigned k=0; k<n_Properties; ++k) {
                properties[k]   =       sorted_Contents(n_Index,k+n_Dimension);
        }
}

/////////////////GLOBAL FUNCTIONS TO BE CALLED FROM MAIN ///////////////////////////
void display(std::string display_String, unsigned n_Locations, unsigned n_Dimension, double* locations, unsigned n_Properties, double* properties) {
        std::cout << display_String << std::endl;

        unsigned count_Location  =       0;
        unsigned count_Property  =       0;

        for (unsigned j=0; j<n_Locations; ++j) {
                for (unsigned k=0; k<n_Dimension; ++k) {
                        std::cout << locations[count_Location] << "\t";
                        ++count_Location;
                }
                for (unsigned k=0; k<n_Properties; ++k) {
                        std::cout << properties[count_Property] << "\t";
                        ++count_Property;
                }
                std::cout << std::endl;
        }
        std::cout << std::endl;
}

void sort_KDTree(unsigned N, unsigned n_Dimension, double* locations, unsigned n_Properties, double* properties, unsigned MinParticlesInLeaf, unsigned nLevels, double* sorted_Locations, double* sorted_Properties, std::vector<std::vector<int> >& boxNumbers, std::vector<int>& NumberOfParticlesInLeaves) {
    KDTree* B                       =       new KDTree(N, n_Dimension, locations, n_Properties, properties, MinParticlesInLeaf, nLevels);
    delete locations;
    delete properties;
    NumberOfParticlesInLeaves.reserve(pow(4,nLevels));
    //      Sorts the locations based on KDTree.
    B->sort_KDTree(NumberOfParticlesInLeaves);
    // converting N ordering to the mirrored C ordering.
    // The sorted locations that the KD Tree outputs, correspond to boxes in N ordering.
    // So here we are generating N ordering sequence of box numbers in terms of the mirrored C ordering.
    // only leaf level ordering is what we want: boxNumbers[nLevels]
    std::vector<int> boxNumbersInLevel0;
    boxNumbersInLevel0.push_back(0);//level 0
    boxNumbers.push_back(boxNumbersInLevel0);
    for (size_t j = 1; j <= nLevels; j++) {
      std::vector<int> boxNumbersInALevel;
      for (size_t k = 0; k < boxNumbers[j-1].size(); k++) {
        boxNumbersInALevel.push_back(boxNumbers[j-1][k]*4 + 0);
        boxNumbersInALevel.push_back(boxNumbers[j-1][k]*4 + 3);
        boxNumbersInALevel.push_back(boxNumbers[j-1][k]*4 + 1);
        boxNumbersInALevel.push_back(boxNumbers[j-1][k]*4 + 2);
      }
      boxNumbers.push_back(boxNumbersInALevel);
    }
    // Obtains the sorted location.
    // sorted_Locations contains the locations sorted as per the KD Tree.(In N ordering of boxes)
    B->get_Location_Properties(sorted_Locations, sorted_Properties);

    // Display the sorted contents.
   // display("Sorted contents: ", N, n_Dimension, sorted_Locations, n_Properties, sorted_Properties);
   // exit(0);
   // Take away from kD tree class: sorted_Locations, sorted_Properties, boxNumbers[nLevels], NumberOfParticlesInLeaves. Give these as inputs to FMM2DTree constructor.
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////

    /* example:
    boxNumbers NumberOfParticlesInLeaves
    0           4
    3           4
    1           4
    2           4
    */
    delete B;
}
