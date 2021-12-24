//
//  KDTree.hpp
//
//
//  Created by Sivaram Ambikasaran on 12/3/13.
//
//

#ifndef __KDTREE_HPP__
#define __KDTREE_HPP__

#include <iostream>
#include <bits/stdc++.h>
#include <Eigen/Dense>

class KDTree {
private:
        unsigned n_Locations;                   //      Number of locations.
        unsigned n_Dimension;                   //      Number of dimensions.
        unsigned n_Properties;                  //      Number of properties.
        unsigned MinParticlesInLeaf;
        unsigned nLevels;

        Eigen::MatrixXd sorted_Contents;        //      Array of sorted contents, that contains both the locations and properties.

        //      Merges two sorted lists.
        //      Usage:
        //      merge_Sorted_Lists(7, 20, 9, 50, 0);
        void merge_Sorted_Lists(unsigned n_Left_Start, unsigned n_Left_Size, unsigned n_Right_Start, unsigned n_Right_Size, unsigned n_Index);

        //      Performs a merge sort on a list given the starting index, size of the set and index based on which the sorting needs to be done.
        //      Usage:
        //      merge_Sort(11, 5, 1);
        void merge_Sort(unsigned n_Start, unsigned n_Size, unsigned n_Index);

        //      Given the starting index, size of the set of points and the index (which has to be less than n_Dimension) based on which to sort, this function performs a sort based on KDTree.
        //      Usage:
        //      sort_KDTree(7, 193, 2);
        int sort_KDTree(unsigned n_Start, unsigned n_Size, unsigned n_Index, unsigned local_level, std::vector<int>& NumberOfParticlesInLeaves);

public:
        //      Constructor.
        //      Usage:
        //      KDTree(10000, 3, locations);
        KDTree(const unsigned n_Locations, const unsigned n_Dimension, double* locations, const unsigned n_Properties, double* properties, const unsigned MinParticlesInLeaf, const unsigned nLevels);

        //      KDTree sorting.
        int sort_KDTree(std::vector<int>& NumberOfParticlesInLeaves);

        //      Obtains all the sorted locations and properties.
        //      Usage:
        //      get_Location(location, properties);
        //      Obtains all the locations and properties in sorted fashion.
        void get_Location_Properties(double* locations, double* properties);

        //      Given an index, obtains the location and property from the sorted location and property.
        //      n_Index <= n_Locations.
        //      Usage:
        //      get_Location_Properties(1234, location, properties);
        //      Obtains the location and properties corresponding to the 1234 index in the sorted fashion.
        void get_Location_Properties(const unsigned n_Index, double* location, double* properties);
};
void display(std::string display_String, unsigned n_Locations, unsigned n_Dimension, double* locations, unsigned n_Properties, double* properties);
void sort_KDTree(unsigned N, unsigned n_Dimension, double* locations, unsigned n_Properties, double* properties, unsigned MinParticlesInLeaf, unsigned nLevels, double* sorted_Locations, double* sorted_Properties, std::vector<std::vector<int> >& boxNumbers, std::vector<int>& NumberOfParticlesInLeaves);
#endif /* defined(__KDTREE_HPP__) */
