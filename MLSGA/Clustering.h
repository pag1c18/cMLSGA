#pragma once
#ifndef CLUSTERING_H
#define CLUSTERING_H

#include "Class.h"
#include <vector>


/*
*k-means Clustering, calculate the labels for the population separation*
@param pop population for which labels will be calculated
*/
std::vector<short> Clustering(const population & pop);
/*
*Copying the data from the population to a vector*
@param pop source population
*/
static void Data_Set(const population & pop);												
/*
*Setting up initial cluster points*
@param pop source population
*/
static void C_Points_Set(const population & pop);											



#endif // !CLUSTERING_H
