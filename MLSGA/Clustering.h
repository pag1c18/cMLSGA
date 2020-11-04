/**
Copyright(C) 2019  Przemyslaw A.Grudniewski and Adam J.Sobey

This file is part of the MLSGA framework

The MLSGA framework is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

The MLSGA framework is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see < https://www.gnu.org/licenses/>. */

#pragma once
#ifndef CLUSTERING_H
#define CLUSTERING_H

#include "Class.h"
#include <vector>


/**
*k-means Clustering, calculate the labels for the population separation*
@param pop population for which labels will be calculated
*/
std::vector<short> Clustering(const population & pop);
/**
*Copying the data from the population to a vector*
@param pop source population
*/
static void Data_Set(const population & pop);												
/**
*Setting up initial cluster points*
@param pop source population
*/
static void C_Points_Set(const population & pop);											



#endif // !CLUSTERING_H
