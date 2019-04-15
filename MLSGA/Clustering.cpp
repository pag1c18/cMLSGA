/*
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

#include "Clustering.h"
#include <fstream>
//#include "Const.h"
//#include "Struct.h"
//#include <algorithm>

static std::vector<std::vector<double>> c_data;						//data for clustering 
static std::vector<std::vector<double>> c_points;					//data for cluster points

extern short n_col;

/*
*k-means Clustering, calculate the labels for population separation*
@param pop population for which labels will be calculated
*/
std::vector<short> Clustering(const population & pop)
{
	//create empty vector for labels
	//std::vector<int> label{ (int)c_data.size(), 0 };				//label vector - output

	std::vector<short> label(pop.Size_Show(), 0);

	std::ifstream file;
	file.open("Input/Labels.txt");
	for (int i = 0; i < pop.Size_Show(); i++)
	{
		file >> label[i];
	}
	if (n_col != 6)
		abort();



	/*
	//copy the data to the data vector
	Data_Set(pop);

	//Set starting points
	C_Points_Set(pop);
	//int c_data_size = c_data.size();

	
	
	//calculate real label values
	*/

	return label;
}
/*
*Copying data from population to vector*
@param pop source population
*/
void Data_Set(const population & pop)
{
	//copy the data of each individual to a vector
	for (int i = 0; i < pop.Size_Show(); i++)
	{
		c_data.push_back(pop.Indiv_Show(i).Code_Show());
	}
}
/*
*Setting up initial cluster points*
@param pop source population
*/
//work on it!
void C_Points_Set(const population & pop)
{
	STRUCTURES::boundaries d_max_min{ 0,100000 };		//boundaries structure with max and min values for stating points
	
	//find min and max values
	for (int i = 0; i < c_data.size(); i++)
	{
		auto max_min = std::minmax_element(c_data[i].begin(), c_data[i].end());
		if (d_max_min.lower > *max_min.first)
			d_max_min.lower = *max_min.first;
		if (d_max_min.upper < *max_min.second)
			d_max_min.upper = *max_min.second;
	}

	//set up stating points
	for (int i = 0; i < n_col; i++)
	{
		std::vector<double> code;			//code of the statring point
		
		//set up code of the stating point (by random)
		for (int j = 0; j < c_data[0].size(); j++)
		{
			double x = Random();
			code.push_back(d_max_min.lower + x * (d_max_min.upper - d_max_min.lower));
		}
		//push the new stating point to the vector
		c_points.push_back(code);
	}
}

