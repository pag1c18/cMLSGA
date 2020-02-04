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

/*The IGD calculation is based on P. A. N. Bosman and D. Thierens, “The balance between proximity and diversity in multiobjective evolutionary algorithms,” IEEE Trans. Evol. Comput., vol. 7, no. 2, pp. 174–188, 2003.*/


#include "IGD.h"
#include "Support_Functions.h"
#include <numeric>
/*
*IGD calculation*
@param PF Pareto Front for which IGD will be calculated
@param real_PF real Pareto Front for IGD calculation
*/
double IGD_calc(pareto_front & PF, std::vector<std::vector<double>> & real_PF)
{
	//copy the pareto front individuals to the storage vector
	std::vector<individual> Pareto = PF.Indiv_Show();		//storage vector of the pareto front individuals
	
	//get the size of the pareto front
	int Pareto_size = Pareto.size();						//size of the pareto front

	if (Pareto_size == 0)
		return INF;

	//get the number of objectives
	//short Obj_num = Pareto[0].Fitness_Show().size();		//number of objectives

	double IGD_val;											//IGD value - output
	
	std::vector<double>di_min_GA;							//storage vector with min distance, between point from real PF and each point from given PF, for each real PF point

	//calculate the min distance
	for (int i = 0; i < real_PF.size(); i++)			//for each real PF point
	{
		//calculate the distance to each point of given PF and get the min one
		for (int j = 0; j < Pareto_size; j++)
		{
			double temp_min_di;					//temp min distance
			//double temp = 0;					//distance of the current point to the real PF point
			
			//calcualte the distance
			/*for (short n = 0; n < Obj_num; n++)
			{
				temp += pow(real_PF[i][n] - Pareto[j].Fitness_Show(n), 2);
			}*/
			temp_min_di = Distance(Pareto[j].Fitness_Show(true), real_PF[i]);
			//chack if it is the min one
			if (j == 0)
				di_min_GA.push_back(temp_min_di);
			else
			{
				if (di_min_GA[i] > temp_min_di)
					di_min_GA[i] = temp_min_di;
			}
		}
	}
	//calculate the IGD
	IGD_val = std::accumulate(di_min_GA.begin(), di_min_GA.end(), 0.0);
	IGD_val /= real_PF.size();

	return IGD_val;
}