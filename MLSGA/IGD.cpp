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
			temp_min_di = Distance(Pareto[j].Fitness_Show(), real_PF[i]);
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