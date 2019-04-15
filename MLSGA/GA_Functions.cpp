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


#include "GA_Functions.h"
#include <numeric>
#include "Define.h"

/*
*Adding values to the vectors and calcualting min and max values*
@param time_val time for the current run
@param GA_time_val GA time for the current run
@param IGD_val IGD for the current run
@param HV_val HV for the current run
@param generation_val generation for the current run
@param IGD2_val IGD for the current run - for dynamic calcualted every gen not only during change like in IGD_val
@param HV2_val HV for the current run - for dynamic calcualted every gen not only during change like in HV_val
@param fitness_val fitness for the current run
*/
template<typename tname>
void GA_data<tname>::Add(tname time_val, tname GA_time_val, tname IGD_val, tname HV_val, int generation_val, int iteration_val, tname IGD2_val, tname HV2_val, tname fitness_val = 0)
{
	//push values to the vectors
	time.push_back(time_val);
	GA_time.push_back(GA_time_val);
	if (ONE_OBJ_OVERRIDE != true)
	{	
		if (IGD_on == true)
		{
			IGD.push_back(IGD_val);
			IGD2.push_back(IGD2_val);
		}
		if (HV_on == true)
		{
			HV.push_back(HV_val);
			HV2.push_back(HV2_val);
		}
	}
	generation.push_back(generation_val);
	iteration.push_back(iteration_val);
	if (ONE_OBJ_OVERRIDE == true)
		fitness.push_back(fitness_val);

	//check if it is first values
	if (time.size() == 1)
	{
		//assign first values to both min and max
		time_struct.max = time_struct.min = time_val;
		GA_time_struct.max = GA_time_struct.min = GA_time_val;
		if (IGD_on == true)
		{
			IGD_struct.max = IGD_struct.min = IGD_val;
			IGD2_struct.max = IGD2_struct.min = IGD2_val;
		}
		if (HV_on == true)
		{
			HV_struct.max = HV_struct.min = HV_val;
			HV2_struct.max = HV2_struct.min = HV2_val;
		}
		generation_struct.max = generation_struct.min = generation_val;
		iteration_struct.max = iteration_struct.min = iteration_val;
		fitness_struct.max = fitness_struct.min = fitness_val;
	}
	else
	{
		//check if value is bigger then max
		if (time_val > time_struct.max)
			//assign as max
			time_struct.max = time_val;
		//check if value is lower than min
		else if (time_val < time_struct.min)
			//assign as min
			time_struct.min = time_val;
		
		//check if value is bigger then max
		if (GA_time_val > GA_time_struct.max)
			//assign as max
			GA_time_struct.max = GA_time_val;
		//check if value is lower than min
		else if (GA_time_val < GA_time_struct.min)
			//assign as min
			GA_time_struct.min = GA_time_val;

		//check if value is bigger then max
		if (generation_val > generation_struct.max)
			//assign as max
			generation_struct.max = generation_val;
		//check if value is lower than min
		else if (generation_val < generation_struct.min)
			//assign as min
			generation_struct.min = generation_val;

		//check if value is bigger then max
		if (iteration_val > iteration_struct.max)
			//assign as max
			iteration_struct.max = iteration_val;
		//check if value is lower than min
		else if (iteration_val < iteration_struct.min)
			//assign as min
			iteration_struct.min = iteration_val;

		//if multi objective calculate values for IGD - no IGD for one objective
		if (ONE_OBJ_OVERRIDE != true)
		{
			if (IGD_on == true)
			{
				//check if value is bigger then max
				if (IGD_val > IGD_struct.max)
					//assign as max
					IGD_struct.max = IGD_val;
				//check if value is lower than min
				else if (IGD_val < IGD_struct.min)
					//assign as min
					IGD_struct.min = IGD_val;
				//do the same for second IGD
				//check if value is bigger then max
				if (IGD2_val > IGD2_struct.max)
					//assign as max
					IGD2_struct.max = IGD2_val;
				//check if value is lower than min
				else if (IGD2_val < IGD2_struct.min)
					//assign as min
					IGD2_struct.min = IGD2_val;
			}
			if (HV_on == true)
			{
				//check if value is bigger then max
				if (HV_val > HV_struct.max)
					//assign as max
					HV_struct.max = HV_val;
				//check if value is lower than min
				else if (HV_val < HV_struct.min)
					//assign as min
					HV_struct.min = HV_val;
				//do the same for second HV
				//check if value is bigger then max
				if (HV2_val > HV2_struct.max)
					//assign as max
					HV2_struct.max = HV2_val;
				//check if value is lower than min
				else if (HV2_val < HV2_struct.min)
					//assign as min
					HV2_struct.min = HV2_val;
			}
		}

		//if one objective calculate values for fitness	- no one fitness value for multi objective
		if (ONE_OBJ_OVERRIDE == true)
		{
			//check if value is bigger then max
			if (fitness_val > fitness_struct.max)
				//assign as max
				fitness_struct.max = fitness_val;
			//check if value is lower than min
			else if (fitness_val < fitness_struct.min)
				//assign as min
				fitness_struct.min = fitness_val;
		}
	}
}

/**Calculate the average values**/
template<typename tname>
void GA_data<tname>::Average_Calc()
{
	//get the size of the data
	int size = time.size();

	//calculate average values
	time_struct.avg = std::accumulate(time.begin(), time.end(), 0.0) / (tname)size;
	GA_time_struct.avg = std::accumulate(GA_time.begin(), GA_time.end(), 0.0) / (tname)size;
	generation_struct.avg = std::accumulate(generation.begin(), generation.end(), 0.0) / (tname)size;
	iteration_struct.avg = std::accumulate(iteration.begin(), iteration.end(), 0.0) / (tname)size;

	//if multi objective calculate values for IGD - no IGD for one objective
	if (ONE_OBJ_OVERRIDE != true)
	{
		if (IGD_on)
		{
			IGD_struct.avg = std::accumulate(IGD.begin(), IGD.end(), 0.0) / (tname)size;
			IGD2_struct.avg = std::accumulate(IGD2.begin(), IGD2.end(), 0.0) / (tname)size;
		}
		if (HV_on)
		{
			HV_struct.avg = std::accumulate(HV.begin(), HV.end(), 0.0) / (tname)size;
			HV2_struct.avg = std::accumulate(HV2.begin(), HV2.end(), 0.0) / (tname)size;
		}
	}
	
	//if one objective calculate values for fitness	- no one fitness value for multi objective
	if (ONE_OBJ_OVERRIDE == true)
		fitness_struct.avg = std::accumulate(fitness.begin(), fitness.end(), 0.0) / (tname)size;

}

template<typename tname>
void GA_data<tname>::Std_Dev_Calculation()
{
	//calculate averages
	Average_Calc();

	//get the size of the data
	int size = time.size();

	tname time_std_temp = 0, GA_time_std_temp = 0, IGD_std_temp = 0, IGD2_std_temp = 0, HV_std_temp = 0, HV2_std_temp = 0, generation_std_temp = 0, iteration_std_temp = 0, fitness_std_temp = 0;		//temporary value for the std deviation calculation
	//calculate the sum of squares
	for (int i = 0; i < size; i++)
	{
		//add the squares to temp values
		time_std_temp += pow(time[i] - time_struct.avg, 2);
		GA_time_std_temp += pow(GA_time[i] - GA_time_struct.avg, 2);
		generation_std_temp += pow(generation[i] - generation_struct.avg, 2);
		iteration_std_temp += pow(iteration[i] - iteration_struct.avg, 2);

		//if multi objective calculate values for IGD - no IGD for one objective
		if (ONE_OBJ_OVERRIDE != true)
		{
			if (IGD_on)
			{
				IGD_std_temp += pow(IGD[i] - IGD_struct.avg, 2);
				IGD2_std_temp += pow(IGD2[i] - IGD2_struct.avg, 2);
			}
			if (HV_on)
			{
				HV_std_temp += pow(HV[i] - HV_struct.avg, 2);
				HV2_std_temp += pow(HV2[i] - HV2_struct.avg, 2);
			}
		}
		//if one objective calculate values for fitness	- no one fitness value for multi objective
		if (ONE_OBJ_OVERRIDE == true)
			fitness_std_temp += pow(fitness[i] - fitness_struct.avg, 2);
	}
	//calculate the standard deviation
	time_struct.std_deviation = (tname)sqrt(time_std_temp / (tname)size);
	GA_time_struct.std_deviation = (tname)sqrt(GA_time_std_temp / (tname)size);
	generation_struct.std_deviation = (tname)sqrt(generation_std_temp / (tname)size);
	iteration_struct.std_deviation = (tname)sqrt(iteration_std_temp / (tname)size);

	//if multi objective calculate values for IGD - no IGD for one objective
	if (ONE_OBJ_OVERRIDE != true)
	{
		if (IGD_on)
		{
			IGD_struct.std_deviation = (tname)sqrt(IGD_std_temp / (tname)size);
			IGD2_struct.std_deviation = (tname)sqrt(IGD2_std_temp / (tname)size);
		}
		if (HV_on)
		{
			HV_struct.std_deviation = (tname)sqrt(HV_std_temp / (tname)size);
			HV2_struct.std_deviation = (tname)sqrt(HV2_std_temp / (tname)size);
		}
	}

	//if one objective calculate values for fitness	- no one fitness value for multi objective
	if (ONE_OBJ_OVERRIDE == true)
		fitness_struct.std_deviation = (tname)sqrt(fitness_std_temp / (tname)size);

}
