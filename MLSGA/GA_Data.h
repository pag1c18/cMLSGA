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
along with this program.If not, see < https://www.gnu.org/licenses/>. 

		GA_DATA header

Storage of additional classes, specific for MLSGA framework (without classes used for evolutioanry processes, as those are stored in class.h). Were implemented to test multiple problems in series and help with statistical data.


*/


#pragma once

//****************************************


#ifndef GA_FUNCTIONS_H
#define GA_FUNCTIONS_H

#include <vector>
#include <iostream>

#include "Const.h"
#include "Struct.h"
/**
Class for storing current GA parameters
*/
class GA_parameters
{
private:
	int pop_size;			///< population size
	int max_gen;			///< maximum generation
	float cross_prob;		///< crossover probability
	double mut_prob;			///< mutation probability
	short di_c;				///< distribution index of the crossover -  for crossover 1
	short di_m;				///< distribution index of the mutation - for mutation 1

public:
	/**Default constructor. Throws ERROR - GA_parameters cannot be empty**/
	GA_parameters()
	{
		std::cout << "ERROR GA_parameters creation!";
		system("pause");
		abort();
	}
	/**
	Default normal constructor
	@param - cprob crossover probability
	@param - mprob mutation probability
	@param - mgen number of maximum gnerations default value defined in Const.h
	@param psize - population size default value defined in Const.h
	@param dic - distribution index of the crossover. For crossover 1 - default value defined in Const.h
	@param dim - distribution index of the mutation. For mutation 1 - default value defined in Const.h
	*/
	GA_parameters(float cprob, double mprob, int psize, int mgen, short dic = Di_c, short dim = Di_m)
	{
		pop_size = psize;
		max_gen = mgen;
		cross_prob = cprob;
		mut_prob = mprob;
		di_c = dic;
		di_m = dim;
	};
	~GA_parameters() {}
	/**Returning population size**/
	int Pop_Size() { return pop_size; }
	/**Returning number of max generations**/
	int Max_Gen() { return max_gen; }
	/**Returning crossover probability**/
	float Cross_Prob() { return cross_prob; }
	/**Returning mutation probability**/
	double Mut_Prob() { return mut_prob; }
	/**Returning distribution index of the crossover for crossover 1**/
	short Di_C() { return di_c; }
	/**Returning distribution index of the mutation for mutation 1**/
	short Di_M() { return di_m; }

};

/**Class for storing output values and statistical data for a number of runs.
*/
template <typename tname>
class GA_data
{
private:
	std::vector<tname> time;	///<time of the run
	std::vector<tname> GA_time;	///<time of the GA
	std::vector<tname> IGD;		///<IGD for the run
	std::vector<tname> IGD2;		///<IGD for the run - for dynamic calcualted every gen not only during change like in IGD_val
	std::vector<tname> HV;		///<HV for the run
	std::vector<tname> HV2;		///<HV for the run - for dynamic calcualted every gen not only during change like in HV_val
	std::vector<tname> fitness;	///<Min fitness of the current run
	std::vector<int> generation;	///<generation at which min fitness was achieved
	std::vector<int> iteration;	///<generation at which min fitness was achieved

	bool IGD_on;	///<if IGD is calculated in the current run
	bool HV_on;		///<if HV is calculated in the current run

	STRUCTURES::min_max_avg_std<tname> time_struct; ///<stucture with min,max and average time
	STRUCTURES::min_max_avg_std<tname> GA_time_struct; ///<stucture with min,max and average time
	STRUCTURES::min_max_avg_std<tname> IGD_struct; ///<stucture with min,max and average IGD
	STRUCTURES::min_max_avg_std<tname> IGD2_struct; ///<stucture with min,max and average IGD2
	STRUCTURES::min_max_avg_std<tname> HV_struct; ///<stucture with min,max and average HV
	STRUCTURES::min_max_avg_std<tname> HV2_struct; ///<stucture with min,max and average HV2
	STRUCTURES::min_max_avg_std<tname> fitness_struct; ///<stucture with min,max and average fitness
	STRUCTURES::min_max_avg_std<int> generation_struct; ///<stucture with min,max and average number of generations
	STRUCTURES::min_max_avg_std<int> iteration_struct; ///<stucture with min,max and average number of generations
protected:
	/**Calculate the avg and std values**/
	void Average_Calc();
public:
	/**Default empty constructor. Throwing error.**/
	GA_data() { abort(); };
	/**Default constructor
	* @param IGD_o - if IGD is calculated during the run and thus should be stored
	* @param HV_o - if HV is calculated during the run and thus should be stored
	*/
	GA_data(bool IGD_o, bool HV_o) { IGD_on = IGD_o; HV_on = HV_o; };
	~GA_data() {};

	/**Adding values to the vectors and recalcualting the current min and max values*
	@param time_val - time for the current run
	@param GA_time_val - GA time for the current run
	@param IGD_val - IGD for the current run
	@param HV_val - HV for the current run
	@param generation_val - generation for the current run
	@param IGD2_val - IGD for the current run. For dynamic problems, where IGD is calcualted every generation and not only during the time change, unlike in IGD_val
	@param HV2_val - HV for the current run. For dynamic problems, where HV is calcualted every generation and not only during the time change, unlike in HV_val
	@param fitness_val - fitness for the current run
	*/
	void Add(tname time_val, tname GA_time_val, tname IGD_val, tname HV_val, int generation_val, int iteration, tname IGD2_val, tname HV2_val, tname fitness_val = 0 );

	
	/**Calculate the standard deviation**/
	void Std_Dev_Calculation();
	/**Returning time struct**/
	STRUCTURES::min_max_avg_std<tname> Show_Time_Struct() const { return time_struct; }
	/**Returning GA time structure**/
	STRUCTURES::min_max_avg_std<tname> Show_GA_Time_Struct() const { return GA_time_struct; }
	/**Returning IGD struct**/
	STRUCTURES::min_max_avg_std<tname> Show_IGD_Struct() const { return IGD_struct; }
	/**Returning IGD2 struct**/
	STRUCTURES::min_max_avg_std<tname> Show_IGD2_Struct() const { return IGD2_struct; }
	/**Returning HV struct**/
	STRUCTURES::min_max_avg_std<tname> Show_HV_Struct() const { return HV_struct; }
	/**Returning HV2 struct**/
	STRUCTURES::min_max_avg_std<tname> Show_HV2_Struct() const { return HV2_struct; }
	/**Returning fitness struct**/
	STRUCTURES::min_max_avg_std<tname> Show_Fitness_Struct() const { return fitness_struct; }
	/**Returning generation structure**/
	STRUCTURES::min_max_avg_std<int> Show_Generation_Struct() const { return generation_struct; }
	STRUCTURES::min_max_avg_std<int> Show_Iteration_Struct() const { return iteration_struct; }

};

template class GA_data<double>;
template class GA_data<float>;
template class GA_data<double long>;
#endif // !GA_FUNCTIONS_H
