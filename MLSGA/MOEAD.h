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

		MOEAD header
		Functions specific to MOEAD algorithm


*/


#pragma once
#ifndef MOEAD_H
#define MOEAD_H
#include "Class.h"

/**Functions specific to MOEAD algorithm*/
namespace MOEAD
{
	/**
	Evolve the individuals using MOEAD algorithm
	@param Col - address of a given collective
	@param iGen - current generation index
	*/
	std::vector<individual> MOEAD_Calc(collective & col, int iGen);


	/**
	Initialise the algorithm specific parameters, for a whole population.
	@param n_obj - current number of objectives
	@param pop - current (overall) population
	*/
	void MOEAD_Init(short n_obj, population & pop);


	/**
	Initialise the algorithm specific parameters, for a single collective.
	@param n_obj - current number of objectives
	@param col - current (overall) population
	*/
	void Population_Init(short n_obj, collective & col);
	/**
	Initialise the neighbourhood of closest solutions, for a single collective.
	@param col - current (overall) population
	*/
	void Neighbourhood_Init(collective & col);

	/**
	* Calculate the utility values for a single collective
	* @param Col - address of a given collective
	* @param iGen - current generation index
	*/
	void Utility_Comp(collective & col, int iGen);

	/**
	* Perform the differential evolution crossover on a given individuals. Result is a single offspring.
	* @param ind0 - parent individual
	* @param ind1 - parent individual
	* @param ind2 - parent individual
	* @param rate - rate of crossover
	* @param fcode - function according to which the process will be conducted
	* @param ccode - crossover class from which the crossover parameters will be taken
	* @param gapara - ga specific parameters
	*/
	individual Diff_Evo_XoverB(individual &ind0, individual &ind1, individual &ind2, double rate, function & fcode, crossover<individual> & ccode, GA_parameters & gapara);

	/**
	* Perform the LL crossover on a given individuals. Result is a single offspring.
	* @param ind1 - parent individual
	* @param ind2 - parent individual
	* @param iGen - current generation index
	* @param fcode - function according to which the process will be conducted
	*/
	individual LL_Crossover( individual &ind1, individual &ind2, int iGen, function & fcode);

	/**
	* Perform the LL mutation on a single individual.
	* @param child - source and output individual
	* @param gapara - ga specific parameter
	* @param fcode - function according to which the process will be conducted
	*/
	void LL_Mutation(individual &child, int iGen, GA_parameters & gapara, function & fcode);

	/**Update the neighbourhood of a given individual
	@param indiv - source individual, according to which the neighbourhood will be updated.
	@param col - collective to be updated
	@param id - index of a subproblem (neighbourhood)
	@param type - update solutions in - neighborhood (1) or whole population (otherwise)
	* @param iGen - current generation index
	*/
	void Problem_Update(individual &indiv, collective & col, int &id, int &type, int iGen);
	/**Update the neighbourhood of a given individual (globally)
	@param indiv - source individual, according to which the neighbourhood will be updated.
	@param col - collective to be updated
	@param id - index of a subproblem (neighbourhood)
	@param type - update solutions in - neighborhood (1) or whole population (otherwise)
	* @param iGen - current generation index
	*/
	void Problem_Update_Global(individual &indiv, collective & col, int &id, int &type, int iGen);
}

#endif // !MOEAD_H
