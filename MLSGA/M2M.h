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

M2M header
		Functions specific to M2M algorithm


*/

#pragma once
#ifndef M2M_H
#define M2M_H
#include "Const.h"
#include "Define.h"
#include "Class.h"
/**Functions specific to M2M algorithm*/
namespace M2M
{
	/**
	Evolve the individuals using MOEAD/M2M algorithm
	@param Col - address of a given collective
	@param iGen - current generation index
	*/
	std::vector<individual> M2M_Calc(collective & col, int iGen);
	/**
	Initialise the algorithm specific parameters, for a whole population.
	@param n_obj - current number of objectives
	@param pop - current (overall) population
	*/
	void M2M_Init(short n_obj, population & pop);

	/**
	Initialise the algorithm specific parameters, for a single collective.
	@param n_obj - current number of objectives
	@param col - current collective
	*/
	void M2M_Population_Init(short n_obj, collective & col);


	/**
	Obtain the algorithm specific parameters.
	@param weights - weight vectors to be calculated - Output.
	@param center - vector storing central points - Output.
	@param pop_size - size of the population for which vectors will be generated.
	@param n_obj - current number of objectives
	*/
	void M2M_Get_Params(std::vector<std::vector<std::vector<double>>> &weights, std::vector<std::vector<double>> &center, int pop_size, short n_obj);
}

#endif // !M2M_H

