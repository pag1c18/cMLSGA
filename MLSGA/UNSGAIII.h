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

UNSGAIII header
Functions specific for the UNSGAIII algorithm

*/
#pragma once
#ifndef UNSGAIII_H
#define UNSGAIII_H
#include "Class.h"

/**Functions specific for the UNSGAIII algorithm*/
namespace UNSGAIII {
	/**
	Evolve the individuals using UNSGAIII algorithm
	@param Col - address of a given collective
	@param iGen - current generation
	*/
	void UNSGAIII_Calc(collective & col);


	/**
	Initialise the algorithm specific parameters, for a whole population.
	@param n_obj - current number of objectives
	@param n_col - current number of collectives
	@param pop - current (overall) population
	*/
	void Init(short n_obj, short n_col, population & pop);
}

#endif // !UNSGAII_H
