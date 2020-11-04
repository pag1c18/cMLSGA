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

		BCE header
		Functions specific to BCE algorithm


*/

#pragma once
#ifndef BCE_H
#define BCE_H
#include "Class.h"
/**
Functions specific to BCE algorithm
*/
namespace BCE
{
	/**
	Evolve the individuals using BCE algorithm
	@param Col - address of a given collective. For which the evolution will occur.
	@param iGen - current generation index
	*/
	std::vector<individual> BCE_Calc(collective & col, int iGen);

	/**
	Initialise the algorithm specific parameters, for a whole population.
	@param n_obj - current number of objectives
	@param n_col - current number of collectives
	@param pop - current (overall) population
	*/
	void BCE_Init(short n_obj, short n_col, population & pop);

	/**
	Initialise the algorithm specific parameters, for a single collective.
	@param n_obj - current number of objectives
	@param n_col - current number of collectives
	@param pop - current collective
	*/
	void BCE_Pop_Init(short n_obj, collective & pop);

	/**Update the external population for the next time step. Occurs after dynamic change of the problem.
	@param fcode - address of the currently evaluated function
	@param iCol - index of the collective. For which the parameters will be updated.
	*/
	void BCE_Time_Update(function &fcode, short iCol);
}


#endif // !BCE_H