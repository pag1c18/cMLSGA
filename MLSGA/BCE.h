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

#pragma once
#ifndef BCE_H
#define BCE_H
#include "Class.h"

namespace BCE
{
	/*
	Calculate individuals using BCE algorithm
	@param Col - address of a given collective
	@param iGen - current generation index
	*/
	std::vector<individual> BCE_Calc(collective & col, int iGen);

	/*
	Initialise the algorithm specific parameters
	@param n_obj - current number of objectives
	@param n_col - current number of collectives
	@param pop - current (overall) population
	*/
	void BCE_Init(short n_obj, short n_col, population & pop);

	/*
	Initialise the algorithm specific parameters
	@param n_obj - current number of objectives
	@param n_col - current number of collectives
	@param pop - current (overall) population
	*/
	void BCE_Pop_Init(short n_obj, collective & pop);

	//Update the external population for the time step
	void BCE_Time_Update(function &fcode);
}


#endif // !BCE_H