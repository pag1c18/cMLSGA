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
#ifndef MOEAD_H
#define MOEAD_H
#include "Class.h"


namespace MOEAD
{
	/*
	Calculate individuals using MOEAD algorithm
	@param Col - address of a given collective
	@param iGen - current generation index
	*/
	std::vector<individual> MOEAD_Calc(collective & col, int iGen);

	void MOEAD_Init(short n_obj, population & pop);
	void Population_Init(short n_obj, collective & col);
	void Neighbourhood_Init(collective & col);
	void Utility_Comp(collective & col, int iGen);

	individual Diff_Evo_XoverB(individual &ind0, individual &ind1, individual &ind2, double rate, function & fcode, crossover<individual> & ccode, GA_parameters & gapara);
	individual LL_Crossover( individual &ind1, individual &ind2, int iGen, function & fcode);
	void LL_Mutation(individual &child, int iGen, GA_parameters & gapara, function & fcode);


	void Problem_Update(individual &indiv, collective & col, int &id, int &type, int iGen);
	void Problem_Update_Global(individual &indiv, collective & col, int &id, int &type, int iGen);
}

#endif // !MOEAD_H
