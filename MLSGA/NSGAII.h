/*CODE BASED ON THE ORIGINAL NSGAII CODE FROM THE WEBSITE*/

#pragma once
#ifndef NSGAII_H
#define NSGAII_H
#include "Class.h"

namespace NSGAII {
	/*
	Calculate individuals using NSGAII algorithm
	@param Col - address of a given collective
	@param iGen - current generation
	*/
	void NSGAII_Calc(collective & col, int iGen);
	void Crowding_Distance_Indices_Assign(std::vector<individual> &pop, int c1, int c2);
	void Rank_Crowding_Distance_Assign(std::vector<individual> & new_pop, std::vector<short>& fit_indexes);
}

#endif // !NSGAII_H
