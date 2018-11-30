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
