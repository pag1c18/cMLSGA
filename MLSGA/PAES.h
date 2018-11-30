#pragma once
#ifndef PAES_H
#define PAES_H

#include "Class.h"

/*
Calculate individuals using PAES algorithm
@param Col - address of a given collective
@param mode - mode which will be used: 1 - choose the best individual and replace the worst indiviiduals; 2 - for each individual; 3 - the same as 1 but with weight vectors; 4 - the same as 2 but with weight vectors
*/
std::vector<individual> PAES_Calc(collective & col, int mode, selection<individual> & scode);

namespace PAES
{
	/*
	*Initialise the PAES parameters*
	@param n_obj - number of objectives
	@param pop - population used for initialisation
	*/
	void PAES_Init(short n_obj, population & pop);
	/*
	*Initialise the PAES population*
	@param n_obj - number of objectives
	@param col - collective used for initialisation
	*/
	void PAES_Population_Init(short n_obj, collective & col);
}

#endif // !PAES_H
