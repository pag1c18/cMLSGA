#pragma once

#ifndef HEIA_H
#define HEIA_H
#include "Class.h"

namespace HEIA
{
	/*
	Calculate individuals using HEIA algorithm
	@param Col - address of a given collective
	@param iGen - current generation index
	*/
	std::vector<individual> HEIA_Calc(collective & col, int iGen);

	/*
	Initialise the algorithm specific parameters
	@param n_col - current number of collectives
	@param pop - current (overall) population
	*/
	void HEIA_Init( short n_col, population & pop);

	/*
	Initialise the algorithm specific parameters
	@param n_obj - current number of objectives
	@param n_col - current number of collectives
	@param pop - current (overall) population
	*/
	void HEIA_Pop_Init(collective & pop);

	//Update the external population for the time step
	void HEIA_Time_Update(function &fcode);
}




#endif //HEIA_H

