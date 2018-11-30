#pragma once
#ifndef IBEA_H
#define IBEA_H
#include "Class.h"

namespace IBEA
{
	/*
	Calculate individuals using IBEA algorithm
	@param Col - address of a given collective
	@param iGen - current generation index
	*/
	std::vector<individual> IBEA_Calc(collective & col, int iGen);

	/*
	Initialise the algorithm specific parameters
	@param n_col - current number of collectives
	@param pop - current (overall) population
	*/
	void IBEA_Init(short n_col, population & pop);

	/*
	Initialise the algorithm specific parameters
	@param n_obj - current number of objectives
	@param col - current (overall) population
	*/
	void IBEA_Pop_Init(short n_obj, collective & col);

	//Update the external population for the time step
	void IBEA_Time_Update(function &fcode);
}


#endif // !IBEA_H