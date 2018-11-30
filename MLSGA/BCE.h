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