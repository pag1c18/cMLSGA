#pragma once
#ifndef MTS_H
#define MTS_H

#include "Class.h"

namespace MTS
{
	/*
	Set the MTS initial parameters
	@param pop - current population
	@param PF - PF address
	@param ncom - number of collectives
	*/
	void MTS_Init(population & pop, pareto_front & PF, short ncol);
	/*
	Set the MTS initial parameters - collective specific
	@param col - current collective
	*/
	void MTS_Init_Col(collective & col);
	/*
	Calculate individuals using MTS algorithm
	@param Col - address of a given collective
	@param iGen - current generation index
	*/
	void MTS_Calc(collective & col);
}

#endif // !MTS_H