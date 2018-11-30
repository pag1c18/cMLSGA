#pragma once
#ifndef M2M_H
#define M2M_H
#include "Const.h"
#include "Define.h"
#include "Class.h"

namespace M2M
{
	/*
	Calculate individuals using MOEAD algorithm
	@param Col - address of a given collective
	@param iGen - current generation index
	*/
	std::vector<individual> M2M_Calc(collective & col, int iGen);

	void M2M_Init(short n_obj, population & pop);

	void M2M_Population_Init(short n_obj, collective & col);

	void M2M_Get_Params(std::vector<std::vector<std::vector<double>>> &weights, std::vector<std::vector<double>> &center, int pop_size, short n_obj);
}

#endif // !M2M_H

