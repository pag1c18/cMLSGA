/**
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
along with this program.If not, see < https://www.gnu.org/licenses/>. 

     MTS header

Containing the functions typical for MTS algorithm

*/


#pragma once
#ifndef MTS_H
#define MTS_H

#include "Class.h"

/**Containing the functions typical for MTS algorithm*/
namespace MTS
{
	/**
	Set the MTS initial parameters
	@param pop - current population
	@param PF - PF address
	@param ncom - number of collectives
	*/
	void MTS_Init(population & pop, pareto_front & PF, short ncol);
	/**
	Set the MTS initial parameters - collective specific
	@param col - current collective
	*/
	void MTS_Init_Col(collective & col);
	/**
	Evolve the population using MTS algorithm
	@param Col - address of a given collective
	@param iGen - current generation index
	*/
	void MTS_Calc(collective & col);
}

#endif // !MTS_H