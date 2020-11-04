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
along with this program.If not, see < https://www.gnu.org/licenses/>. */


#pragma once
#ifndef NSGAII_H
#define NSGAII_H
#include "Class.h"

namespace NSGAII {
	/**
	Evolve the population using NSGAII algorithm
	@param Col - address of a given collective
	*/
	void NSGAII_Calc(collective & col);



	/**
	Calculate the crowding distance for a given population.  Routine to compute crowding distance based on objective function values when the population in in the form of an array
	@param pop - individuals vector for which the crowding will be calculated
	@param c1 - index of the first individual
	@param c2 - index of the last individual
	*/
	void Crowding_Distance_Indices_Assign(std::vector<individual> &pop, int c1, int c2);


	/**
	Calculate the crowding distance for a given population and assign ranks to each individual
	@param now_pop - individuals vector for which the crowding will be calculated
	@param fit_indexes - indexes of objectives which will be used for rank calculation
	*/
	void Rank_Crowding_Distance_Assign(std::vector<individual> & new_pop, std::vector<short>& fit_indexes);


	/**
	Routine to select the final offspring population. The final population is taken from mixed population (old and new one) basing on a individuals ranks.
	@param mixed_pop - vector containing old and new populations.
	@param now_pop - final offspring population (output).
	@param fit_indexes - indexes of objectives which will be used for rank calculation
	*/
	void Nondominated_Sort_Fill(std::vector<individual> & mixed_pop, std::vector<individual> & new_pop, std::vector<short>& fit_indexes);


	/**
	Routine to select the individuals for crossover. Returns the vector of selected individuals.
	@param old_pop - parents population
	@param fit_indexes - indexes of objectives which will be used for selection
	*/
	std::vector<individual> Select(std::vector<individual> & old_pop, std::vector<short> &fit_indexes);

	struct list
	{
		int index;
		list *parent;
		list *child;
	};

	void Insert(list *node, int x);
	list* Del(list *node);
}

#endif // !NSGAII_H
