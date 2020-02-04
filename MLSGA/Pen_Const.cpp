/*
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

#include "Pen_Const.h"

namespace Pen_const {
	short n_cons;
	short n_obj;
	std::vector<std::vector<double>> c_max;
	std::vector<double> min_v ;
	std::vector<double> max_v;
	std::vector<double> r_f;
	std::vector<int> pop_size;
}

void Pen_const::Initialise(short ncons, short n_col, short nobj)
{
	//get number of constraints
	n_cons = ncons;

	//return for unconstrained problems
	if (n_cons == 0)
		return;

	//get number of objectives
	n_obj = nobj;

	c_max = std::vector<std::vector<double>>(n_col,std::vector<double>(n_cons, 0));
	min_v = std::vector<double>(n_col,INF);
	max_v = std::vector<double>(n_col, 0);
	r_f = std::vector<double>(n_col, 0);
	pop_size =  std::vector<int>(n_col, 0);

}
void Pen_const::Initialise_Col(int psize, short ix)
{
	if (n_cons == 0)
		return;

	pop_size[ix] = psize;
}


void Pen_const::Fitness_Recalc(individual& ind, short ix)
{
	//return for unconstrained problems
	if (n_cons == 0)
		return;

	if (!ind.Cons_Viol_Show())
		return;

	std::vector<double> c(n_cons, 0);
	//Calculate c_i(x)


	//copy the constraint vector
	std::vector<double> const_vect = ind.Cons_Show();

	for (short j = 0; j < n_cons; j++)
	{
		if (const_vect[j] < 0)
			c[j] = abs(const_vect[j]);
	}

	//find the max for every constraint
	for (short j = 0; j < n_cons; j++)
	{
		if (c[j] > c_max[ix][j])
			c_max[ix][j] = c[j];
	}
	//normalise c_i(x)


	for (short j = 0; j < n_cons; j++)
	{
		c[j] /= c_max[ix][j];
	}

	//calculate v(x)
	double v = 0;;


	for (short j = 0; j < n_cons; j++)
	{
		v += c[j];
	}
	v /= n_cons;

	//find min/max v

	if (v < min_v[ix])
		min_v[ix] = v;
	else if (v > max_v[ix])
		max_v[ix] = v;

	//recalculate v'(x)

	v = (v - min_v[ix]) / (max_v[ix] - min_v[ix]);



	//calculate penalty and add to fitness
	//store fitness
	std::vector<double> fitness = ind.Fitness_Show();
	for (short j = 0; j < n_obj; j++)
	{
		double penalty = v * (1 - fitness[j]);
		if ((1 - fitness[j]) > (1 - v))
			penalty *= r_f[ix];
		else
			penalty += 1;
		fitness[j] += penalty;
	}
	//set the new fitness value
	ind.Fitness_Set(fitness);

}
void Pen_const::R_f_update(bool val, short ix)
{
	//return for unconstrained problems
	if (n_cons == 0)
		return;

	if (val)
		r_f[ix] -= 1. / pop_size[ix];
	else
		r_f[ix] += 1. / pop_size[ix];
}

void Pen_const::Fitness_Recalc(std::vector<individual>& ind, short ix, bool u)
{


	//return for unconstrained problems
	if (n_cons == 0)
		return;

	int pop_size_temp = ind.size();

	std::vector <std::vector<double>> c(n_cons, std::vector<double>(pop_size_temp, 0));
	//Calculate c_i(x) for every individual
	for (int i = 0; i < pop_size_temp; i++)
	{

		//copy the constraint vector
		std::vector<double> const_vect = ind[i].Cons_Show();

		for (short j = 0; j < n_cons; j++)
		{
			if (const_vect[j] < 0)
				c[j][i] = abs(const_vect[j]);
		}
	}
	//find the max for every constraint
	if (u)
	{
		for (short j = 0; j < n_cons; j++)
		{
			double c_max_temp = *std::max_element(c[j].begin(), c[j].end());
			if (c_max_temp > c_max[ix][j])
				c_max[ix][j] = c_max_temp;
		}
	}
	//normalise c_i(x)
	for (int i = 0; i < pop_size_temp; i++)
	{

		for (short j = 0; j < n_cons; j++)
		{
			c[j][i] /= c_max[ix][j];
		}
	}
	//calculate v(x)
	std::vector<double> v(pop_size_temp, 0);
	for (int i = 0; i < pop_size_temp; i++)
	{

		for (short j = 0; j < n_cons; j++)
		{
			v[i] += c[j][i];
		}
		v[i] /= n_cons;
	}
	//find min/max v
	if (u)
	{
		auto minmax_v_temp = std::minmax_element(v.begin(), v.end());
		if (min_v[ix] > * minmax_v_temp.first)
			min_v[ix] = *minmax_v_temp.first;
		if (max_v[ix] < *minmax_v_temp.second)
			max_v[ix] = *minmax_v_temp.second;
	}
	//recalculate v'(x)
	for (int i = 0; i < pop_size_temp; i++)
	{
		v[i] = (v[i] - min_v[ix]) / (max_v[ix] - min_v[ix]);
	}
	if (u)
	{
		//calculate r_f
		r_f[ix] = 0;
		for (int i = 0; i < pop_size[ix]; i++)
		{
			if (!ind[i].Cons_Viol_Show())
				r_f[ix]++;
		}
		r_f[ix] /= pop_size[ix];
	}
	//calculate penalty and add to fitness
	for (int i = 0; i < pop_size_temp; i++)
	{
		if (!ind[i].Cons_Viol_Show())
			continue;

		//store fitness
		std::vector<double> fitness = ind[i].Fitness_Show();
		for (short j = 0; j < n_obj; j++)
		{
			double penalty = v[i] * (1 - fitness[j]);
			if ((1 - fitness[j]) > (1 - v[i]))
				penalty *= r_f[ix];
			else
				penalty += 1;
			fitness[j] += penalty;
		}
		//set the new fitness value
		ind[i].Fitness_Set(fitness);
	}
}
void Pen_const::Fitness_Recalc(collective& col, bool u)
{
	Fitness_Recalc(col.Indiv_Set(), col.Index_Show() - 1, u);
}