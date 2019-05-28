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



/*This code is based on java code by Antonio J. Nebro, Juan J. Durillo [1]
Translated to C++ and modified for the purposes of the MLSGA framework by Przemyslaw A.Grudniewski (2019)*/


/*[1] Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.


Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>

//Originally published: Q. Lin et al., “A Hybrid Evolutionary Immune Algorithm for Multiobjective Optimization Problems,” IEEE Trans. Evol. Comput., vol. 20, no. 5, pp. 711–729, 2016.*/


#define _CRT_SECURE_NO_WARNINGS

#include "HEIA.h"
#include "NSGAII.h"
#include "MOEAD.h"
#include <time.h>
#include "Support_Functions.h"


extern int nfes;
extern time_t elit_t;
extern time_t col_t;


std::vector<std::vector<individual>> clonepopulation; //Separate population for each collective

std::vector<std::vector<individual>> Archive;	//Archive of individuals for each collective

float DE_rate = 0.5f;
short col_ix;	//Index of the collective
short sort_fit_ix;

namespace HEIA
{
	//Create offspring population
	std::vector<individual> Clone_Operation(collective& col);

	//Update DE population
	void DE_Update(std::vector<individual> & DEpop, collective & col);

	//Update SBX population
	void SBX_Update(std::vector<individual> & SBXpop, collective & col);

	//Update the Archive population
	void Archive_Update(std::vector<individual> & merge_pop, int popsize, int clonesize, std::vector<short> fit_indexes);
	
	//Update the crowding distance of current population
	void Crowding_Dist_Update(std::vector<individual> &front, std::vector<short> fit_indexes);

	//Random Permutation
	void Random_Perm(std::vector<int> & perm, int size);

	//Suppress individuals
	void Suppress(std::vector<individual> &pop, std::vector<short> fit_indexes);

	//Sort individuals according to fitness
	static bool Sort_Fit(const individual& c1, const individual& c2);
	//Sort individuals according to crowding distance
	static bool Sort_Crowd(const individual& c1, const individual& c2);
}


/*
Calculate individuals using HEIA algorithm
@param Col - address of a given collective
@param iGen - current generation index
*/
std::vector<individual> HEIA::HEIA_Calc(collective & col, int iGen)
{
	//get the index of the collective (-1 for the ease of operations)
	col_ix = col.Index_Show() - 1;
	int clonesize = col.Size_Show() / 5;


	if (col.Was_Erased())
	{
		//Update and create new archive
		//Archive_Update(col.Indiv_Show(), col.Size_Show(), clonesize);
		HEIA_Pop_Init(col);
		col.Clear_Erased();
	}



	time_t elite_t_temp = clock();				//starting time of HEIA

	//initialise the temporary populations
	std::vector<individual> SBXpop;
	std::vector<individual> DEpop;
	std::vector<individual> offspringpopulation;
	
	//Do HEIA routine
	offspringpopulation = Clone_Operation(col);

	for (int i = 0; i < offspringpopulation.size(); i++)
	{
		if (HEIA_crossover_prob < Random())
			SBXpop.push_back(offspringpopulation[i]);
		else
			DEpop.push_back(offspringpopulation[i]);
	}

	DE_Update(DEpop, col);
	SBX_Update(SBXpop, col);
	
	//merge populations back
	offspringpopulation.clear();
	offspringpopulation = DEpop;
	offspringpopulation.insert(offspringpopulation.end(), SBXpop.begin(), SBXpop.end());

	Archive_Update(offspringpopulation, col.Size_Show(), clonesize, col.fit_index[0]);

	//Check if the size is maintained
	if (offspringpopulation.size() != col.Size_Show())
		abort();

	//copy the offsping population to the collective
	col.Indiv_Set() = offspringpopulation;

	//calculate time
	elit_t += clock() - elite_t_temp;

	//save the inidividuals
	col.save();

	time_t col_t_temp = clock();				//Starting time of collective operations

												//Calculate fitness of the collective
	col.Fitness_Calc();

	//save col time
	col_t += clock() - col_t_temp;

	return Archive[col_ix];
}

/*
Initialise the algorithm specific parameters for whole population
@param n_obj - current number of objectives
@param n_col - current number of collectives
@param pop - current (overall) population
*/
void HEIA::HEIA_Init(short n_col, population & pop)
{
	Archive.clear();
	clonepopulation.clear();

	Archive = std::vector<std::vector<individual>>(n_col, std::vector<individual>());
	clonepopulation = Archive;

	//Calculate fitness of the population
	pop.Fitness_Calc();

	nfes += pop.Size_Show();

}

/*
Initialise the algorithm specific parameters for each collective
@param n_obj - current number of objectives
@param n_col - current number of collectives
@param pop - current (overall) population
*/
void HEIA::HEIA_Pop_Init(collective & col)
{
	time_t elite_t_temp = clock();				//starting time of HEIA

	col_ix = col.Index_Show() - 1;
	//Clear the precious archives
	Archive[col_ix].clear();
	clonepopulation[col_ix].clear();
	
	int clonesize = col.Size_Show() / 5;

	//copy the population for ease of operations
	std::vector<individual> col_pop = col.Indiv_Show();

	//Calculate the ranks and crowding distances
	NSGAII::Rank_Crowding_Distance_Assign(col_pop, col.fit_index[0]);

	int index = 0;
	sort_fit_ix = 0;
	std::sort(col_pop.begin(), col_pop.end(), Sort_Fit);
	//copy the first front to archive
	
		for (int i = 0; i < col_pop.size(); i++)
		{
			if (col_pop[i].Rank_Show() == 1)
			{
				col_pop[i].table.clear();
				col_pop[i].table.push_back(index);
				index++;
				Archive[col_ix].push_back(col_pop[i]);
			}
		}




	std::vector<individual> temp_arch = Archive[col_ix];
	//sort archive according to crowding distance
	std::sort(temp_arch.begin(), temp_arch.end(), Sort_Crowd);
	//Copy to clonepopultion
	for (int i = 0; i < Archive[col_ix].size() && i < clonesize; i++)
	{
		clonepopulation[col_ix].push_back(temp_arch[i]);
	}

	temp_arch.clear();

	col.Indiv_Set() = col_pop;

	//calculate time
	elit_t += clock() - elite_t_temp;
}

//Update the external population for the time step
void HEIA::HEIA_Time_Update(function &fcode)
{
	//Update the external populations
	for (int i_col = 0; i_col < Archive.size(); i_col++)
	{
		//copy the temp parameters
		int arch_size = Archive[i_col].size();
		
		std::vector<individual> temp_vect = Archive[i_col];

		//Update the population
		for (int i_ind = 0; i_ind < arch_size; i_ind++)
		{
			temp_vect[i_ind].Fitness_Calc(fcode);
			nfes++;
			temp_vect[i_ind].save();
		}
		Archive[i_col] = temp_vect;

		
		//copy the temp parameters
		int clone_size = clonepopulation[i_col].size();

		temp_vect = clonepopulation[i_col];

		//Update the population
		for (int i_ind = 0; i_ind < clone_size; i_ind++)
		{
			temp_vect[i_ind].Fitness_Calc(fcode);
			nfes++;
			temp_vect[i_ind].save();
		}
		clonepopulation[i_col] = temp_vect;
	}
}

//Create offspring population
std::vector<individual> HEIA::Clone_Operation(collective& col)
{
	int size = col.Size_Show();
	std::vector<individual> offspring;		//output
	std::vector<individual> parents = clonepopulation[col_ix];

	double min_distance = 0.0;
	double max_distance = 1.0;
	double sum_distance = 0.0;
	int k = 0;
	for (k = 0; k < parents.size(); k++) 
	{
		if (parents[k].Crowd_Dist_Show() != INF) 
		{
			max_distance = parents[k].Crowd_Dist_Show();
			min_distance = parents[parents.size() - 1].Crowd_Dist_Show();
			for (int l = 0; l<k; l++) 
			{
				parents[l].Crowd_Dist_Set(parents[k].Crowd_Dist_Show());
			}
			break;
		}
	} 
	if (parents[0].Crowd_Dist_Show() == INF) 
	{
		for (int l = 0; l<parents.size(); l++) 
		{
			parents[l].Crowd_Dist_Set(1.0);
		}
	}
	for (k = 0; k < parents.size(); k++) 
	{
		sum_distance += parents[k].Crowd_Dist_Show();
	}
	std::vector<double>clones(parents.size(), 0.);
	for (k = 0; k<parents.size(); k++) 
	{
		clones[k] = std::ceil(size*parents[k].Crowd_Dist_Show() / sum_distance);
		if (sum_distance == 0) 
		{
			clones[k] = std::ceil((double)size / parents.size());
		}
		
	}
	int remain = size;
	int i = 0;
	for (k = 0; k<parents.size(); k++) 
	{
		for (int l = 0; l<clones[k]; l++) 
		{
			if (remain>0) 
			{
				offspring.push_back(parents[k]);
				remain--;
			}
			i++;
		}
		if (remain == 0)
			break;
	
	}




	return offspring;
}

void HEIA::DE_Update(std::vector<individual> & DEpop, collective & col)
{
	individual offspring(col.FCode_Show()[0]);
	int clonepop_size = clonepopulation[col_ix].size();

	for (int i = 0; i < DEpop.size(); i++)
	{
		if (Termination_Check(col.FCode_Show()[0].Time_Dep()))
			break;

		std::vector<individual> parents;
		parents.push_back(DEpop[i]);

		if (clonepop_size < 20)
		{
			if (clonepop_size > 1)
			{
				std::vector<int> perm;
				Random_Perm(perm, clonepop_size);
				parents.push_back(clonepopulation[col_ix][perm[0]]);
				parents.push_back(clonepopulation[col_ix][perm[1]]);
			}
			else
			{
				parents.push_back(clonepopulation[col_ix][0]);
				parents.push_back(clonepopulation[col_ix][0]);
			}
		}
		else
		{
			if (0.1 < Random())
			{
				int neighbours = DEpop[i].table[0];
				std::vector<int> perm;
				Random_Perm(perm, 20);
				int selected = perm[0];
				int selected2 = perm[1];
				int arch_size = Archive[col_ix].size();
				if (neighbours < 10)
				{
					parents.push_back(Archive[col_ix][selected2]);
					parents.push_back(Archive[col_ix][selected]);
				}
				else if (neighbours > (arch_size - 10))
				{
					parents.push_back(Archive[col_ix][arch_size - 20 + selected2]);
					parents.push_back(Archive[col_ix][arch_size - 20 + selected]);
				}
				else
				{
					parents.push_back(Archive[col_ix][neighbours - 10 + selected2]);
					parents.push_back(Archive[col_ix][neighbours - 10 + selected]);
				}
			}
			else
			{
				std::vector<int> perm;
				Random_Perm(perm, clonepop_size);
				parents.push_back(clonepopulation[col_ix][perm[0]]);
				parents.push_back(clonepopulation[col_ix][perm[1]]);
			}
		}

		offspring = MOEAD::Diff_Evo_XoverB(parents[0], parents[1], parents[2], DE_rate, col.FCode_Show()[0], col.CCode_Show()[0], col.GAPara_Show()[0]);
		col.Mutation(offspring);
		offspring.Fitness_Calc(col.FCode_Show()[0]);
		nfes++;
		DEpop[i] = offspring;
	}
}

void HEIA::SBX_Update(std::vector<individual> & SBXpop, collective & col)
{
	individual offspring = individual(col.FCode_Show()[0]);

	for (int i = 0; i < SBXpop.size(); i++)
	{
		if (Termination_Check(col.FCode_Show()[0].Time_Dep()))
			break;

		std::vector<individual> parents;
		parents.push_back(SBXpop[i]);
		int random_index = Random_I(0, clonepopulation[col_ix].size() - 1);
		parents.push_back(clonepopulation[col_ix][random_index]);
		offspring = col.Crossover_NSGAII(parents)[0];
		col.Mutation(offspring);
		offspring.Fitness_Calc(col.FCode_Show()[0]);
		nfes++;
		SBXpop[i] = offspring;
	}
}

//Update the Archive[col_ix] population with two given populations
void HEIA::Archive_Update(std::vector<individual> & merge_pop, int popsize, int clonesize, std::vector<short> fit_indexes)
{
	//create the union of solutionset and offsprings
	std::vector<individual> solution_union = Archive[col_ix];
	solution_union.insert(solution_union.end(), merge_pop.begin(), merge_pop.end());

	Suppress(solution_union, fit_indexes);

	//Ranking the union
	//Calculate the ranks and crowding distances
	NSGAII::Rank_Crowding_Distance_Assign(solution_union, fit_indexes);

	Archive[col_ix].clear();
	clonepopulation[col_ix].clear();

	//create temporary front
	std::vector<individual> front;
	for (int i = 0; i < solution_union.size(); i++)
	{
		if (solution_union[i].Rank_Show() == 1)
		{
			front.push_back(solution_union[i]);
		}
	}

	Suppress(front, fit_indexes);
	std::sort(front.begin(), front.end(), Sort_Crowd);

	while (front.size() > popsize)
	{
		front.pop_back();
		Crowding_Dist_Update(front, fit_indexes);
		std::sort(front.begin(), front.end(), Sort_Crowd);
	}

	sort_fit_ix = 0;
	std::sort(front.begin(), front.end(), Sort_Fit);	
	//copy the first front to archive
	for (int i = 0; i < front.size(); i++)
	{
		front[i].table.clear();
		front[i].table.push_back(i);
		Archive[col_ix].push_back(front[i]);
	}

	//sort front according to crowding distance
	std::sort(front.begin(), front.end(), Sort_Crowd);
	for (int i = 0; i < front.size() && i < clonesize; i++)
	{
		clonepopulation[col_ix].push_back(front[i]);
	}
}

//Update the crowding distance of current population
void HEIA::Crowding_Dist_Update(std::vector<individual> &front, std::vector<short> fit_indexes)
{
	int nobj = fit_indexes.size();
	int size = front.size();

	if (size <= 2)
	{
		for (int i = 0; i < size; i++)
			front[i].Crowd_Dist_Set(INF);
		return;
	}

	for (int i = 0; i < size; i++)
		front[i].Crowd_Dist_Set(0.0);


	double objetiveMaxn;
	double objetiveMinn;
	double distance;

	for (int i = 0; i<nobj; i++) 
	{
		short fit_ix = fit_indexes[i] - 1;
		// Sort the population by Obj n     
		sort_fit_ix = fit_ix;
		std::sort(front.begin(), front.end(), Sort_Fit);
		objetiveMinn = front[0].Fitness_Show(fit_ix);
		objetiveMaxn = front[front.size() - 1].Fitness_Show(fit_ix);

		//Set de crowding distance            
		front[0].Crowd_Dist_Set(INF);
		front[size - 1].Crowd_Dist_Set(INF);


		for (int j = 1; j < size - 1; j++) 
		{
			distance = front[j + 1].Fitness_Show(fit_ix) - front[j - 1].Fitness_Show(fit_ix);
			if (objetiveMaxn != objetiveMinn)
				distance = distance / (objetiveMaxn - objetiveMinn);
			distance += front[j].Crowd_Dist_Show();
			front[j].Crowd_Dist_Set(distance);
		} // for
	} // for  

}

void HEIA::Random_Perm(std::vector<int> & perm, int size)
{
	std::vector<int> index;
	std::vector<bool> flag;

	for (int i = 0; i < size; i++)
	{
		index.push_back(i);
		flag.push_back(true);
		perm.push_back(0);
	}
	int num = 0;

	while (num < size)
	{
		int start = Random_I(0, size - 1);
		while (true)
		{
			if (flag[start])
			{
				perm[num] = index[start];
				flag[start] = false;
				num++;
				break;
			}
			if (start == (size - 1))
				start = 0;
			else
				start++;
		}
	}
}

void HEIA::Suppress(std::vector<individual> &pop, std::vector<short> fit_indexes)
{
	int nobj = fit_indexes.size();
	double diff;
	int pop_size = pop.size();

	for (int k = 0; k < pop_size; k++) 
	{
		for (int l = k + 1; l < pop_size; l++) 
		{
			int m = 0;
			for (m = 0; m < nobj; m++) 
			{
				short fit_ix = fit_indexes[m] - 1;
				diff = pop[k].Fitness_Show(fit_ix) - pop[l].Fitness_Show(fit_ix);
				if (diff<0)
					diff = -diff;
				if (diff>0.000001) 
				{
					break;
				}
			}
			if (m == nobj) 
			{
				pop.erase(pop.begin() + l);
				l--;
				pop_size--;
			}
		}
	}
}

//Sort individuals according to fitness
static bool HEIA::Sort_Fit(const individual& c1, const individual& c2)
{
	return c1.Fitness_Show(sort_fit_ix)< c2.Fitness_Show(sort_fit_ix);
};

//Sort individuals according to crowding distance
static bool HEIA::Sort_Crowd(const individual& c1, const individual& c2)
{
	return c1.Crowd_Dist_Show() > c2.Crowd_Dist_Show();
};

