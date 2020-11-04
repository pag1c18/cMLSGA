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



/*This code is based on C++ code by Miqing Li[1]
Modified for the purposes of the MLSGA framework by Przemyslaw A.Grudniewski (2019)*/

/*==========================================================================
//  [1] Implementation of the Bi-Criterion Evolution (BCE) framework with MOEA/D
//  Last update 25, MAR, 2016
//
//  Please find details of BCE in the following paper
//  M. Li, S. Yang, and X. Liu. Pareto or non-Pareto: Bi-criterion evolution in multi-objective optimization.
//  IEEE Transactions on Evolutionary Computation, vol. 20, no. 5, pp. 645-665, 2016.
//
//  Please find details of this MOEA/D version in the following paper
//  H. Li and Q. Zhang. Multiobjective optimization problems with complicated Pareto sets, MOEA/D and NSGA-II.
//  IEEE Transactions on Evolutionary Computation, vol. 13, no. 2, pp. 284¨C302, 2009.
//
//  The source code of was implemented by Miqing Li (http://www.cs.bham.ac.uk/~limx)
//
//  The codes are free for reserach work.
//  If you have any problem with the source codes, please contact
//  Miqing Li at limitsing@gmail.com
===========================================================================*/



/*Modified by Przemyslaw Grudniewski for the purposed of MLSGA-framework, 2019*/

#define _CRT_SECURE_NO_WARNINGS

#include "BCE.h"
#include "MOEAD.h"
#include "Support_Functions.h"
#include "Pen_Const.h"
#include <time.h>
std::vector<std::vector<individual>> PC_pop;

extern int nfes;
extern time_t elit_t;
extern time_t col_t;
#define INF 1.0e14
extern int pop_size;

int BCE_PC_capacity;

namespace BCE
{
	void enter_PCpop(std::vector<individual> &PCp, collective &NPCp, std::vector<short>& fit_indexes);
	void update_PCpop(std::vector<individual> &PCp, individual &ind, std::vector<short>& fit_indexes);
	void explore_PCpop(std::vector<individual> &PCp, collective &NPCp,  int iGen, std::vector<short>& fit_indexes);
	void NPC_evolution(collective &NPCp, std::vector<individual> &PCp, int iGen);
	void maintain_PCpop(std::vector<individual> &PCp, std::vector<short>& fit_indexes);
	
	void normalize_bothpop(std::vector<individual> &pop, int PC_popsize, population &NPCp, int NPC_popsize);
	void normalize_pop(std::vector<individual> &tempp, int PC_popsize, std::vector<short>& fit_indexes);
	double determine_radius(std::vector<std::vector<double>> & c, int PC_popsize);
}


std::vector<individual> BCE::BCE_Calc(collective & col, int iGen)
{
	//get the index of the collective (-1 for the ease of operations)
	int ix = col.Index_Show() - 1;

	if (col.Was_Erased())
	{
		time_t elite_t_temp = clock();				//starting time of whole operation
													//create the neighbourhood
		MOEAD::Neighbourhood_Init(col);

		enter_PCpop(PC_pop[ix], col, col.fit_index[0]);

		//set the marker to not erased
		col.Clear_Erased();

		//calculate time
		elit_t += clock() - elite_t_temp;
	}

	BCE_PC_capacity = col.Size_Show();


	time_t elite_t_temp = clock();				//starting time of BCE
	// individual exploration in the PC evolution
	explore_PCpop(PC_pop[ix], col, iGen, col.fit_index[0]);

	if (!(Termination_Check(col.FCode_Show()[0].Time_Dep())))
	// NPC evolution
		NPC_evolution(col, PC_pop[ix], iGen);

	if (PC_pop[ix].size() > BCE_PC_capacity)
	{
		// population maintenance operation in the PC evolution
		maintain_PCpop(PC_pop[ix], col.fit_index[0]);
	}


	//calculate time
	elit_t += clock() - elite_t_temp;

	//save the inidividuals
	//col.save();

	time_t col_t_temp = clock();				//Starting time of collective operations

												//Calculate fitness of the collective
	col.Fitness_Calc();

	//save col time
	col_t += clock() - col_t_temp;

	return PC_pop[ix];
}

void BCE::BCE_Init(short n_obj, short n_col, population & pop)
{
	//clear the vectors
	PC_pop.clear();
	
	//initialise the vector
	PC_pop = std::vector<std::vector<individual>>(n_col, std::vector<individual>{});

	//initialise the MOEAD
	MOEAD::MOEAD_Init(n_obj, pop);
}

void BCE::BCE_Pop_Init(short n_obj, collective & col)
{

	MOEAD::Population_Init(n_obj, col);

	// put the nondominated individuals of the NPC population into the PC population	
	enter_PCpop(PC_pop[col.Index_Show() - 1], col, col.fit_index[0]);
}

void BCE::BCE_Time_Update(function& fcode, short i_col)
{
	//Update the external populations

		//copy the temp parameters
	int size = PC_pop[i_col].size();

	//Update the population
	for (int i_ind = 0; i_ind < size; i_ind++)
	{
		PC_pop[i_col][i_ind].Fitness_Calc(fcode);
		nfes++;
		PC_pop[i_col][i_ind].save();
	}
	if (PENALTY_BASED_CONSTRAINTS)
		Pen_const::Fitness_Recalc(PC_pop[i_col], i_col, false);

}
void BCE::enter_PCpop(std::vector<individual> &PCp, collective &NPCp, std::vector<short>& fit_indexes)
{
	if (PCp.size() == 0)
		PCp.push_back(NPCp.Indiv_Show(0));

	for (int i_ind = 1; i_ind < NPCp.Size_Show(); i_ind++)
	{
		update_PCpop(PCp, NPCp.Indiv_Show(i_ind), fit_indexes);
	}
}

void BCE::update_PCpop(std::vector<individual> &PCp, individual &ind, std::vector<short>& fit_indexes)
{
	int flag;
	for (int i_ind = 0; i_ind < PCp.size(); i_ind++)
	{
		flag = Dominance_Check(ind, PCp[i_ind], fit_indexes);
		switch (flag)
		{
		case 1:
			PCp.erase(PCp.begin() + i_ind);
			break;
		case -1:
			return;
		}
	}

	PCp.push_back(ind);
}

void BCE::explore_PCpop(std::vector<individual> &PCp, collective &NPCp, int iGen, std::vector<short>& fit_indexes)
{

	int i, j, k;
	
	std::vector<int> promising;			// the index of promising individuals in the PC population
	int promising_num;		// the number of promsing individuals in the PC population
	int nobj = NPCp.FCode_Show()[0].Objs();


	double distance, radius;
	int count, original_size;
	int rand, rand2;
	std::vector<individual> child{ NPCp.Indiv_Show(0),NPCp.Indiv_Show(1) };


	// copy PCp into tempp for ease of operation 	
	std::vector<individual>tempp = PCp;
	
	//Copy NPCp for ease of operation
	//std::vector<individual> temp_NPCp = NPCp.Indiv_Show();

	//	normalise both poulations according to the PC individuals 
	normalize_bothpop(tempp, PCp.size(), NPCp, NPCp.Size_Show());

	std::vector<std::vector<double>> c(PCp.size(),std::vector<double>(PCp.size(),0));				// for recording the Euclidean distance between any two individuals in the PC population

	//  calculate the Euclidean distance among individuals
	for (i = 0; i<PCp.size(); i++)
	{
		for (j = i + 1; j<PCp.size(); j++)
		{
			distance = 0;
			for (k = 0; k<nobj; k++)
			{
				distance += (tempp[i].saved_fitness[k] - tempp[j].saved_fitness[k]) * (tempp[i].saved_fitness[k] - tempp[j].saved_fitness[k]);
			}
			distance = sqrt(distance);
			c[i][j] = distance;
			c[j][i] = distance;
		}
		c[i][i] = INF;
	}

	//	calculate the radius for individual exploration
	if (PCp.size() == 1)
	{
		radius = 0;
	}
	else
	{
		radius = determine_radius(c, PCp.size());
	}

	// find the promising individuals for exploration
	promising_num = 0;
	for (i = 0; i<PCp.size(); i++)
	{
		count = 0;					// record how many NPC individuals are located in this PC individual's niche
		for (j = 0; j<NPCp.Size_Show(); j++)
		{
			distance = 0;
			for (k = 0; k<nobj; k++)
			{
				distance += (tempp[i].saved_fitness[k] - NPCp.Indiv_Show(j).saved_fitness[k]) * (tempp[i].saved_fitness[k] - NPCp.Indiv_Show(j).saved_fitness[k]);
			}
			distance = sqrt(distance);
			if (distance <= radius * (double)PCp.size() / (double)BCE_PC_capacity)
			{
				count++;
			}
		}
		if (count < 2)				// when the niche has no NPC individual or has only one NPC individual
		{
			promising.push_back(i);
			promising_num++;
		}
	}

	// explore these promising individuals
	original_size = PCp.size();
	if (promising_num > 0)
	{
		//child = (individual *)malloc(sizeof(individual));
		//allocate_memory_ind(child);

		//child2 = (individual *)malloc(sizeof(individual));
		//allocate_memory_ind(child2);

		for (i = 0; i<promising_num; i++)
		{
			if (original_size > 1)			// when the number of the PC individuals is larger than 1 
			{
				if (BCE_mode == 1)				// SBX crossover
				{
					do
					{
						rand = Random_I(0, original_size - 1);
					} while (promising[i] == rand);
					//create the selected vector
					std::vector<individual> selected{ tempp[promising[i]], tempp[rand] };
					//crossover(&PCp->ind[promising[i]], &PCp->ind[rand], child, child2);
					child = NPCp.Crossover_NSGAII(selected);
				}
				if (BCE_mode == 2)				// DE
				{
					rand = Random_I(0, original_size - 1);
					do
					{
						rand2 = Random_I(0, original_size - 1);
					} while (rand2 == rand);
					if (MOEAD_OVERRIDE)
					{
						short index = 0;
						if (Random() < 0.5)
							index = 1;
						child[0] = NPCp.CCode_Show()[0].Crossover(std::vector<individual>{  tempp[promising[i]], tempp[rand] }, NPCp.GAPara_Show()[0], NPCp.FCode_Show()[0])[index];
					}
					else
						child[0] = MOEAD::Diff_Evo_XoverB(tempp[promising[i]], tempp[rand], tempp[rand2], 0.5, NPCp.FCode_Show()[0], NPCp.CCode_Show()[0], NPCp.GAPara_Show()[0]);
				}
			}
			else
			{
				child[0] = PCp[0];
			}
			//copy dummy values
			child[0].saved_fitness = PCp[0].saved_fitness;



			NPCp.Mutation(child[0]);		// PM mutation

			NPCp.population::Fitness_Calc(child[0]);
			nfes++;

			if (PENALTY_BASED_CONSTRAINTS)
				Pen_const::Fitness_Recalc(child[0],NPCp.Index_Show()-1);

			child[0].save();


			int type = 3;
			MOEAD::Problem_Update(child[0], NPCp, type, type, iGen);			// update the NPC population by "child" (i.e., NPC selection)

			update_PCpop(PCp, child[0], fit_indexes);			// update the PC population by "child"

			

			

			if (Termination_Check(NPCp.FCode_Show()[0].Time_Dep()))
				return;
		}
	}


}

void BCE::NPC_evolution(collective &NPCp, std::vector<individual> &PCp, int iGen)
{

	int i;
	int r1, r2, p1, p2;
	std::vector<individual> child{ NPCp.Indiv_Show(0),NPCp.Indiv_Show(1) };
	int type;


	//Copy NPCp for ease of operation
	//std::vector<individual> temp_NPCp = NPCp.Indiv_Show();

	for (i = 0; i < NPCp.Size_Show(); i++)
	{
		if (Random() < MOEAD_mating_chance)
			type = 1;
		else
			type = 2;

		if (type == 1)
		{
			int niche = NPCp.Indiv_Show(i).table.size();
			r1 = Random_I(0, niche - 1);
			do {
				r2 = Random_I(0, niche - 1);
			} while (r2 == r1);

			p1 = NPCp.Indiv_Show(i).table[r1];
			p2 = NPCp.Indiv_Show(i).table[r2];
		}
		else
		{
			p1 = Random_I(0, NPCp.Size_Show() - 1);
			do
			{
				p2 = Random_I(0, NPCp.Size_Show() - 1);
			} while (p2 == p1);
		}

		if (BCE_mode == 1)				// SBX crossover
		{
			if (p1 != i)
			{
				//create the selected vector
				std::vector<individual> selected{ NPCp.Indiv_Show(i), NPCp.Indiv_Show(p1) };
				//crossover(&NPC_pop->ind[i], &NPC_pop->ind[p1], child, child2);
				child = NPCp.Crossover_NSGAII(selected);
			}
			else
			{
				//create the selected vector
				std::vector<individual> selected{ NPCp.Indiv_Show(i), NPCp.Indiv_Show(p2) };
				//crossover(&NPC_pop->ind[i], &NPC_pop->ind[p2], child, child2);
				child = NPCp.Crossover_NSGAII(selected);
			}
				
		}
		if (BCE_mode == 2)				// DE
		{
			//DE(&NPC_pop->ind[i], &NPC_pop->ind[p1], &NPC_pop->ind[p2], child);
			if (MOEAD_OVERRIDE)
			{
				short index = 0;
				if (Random() < 0.5)
					index = 1;
				child[0] = NPCp.CCode_Show()[0].Crossover(std::vector<individual>{ NPCp.Indiv_Show(p1), NPCp.Indiv_Show(p2) }, NPCp.GAPara_Show()[0], NPCp.FCode_Show()[0])[index];
			}
			else
				child[0] = MOEAD::Diff_Evo_XoverB(NPCp.Indiv_Show(i), NPCp.Indiv_Show(p1), NPCp.Indiv_Show(p2), 0.5, NPCp.FCode_Show()[0], NPCp.CCode_Show()[0], NPCp.GAPara_Show()[0]);
		}
		//copy dummy values
		child[0].saved_fitness = NPCp.Indiv_Show(i).saved_fitness;
		
		NPCp.Mutation(child[0]);		// PM mutation

		NPCp.population::Fitness_Calc(child[0]);
		nfes++;

		if (PENALTY_BASED_CONSTRAINTS)
			Pen_const::Fitness_Recalc(child[0], NPCp.Index_Show() - 1);

		child[0].save();
		
		MOEAD::Problem_Update(child[0], NPCp, i, type, iGen);			// update the NPC population by "child" (i.e., NPC selection)	

		update_PCpop(PCp, child[0], NPCp.fit_index[0]);			// update the PC population by "child"

		

		if (Termination_Check(NPCp.FCode_Show()[0].Time_Dep()))
			break;
	}
}

void BCE::maintain_PCpop(std::vector<individual> &PCp, std::vector<short>& fit_indexes)
{

	int i, j, k;
	int current_size, original_size;
	int index;
	double radius, distance, max_crowding;
	int nobj = fit_indexes.size();

	
	

	//	normalise the PC poulation 
	normalize_pop(PCp, PCp.size(), fit_indexes);

	std::vector<std::vector<double>> c(PCp.size(), std::vector<double>(PCp.size(), 0));				// for recording the Euclidean distance between any two individuals in the PC population

	for (i = 0; i<PCp.size(); i++)
	{
		for (j = i + 1; j<PCp.size(); j++)
		{
			distance = 0;
			for (k = 0; k<nobj; k++)
			{
				short fit_ix = fit_indexes[k] - 1;
				distance += (PCp[i].saved_fitness[fit_ix] - PCp[j].saved_fitness[fit_ix]) * (PCp[i].saved_fitness[fit_ix] - PCp[j].saved_fitness[fit_ix]);
				//distance += (PCp->ind[i].obj_norm[fit_ix] - PCp->ind[j].obj_norm[fit_ix]) * (PCp->ind[i].obj_norm[fit_ix] - PCp->ind[j].obj_norm[fit_ix]);
			}
			distance = sqrt(distance);
			c[i][j] = distance;
			c[j][i] = distance;
		}
		c[i][i] = INF;
	}

	//	calculate the radius for population maintenance
	radius = determine_radius(c, PCp.size());

	//	use "mark" to record which idividuals should be removed
	for (i = 0; i<PCp.size(); i++)
	{
		PCp[i].Rank_Set(1);								// "1" means that the individual is in the current PC population
		PCp[i].Crowd_Dist_Set(1);				// initialisation of PC individuals' crowding degree
		/*temp = PCp->ind[i].nicheNeighbor->child;	// initialisation of PC individuals' neighbor in their niche
		while (temp != NULL)
		{
			temp = del(temp);
			temp = temp->child;
		}*/
		PCp[i].table.clear();
	}

	// find neighbors and calculate the crowding degree
	for (i = 0; i<PCp.size(); i++)
	{
		for (j = i + 1; j<PCp.size(); j++)
		{
			if (c[i][j] < radius)
			{
				PCp[i].Crowd_Dist_Set(PCp[i].Crowd_Dist_Show() * c[i][j] / radius);
				PCp[j].Crowd_Dist_Set(PCp[j].Crowd_Dist_Show() * c[i][j] / radius);

				PCp[i].table.push_back(j);
				PCp[j].table.push_back(i);
				//insert(PCp->ind[i].nicheNeighbor, j);
				//insert(PCp->ind[j].nicheNeighbor, i);
			}
		}
	}
	for (i = 0; i<PCp.size(); i++)
	{
		PCp[i].Crowd_Dist_Set(1.0 - PCp[i].Crowd_Dist_Show());
		//PCp->ind[i].crowd_degree = 1.0 - PCp->ind[i].crowd_degree;
	}


	current_size = PCp.size();
	do
	{
		// find the individual with the highest crowding degree in the current PC population
		max_crowding = -1;
		for (i = 0; i<PCp.size(); i++)
		{
			if (PCp[i].Rank_Show() == 1)
			{
				if (PCp[i].Crowd_Dist_Show() > max_crowding)
				{
					max_crowding = PCp[i].Crowd_Dist_Show();
					index = i;
				}
			}
		}

		if (max_crowding == 0)			// this means that all the remaining individuals are not neighboring to each other
		{
			// in this case, randomly remove some until the PC size reduces to the capacity
			while (current_size > BCE_PC_capacity)
			{
				do
				{
					index = Random_I(0, PCp.size() - 1);
				} while (PCp[index].Rank_Show() == 0);
				//mark[index] = 0;
				PCp[index].Rank_Set(0);
				current_size--;
			}
		}
		else
		{
			PCp[index].Rank_Set(0);					// "0" means that individual "index" is removed from the PC population	
												// renew the information of the neighbors of individual "index"; 
												// this includes removing individual "index" from their neighbor list and adjusting their crowding degree  
			for (int ix = 0; ix < PCp[index].table.size(); ix++)
			{
				int n_index = PCp[index].table[ix];
				// remove individual "index" from the neighbor list
				for (int n_ix = 0; n_ix < PCp[n_index].table.size(); n_ix++)
				{
					if (PCp[n_index].table[n_ix] == index)
					{
						PCp[n_index].table.erase(PCp[n_index].table.begin() + n_ix);
						break;
					}
				}

				// adjust the crowding degree for two situations, i.e. whether the two individuals are overlapping
				if (c[index][n_index] != 0)
				{
					PCp[n_index].Crowd_Dist_Set(1.0 - (1.0 - PCp[n_index].Crowd_Dist_Show()) / (c[index][n_index] / radius));
				}
				else
				{
					PCp[n_index].Crowd_Dist_Set(1.0);
					for (int n_ix = 0; n_ix < PCp[n_index].table.size(); n_ix++)
					{
						int n_index2 = PCp[n_index].table[n_ix];
						PCp[n_index].Crowd_Dist_Set(PCp[n_index].Crowd_Dist_Show()*c[n_index2][n_index] / radius);
					}
					PCp[n_index].Crowd_Dist_Set(1.0 - PCp[n_index].Crowd_Dist_Show());
				}
			}
			current_size--;

			/*temp = PCp->ind[index].nicheNeighbor->child;
			while (temp != NULL)
			{
				// remove individual "index" from the neighbor list
				temp2 = findnode(PCp->ind[temp->index].nicheNeighbor->child, index);
				temp2 = del(temp2);

				// adjust the crowding degree for two situations, i.e. whether the two individuals are overlapping
				if (c[index][temp->index] != 0)
				{
					PCp->ind[temp->index].crowd_degree = 1.0 - (1.0 - PCp->ind[temp->index].crowd_degree) / (c[index][temp->index] / radius);
				}
				else
				{
					PCp->ind[temp->index].crowd_degree = 1.0;
					temp2 = PCp->ind[temp->index].nicheNeighbor->child;
					while (temp2 != NULL)
					{
						PCp->ind[temp->index].crowd_degree *= c[temp2->index][temp->index] / radius;
						temp2 = temp2->child;
					}
					PCp->ind[temp->index].crowd_degree = 1.0 - PCp->ind[temp->index].crowd_degree;
				}
				temp = temp->child;
			}*/
		}
	} while (current_size > BCE_PC_capacity);

	// renew the PC population
	original_size = PCp.size();
	for (i = 0; i < original_size; i++)
	{
		if (PCp[i].Rank_Show() == 0)
		{
			PCp.erase(PCp.begin() + i);
			i--;
			original_size--;
		}
	}

	return;
}

void BCE::normalize_bothpop(std::vector<individual> &PCp, int PC_popsize, population &NPCp, int NPC_popsize)
{
	int i, j;
	std::vector<short> fit_indexes = NPCp.fit_index[0];
	int nobj = fit_indexes.size();

	std::vector<double> max(nobj, INF);
	std::vector<double> min(nobj, -INF);

	for (j = 0; j<nobj; j++)
	{
		short fit_ix = fit_indexes[j] - 1;
		for (i = 0; i<PC_popsize; i++)
		{
			if (PCp[i].Fitness_Show(fit_ix)<min[j])
			{
				min[j] = PCp[i].Fitness_Show(fit_ix);
			}
			if (PCp[i].Fitness_Show(fit_ix)>max[j])
			{
				max[j] = PCp[i].Fitness_Show(fit_ix);
			}
		}

		if (max[j] == min[j])
		{
			for (i = 0; i<PC_popsize; i++)
			{
				PCp[i].saved_fitness[fit_ix] = 0.0;
			}
			for (i = 0; i<NPC_popsize; i++)
			{
				NPCp.Indiv_Set(i).saved_fitness[fit_ix] = NPCp.Indiv_Show(i).Fitness_Show(fit_ix) - min[j];
			}
		}
		else
		{
			for (i = 0; i<PC_popsize; i++)
			{
				PCp[i].saved_fitness[fit_ix] = (PCp[i].Fitness_Show(fit_ix) - min[j]) / (max[j] - min[j]);
			}
			for (i = 0; i<NPC_popsize; i++)
			{
				NPCp.Indiv_Set(i).saved_fitness[fit_ix] = (NPCp.Indiv_Show(i).Fitness_Show(fit_ix) - min[j]) / (max[j] - min[j]);
			}
		}
	}
}

void BCE::normalize_pop(std::vector<individual> &PCp, int PC_popsize, std::vector<short>& fit_indexes)
{
	int i, j;

	int nobj = fit_indexes.size();
	std::vector<double> max(nobj, INF);
	std::vector<double> min(nobj, -INF);

	for (j = 0; j<nobj; j++)
	{
		short fit_ix = fit_indexes[j] - 1;
		min[j] = INF;
		max[j] = -INF;
		for (i = 0; i<PC_popsize; i++)
		{
			if (PCp[i].Fitness_Show(fit_ix)<min[j])
			{
				min[j] = PCp[i].Fitness_Show(fit_ix);
			}
			if (PCp[i].Fitness_Show(fit_ix)>max[j])
			{
				max[j] = PCp[i].Fitness_Show(fit_ix);
			}
		}

		if (max[j] == min[j])
		{
			for (i = 0; i<PC_popsize; i++)
			{
				PCp[i].saved_fitness[j] = 0.0;
			}
		}
		else
		{
			for (i = 0; i<PC_popsize; i++)
			{
				PCp[i].saved_fitness[fit_ix] = (PCp[i].Fitness_Show(fit_ix) - min[j]) / (max[j] - min[j]);
			}
		}
	}
	return;
}

double BCE::determine_radius(std::vector<std::vector<double>> & c, int PC_popsize)
{
	int i, j, l, m;
	double sum;
	int k_closest = 3;				// set k = 3
	std::vector<double> ave_dist(k_closest, 0);				// record the average distance of an individual to its kth closest individual in the population
	double radius;
	std::vector<std::vector<double>> closest_dist(PC_popsize, std::vector<double>(k_closest, 0));			// record the distance of the individual to its kth closest individual 

	for (i = 0; i<PC_popsize; i++)
	{
		for (j = 0; j<k_closest; j++)
		{
			closest_dist[i][j] = INF;
		}
	}

	for (i = 0; i<PC_popsize; i++)
	{
		for (j = i + 1; j<PC_popsize; j++)
		{
			for (l = 0; l<k_closest; l++)
			{
				if (closest_dist[i][l] > c[i][j])
				{
					for (m = k_closest - 1; m>l; m--)
					{
						closest_dist[i][m] = closest_dist[i][m - 1];
					}
					closest_dist[i][l] = c[i][j];
					break;
				}
			}
			for (l = 0; l<k_closest; l++)
			{
				if (closest_dist[j][l] > c[i][j])
				{
					for (m = k_closest - 1; m>l; m--)
					{
						closest_dist[j][m] = closest_dist[j][m - 1];
					}
					closest_dist[j][l] = c[i][j];
					break;
				}
			}
		}
	}

	for (j = 0; j<k_closest; j++)
	{
		sum = 0;
		for (i = 0; i<PC_popsize; i++)
		{
			sum += closest_dist[i][j];
		}
		ave_dist[j] = sum / PC_popsize;
	}

	j = k_closest - 1;
	do
	{
		radius = ave_dist[j];
		j--;
	} while (radius == INF);

	return radius;
}