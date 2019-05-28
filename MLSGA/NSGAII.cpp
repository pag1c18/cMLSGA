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



/*This code is based on C code by  by the Kanpur Genetic Algorithms Laboratory (https://www.iitk.ac.in/kangal/codes.shtml) [1]
Translated to C++ and modified for the purposes of the MLSGA framework by Przemyslaw A.Grudniewski (2019)*/


/*[1] Copyright (c) 2011  Kanpur Genetic Algorithms Laboratory*/



#include "NSGAII.h"
#include "Define.h"
#include "Const.h"
#include "Random.h"
#include "Support_Functions.h"
#include <ctime>

//NSGAII have to give output as a vector of individuals or population (vector of individuals should be easier)

//!!!!! CHeck if changing saving to after elitism will give different results (the same with PF checking)
//!!!! write (if was erased) and copy original parent function



extern time_t selec_t;										//Time of Selection
extern time_t mut_t;					//Time of mutation
extern time_t col_t;					//time of collective
extern time_t elit_t;					//Time of NSGAII
//extern short MLSt;						//MLSt type

namespace NSGAII {



	


	individual Tournament(individual & ind1, individual & ind2, std::vector<short> & fit_indexes);
	void Crowding_Fill(std::vector<individual> & mixed_pop, std::vector<individual> & new_pop, int count, int front_size, list *elite, int pop_size);
	
	void Crowding_Distance_List_Assign(std::vector<individual> &pop, list *lst, int front_size);
	void Crowding_Distance_Assign(std::vector<individual> &pop, std::vector<int> &dist, std::vector<std::vector<int>> &obj_array, int front_size);
	void Front_Obj_Quicksort(std::vector<individual> &pop, int objcount, std::vector<int> &obj_array, int obj_array_size);
	void Front_Obj_Quicksort_Actual(std::vector<individual> &pop, int objcount, std::vector<int> &obj_array, int left, int right);
	void Dist_Quicksort(std::vector<individual> &pop, std::vector<int> &dist, int front_size);
	void Dist_Quicksort_Actual(std::vector<individual> &pop, std::vector<int> &dist, int left, int right);
}
/*
*Calculate new collective using NSGAII algorithm*
@param Col - address of a given collective
*/
void NSGAII::NSGAII_Calc(collective & col)
{
	//Copy the parent vector
	std::vector<individual> parent_pop = col.Indiv_Show();		//parent population



	//Check if collective need initial crowding distance calculation and ranking
	if (col.Was_Erased())
	{
		time_t elite_t_temp = clock();				//starting time of whole operation
		
		//calculate crowding distance and ranking
		Rank_Crowding_Distance_Assign(parent_pop, col.fit_index[0]);

		//set flag to has not been erased
		col.Clear_Erased();

		//calculate time
		elit_t += clock() - elite_t_temp;
	}


	//Create the offspring vector
	std::vector<individual> off_pop = col.Crossover_NSGAII(Select(parent_pop, col.fit_index[0])); //offspring population

	time_t mut_t_temp = clock();										//starting time of the mutation
	//Mutate the offspring population and evaluate fitness
	col.Mutation_NSGAII(off_pop);
	//calculate the time of the mutation
	mut_t += clock() - mut_t_temp;

	

	//Add two vectors together
	std::vector<individual> mix_pop = parent_pop;
	mix_pop.insert(mix_pop.end(), off_pop.begin(), off_pop.end());

	time_t elite_t_temp = clock();				//Starting time of nondominated sort
		

	//Sort, rank and fill archive
	Nondominated_Sort_Fill(mix_pop, parent_pop, col.fit_index[0]);
	
	//save the nondiminated sort time
	elit_t += clock() - elite_t_temp;

	time_t col_t_temp = clock();				//Starting time of collective operations
	//Copy new individuals into collective
	col = collective::collective(parent_pop, col);
	//save col time
	col_t += clock() - col_t_temp;

	//save individuals to the file
	col.save();

	col_t_temp = clock();				//Starting time of collective operations
	//Calculate fitness of the collective
	col.Fitness_Calc();

	//save col time
	col_t += clock() - col_t_temp;
}

/*
*Select individuals for crossover*
@param old_pop - address of the old population
*/
std::vector<individual> NSGAII::Select(std::vector<individual> & old_pop, std::vector<short> &fit_indexes)
{
	time_t select_t_temp = clock();										//starting time of the selection
	std::vector<int> a1, a2;				//Value vector
	std::vector<individual> off_pop;			//offspring population - output
	int pop_size = old_pop.size();			//size of the old population
	

	//Assign the values to the value vector
	for (int i = 0; i<pop_size; i++)
	{
		a1.push_back(i);
		a2.push_back(i);
	}

	//increase the size of the value vector to the multiply of 4
	int mult_size = 0;			//How many number have been added
	while ((a1.size() % 4) != 0)
	{
		a1.push_back(Random_I(0, pop_size - 1));
		a2.push_back(Random_I(0, pop_size - 1));
		mult_size++;
	}
	//Mix the values vectors
	int a_size = a1.size();
	for (int i = 0; i < a_size; i++)
	{
		int temp;		//temporary value of a
		int rand;		//random value
		rand = Random_I(i, a_size - 1);
		temp = a1[rand];
		a1[rand] = a1[i];
		a1[i] = temp;
		rand = Random_I(i, a_size - 1);
		temp = a2[rand];
		a2[rand] = a2[i];
		a2[i] = temp;
	}
	//Create the offspring population
	for (int i = 0; i < pop_size; i += 4)
	{  
			off_pop.push_back(Tournament(old_pop[a1[i]], old_pop[a1[i + 1]], fit_indexes));
			off_pop.push_back(Tournament(old_pop[a1[i + 2]], old_pop[a1[i + 3]], fit_indexes));
			off_pop.push_back(Tournament(old_pop[a2[i]], old_pop[a2[i + 1]], fit_indexes));
			off_pop.push_back(Tournament(old_pop[a2[i + 2]], old_pop[a2[i + 3]], fit_indexes));
	}

	//remove the additional individuals
	for (int i = 0; i < mult_size; i++)
		off_pop.pop_back();

	//calculate the time of the selection
	selec_t += clock() - select_t_temp;

	if (off_pop.size() != old_pop.size())
		abort();
	return off_pop;
}

/*
*Do binary Tournament*
@param ind1 - address of 1st individual
@param ind2 - address of 2nd individual
*/
individual NSGAII::Tournament(individual & ind1, individual & ind2, std::vector<short>& fit_indexes)
{
	int flag;
	//use check dominance specific to selected mode
	
	flag = Dominance_Check(ind1, ind2, fit_indexes);
	

	if (flag == 1)
	{
		return (ind1);
	}
	if (flag == -1)
	{
		return (ind2);
	}
	if (ind1.Crowd_Dist_Show() > ind2.Crowd_Dist_Show())
	{
		return(ind1);
	}
	if (ind2.Crowd_Dist_Show() > ind1.Crowd_Dist_Show())
	{
		return(ind2);
	}
	if ((Random()) <= 0.5)
	{
		return(ind1);
	}
	else
	{
		return(ind2);
	}
}

/* Function to assign rank and crowding distance to a population of size pop_size*/
void NSGAII::Rank_Crowding_Distance_Assign(std::vector<individual> & new_pop, std::vector<short>& fit_indexes)
{
	int flag;
	int i;
	int end;
	int front_size;
	int pop_size = new_pop.size();
	int rank = 1;
	list *orig;
	list *cur;
	list *temp1, *temp2;
	orig = (list *)malloc(sizeof(list));
	cur = (list *)malloc(sizeof(list));
	front_size = 0;
	orig->index = -1;
	orig->parent = NULL;
	orig->child = NULL;
	cur->index = -1;
	cur->parent = NULL;
	cur->child = NULL;
	temp1 = orig;
	for (i = 0; i<pop_size; i++)
	{
		Insert(temp1, i);
		temp1 = temp1->child;
	}
	do
	{
		if (orig->child->child == NULL)
		{
			new_pop[orig->child->index].Rank_Set(rank);
			new_pop[orig->child->index].Crowd_Dist_Set(INF);
			break;
		}
		temp1 = orig->child;
		Insert(cur, temp1->index);
		front_size = 1;
		temp2 = cur->child;
		temp1 = Del(temp1);
		temp1 = temp1->child;
		do
		{
			temp2 = cur->child;
			do
			{
				end = 0;
				flag = Dominance_Check(new_pop[temp1->index], new_pop[temp2->index], fit_indexes);
				if (flag == 1)
				{
					Insert(orig, temp2->index);
					temp2 = Del(temp2);
					front_size--;
					temp2 = temp2->child;
				}
				if (flag == 0)
				{
					temp2 = temp2->child;
				}
				if (flag == -1)
				{
					end = 1;
				}
			} while (end != 1 && temp2 != NULL);
			if (flag == 0 || flag == 1)
			{
				Insert(cur, temp1->index);
				front_size++;
				temp1 = Del(temp1);
			}
			temp1 = temp1->child;
		} while (temp1 != NULL);
		temp2 = cur->child;
		do
		{
			new_pop[temp2->index].Rank_Set(rank);
			temp2 = temp2->child;
		} while (temp2 != NULL);
		Crowding_Distance_List_Assign(new_pop, cur->child, front_size);
		temp2 = cur->child;
		do
		{
			temp2 = Del(temp2);
			temp2 = temp2->child;
		} while (cur->child != NULL);
		rank += 1;
	} while (orig->child != NULL);
	free(orig);
	free(cur);
	return;
}

/* Routine to perform non-dominated sorting */
void NSGAII::Nondominated_Sort_Fill(std::vector<individual> & mixed_pop, std::vector<individual> & new_pop, std::vector<short>& fit_indexes)
{
	int flag;
	int i, j;
	int end;
	int front_size;
	int archieve_size;
	int rank = 1;
	int pop_size = new_pop.size();
	list *pool;
	list *elite;
	list *temp1, *temp2;
	pool = (list *)malloc(sizeof(list));
	elite = (list *)malloc(sizeof(list));
	front_size = 0;
	archieve_size = 0;
	pool->index = -1;
	pool->parent = NULL;
	pool->child = NULL;
	elite->index = -1;
	elite->parent = NULL;
	elite->child = NULL;
	temp1 = pool;
	for (i = 0; i<2 * pop_size; i++)
	{
		Insert(temp1, i);
		temp1 = temp1->child;
	}
	i = 0;
	do
	{
		temp1 = pool->child;
		Insert(elite, temp1->index);
		front_size = 1;
		temp2 = elite->child;
		temp1 = Del(temp1);
		temp1 = temp1->child;
		do
		{
			temp2 = elite->child;
			if (temp1 == NULL)
			{
				break;
			}
			do
			{
				end = 0;
				flag = Dominance_Check(mixed_pop[temp1->index], mixed_pop[temp2->index], fit_indexes);
				if (flag == 1)
				{
					Insert(pool, temp2->index);
					temp2 = Del(temp2);
					front_size--;
					temp2 = temp2->child;
				}
				else if (flag == 0)
				{
					temp2 = temp2->child;
				}
				else if (flag == -1)
				{
					end = 1;
				}
			} while (end != 1 && temp2 != NULL);
			if (flag == 0 || flag == 1)
			{
				Insert(elite, temp1->index);
				front_size++;
				temp1 = Del(temp1);
			}
			temp1 = temp1->child;
		} while (temp1 != NULL);
		temp2 = elite->child;
		j = i;
		if ((archieve_size + front_size) <= pop_size)
		{
			do
			{
				new_pop[i] = mixed_pop[temp2->index];
				new_pop[i].Rank_Set(rank);
				archieve_size += 1;
				temp2 = temp2->child;
				i += 1;
			} while (temp2 != NULL);
			Crowding_Distance_Indices_Assign(new_pop, j, i - 1);
			rank += 1;
		}
		else
		{
			Crowding_Fill(mixed_pop, new_pop, i, front_size, elite, pop_size);
			archieve_size = pop_size;
			for (j = i; j<pop_size; j++)
			{
				new_pop[j].Rank_Set(rank);
			}
		}
		temp2 = elite->child;
		do
		{
			temp2 = Del(temp2);
			temp2 = temp2->child;
		} while (elite->child != NULL);
	} while (archieve_size < pop_size);
	while (pool != NULL)
	{
		temp1 = pool;
		pool = pool->child;
		free(temp1);
	}
	while (elite != NULL)
	{
		temp1 = elite;
		elite = elite->child;
		free(temp1);
	}
	return;
}

/* Routine to fill a population with individuals in the decreasing order of crowding distance */
void NSGAII::Crowding_Fill(std::vector<individual> & mixed_pop, std::vector<individual> & new_pop, int count, int front_size, list *elite, int pop_size)
{
	std::vector<int> dist(front_size,0);
	list *temp;
	Crowding_Distance_List_Assign(mixed_pop, elite->child, front_size);
	temp = elite->child;
	for (int j = 0; j<front_size; j++)
	{
		dist[j] = temp->index;
		temp = temp->child;
	}
	Dist_Quicksort(mixed_pop, dist, front_size);
	int i, j;
	for (i = count, j = front_size - 1; i<pop_size; i++, j--)
	{
		new_pop[i] = mixed_pop[dist[j]];
	}
	dist.clear();
	return;
}

/* Insert an element X into the list at location specified by NODE */
void NSGAII::Insert(list *node, int x)
{
	list *temp;
	if (node == NULL)
	{
		printf("\n Error!! asked to enter after a NULL pointer, hence exiting \n");
		exit(1);
	}
	temp = (list *)malloc(sizeof(list));
	temp->index = x;
	temp->child = node->child;
	temp->parent = node;
	if (node->child != NULL)
	{
		node->child->parent = temp;
	}
	node->child = temp;
	return;
}

/* Delete the node NODE from the list */
NSGAII::list* NSGAII::Del(list *node)
{
	list *temp;
	if (node == NULL)
	{
		printf("\n Error!! asked to delete a NULL pointer, hence exiting \n");
		exit(1);
	}
	temp = node->parent;
	temp->child = node->child;
	if (temp->child != NULL)
	{
		temp->child->parent = temp;
	}
	free(node);
	return (temp);
}


/* Routine to compute crowding distance based on ojbective function values when the population in in the form of a list */
void NSGAII::Crowding_Distance_List_Assign(std::vector<individual> &pop, list *lst, int front_size)
{
	std::vector<std::vector<int>> obj_array;
	std::vector<int> dist;
	int nobj = pop[0].Fitness_Show().size();
	int i, j;
	list *temp;
	temp = lst;
	if (front_size == 1)
	{
		pop[lst->index].Crowd_Dist_Set(INF);
		return;
	}
	if (front_size == 2)
	{
		pop[lst->index].Crowd_Dist_Set(INF);
		pop[lst->child->index].Crowd_Dist_Set(INF);
		return;
	}
	//obj_array = (int **)malloc(nobj * sizeof(int*));
	//dist = (int *)malloc(front_size * sizeof(int));
	//for (i = 0; i<nobj; i++)
	{
		//obj_array[i] = (int *)malloc(front_size * sizeof(int));
	}
	for (j = 0; j<front_size; j++)
	{
		dist.push_back(temp->index);
		temp = temp->child;
	}
	Crowding_Distance_Assign(pop, dist, obj_array, front_size);
	//free(dist);
	for (i = 0; i<nobj; i++)
	{
		//free(obj_array[i]);
	}
	//free(obj_array);
	return;
}

/* Routine to compute crowding distance based on objective function values when the population in in the form of an array */
void NSGAII::Crowding_Distance_Indices_Assign(std::vector<individual> &pop, int c1, int c2)
{
	std::vector<std::vector<int>> obj_array;
	std::vector<int> dist;
	int i, j;
	int nobj = pop[0].Fitness_Show().size();
	int front_size;
	front_size = c2 - c1 + 1;
	if (front_size == 1)
	{
		pop[c1].Crowd_Dist_Set(INF);
		return;
	}
	if (front_size == 2)
	{
		pop[c1].Crowd_Dist_Set(INF);
		pop[c2].Crowd_Dist_Set(INF);
		return;
	}
	//obj_array = (int **)malloc(nobj * sizeof(int*));
	//dist = (int *)malloc(front_size * sizeof(int));
	//for (i = 0; i<nobj; i++)
	//{
		//obj_array[i] = (int *)malloc(front_size * sizeof(int));
	//}
	for (j = 0; j<front_size; j++)
	{
		dist.push_back(c1++);
	}
	Crowding_Distance_Assign(pop, dist, obj_array, front_size);
	//free(dist);
	for (i = 0; i<nobj; i++)
	{
		//free(obj_array[i]);
	}
	//free(obj_array);
	return;
}

/* Routine to compute crowding distances */
void NSGAII::Crowding_Distance_Assign(std::vector<individual> &pop, std::vector<int> &dist, std::vector<std::vector<int>> &obj_array, int front_size)
{
	int i, j;
	int nobj = pop[0].Fitness_Show().size();
	for (i = 0; i<nobj; i++)
	{
		
		obj_array.push_back(dist);

		Front_Obj_Quicksort(pop, i, obj_array[i], front_size);
	}
	for (j = 0; j<front_size; j++)
	{
		pop[dist[j]].Crowd_Dist_Set(0.0);
	}
	for (i = 0; i<nobj; i++)
	{
		pop[obj_array[i][0]].Crowd_Dist_Set(INF);
	}
	for (i = 0; i<nobj; i++)
	{
		for (j = 1; j<front_size - 1; j++)
		{
			if (pop[obj_array[i][j]].Crowd_Dist_Show() != INF)
			{
				if (pop[obj_array[i][front_size - 1]].Fitness_Show(i) != pop[obj_array[i][0]].Fitness_Show(i))
				{
					//pop[obj_array[i][j]].Crowd_Dist_Set(pop[obj_array[i][j]].Crowd_Dist_Show() + 0.0);
				
				//else
					double temp;
					if (TGM == false)
						temp = pop[obj_array[i][j]].Crowd_Dist_Show() + (pop[obj_array[i][j + 1]].Fitness_Show(i) - pop[obj_array[i][j - 1]].Fitness_Show(i)) / (pop[obj_array[i][front_size - 1]].Fitness_Show(i) - pop[obj_array[i][0]].Fitness_Show(i));
					else
						temp = pop[obj_array[i][j]].Crowd_Dist_Show() + (pop[obj_array[i][j + 1]].TGM_fitness[0][i] - pop[obj_array[i][j - 1]].TGM_fitness[0][i]) / (pop[obj_array[i][front_size - 1]].TGM_fitness[0][i] - pop[obj_array[i][0]].TGM_fitness[0][i]);
					pop[obj_array[i][j]].Crowd_Dist_Set(temp);
				}
			}
		}
	}
	for (j = 0; j<front_size; j++)
	{
		if (pop[dist[j]].Crowd_Dist_Show() != INF)
		{
			pop[dist[j]].Crowd_Dist_Set((pop[dist[j]].Crowd_Dist_Show()) / nobj);
		}
	}
	return;
}

/* Randomized quick sort routine to sort a population based on a particular objective chosen */
void NSGAII::Front_Obj_Quicksort(std::vector<individual> &pop, int objcount, std::vector<int> &obj_array, int obj_array_size)
{
	Front_Obj_Quicksort_Actual(pop, objcount, obj_array, 0, obj_array_size - 1);
	return;
}

/* Actual implementation of the randomized quick sort used to sort a population based on a particular objective chosen */
void NSGAII::Front_Obj_Quicksort_Actual(std::vector<individual> &pop, int objcount, std::vector<int> &obj_array, int left, int right)
{
	int index;
	int temp;
	int i, j;
	double pivot;
	if (left<right)
	{
		index = Random_I(left, right);
		temp = obj_array[right];
		obj_array[right] = obj_array[index];
		obj_array[index] = temp;
		pivot = pop[obj_array[right]].Fitness_Show(objcount);
		i = left - 1;
		for (j = left; j<right; j++)
		{
			if (pop[obj_array[j]].Fitness_Show(objcount) <= pivot)
			{
				i += 1;
				temp = obj_array[j];
				obj_array[j] = obj_array[i];
				obj_array[i] = temp;
			}
		}
		index = i + 1;
		temp = obj_array[index];
		obj_array[index] = obj_array[right];
		obj_array[right] = temp;
		Front_Obj_Quicksort_Actual(pop, objcount, obj_array, left, index - 1);
		Front_Obj_Quicksort_Actual(pop, objcount, obj_array, index + 1, right);
	}
	return;
}

/* Randomized quick sort routine to sort a population based on crowding distance */
void NSGAII::Dist_Quicksort(std::vector<individual> &pop, std::vector<int> &dist, int front_size)
{
	Dist_Quicksort_Actual(pop, dist, 0, front_size - 1);
	return;
}

/* Actual implementation of the randomized quick sort used to sort a population based on crowding distance */
void NSGAII::Dist_Quicksort_Actual(std::vector<individual> &pop, std::vector<int> &dist, int left, int right)
{
	int index;
	int temp;
	int i, j;
	double pivot;
	if (left<right)
	{
		index = Random_I(left, right);
		temp = dist[right];
		dist[right] = dist[index];
		dist[index] = temp;
		pivot = pop[dist[right]].Crowd_Dist_Show();
		i = left - 1;
		for (j = left; j<right; j++)
		{
			if (pop[dist[j]].Crowd_Dist_Show() <= pivot)
			{
				i += 1;
				temp = dist[j];
				dist[j] = dist[i];
				dist[i] = temp;
			}
		}
		index = i + 1;
		temp = dist[index];
		dist[index] = dist[right];
		dist[right] = temp;
		Dist_Quicksort_Actual(pop, dist, left, index - 1);
		Dist_Quicksort_Actual(pop, dist, index + 1, right);
	}
	return;
}