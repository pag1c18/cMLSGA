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



/*This code is based on Matlab code by Proteek Chandan Roy[1]
Translated to C++ and modified for the purposes of the MLSGA framework by Przemyslaw A.Grudniewski (2019)*/

/* [1] A MATLAB implementation of both NSGA-II and U-NSGA-III.

Copyright [2017] [Proteek Chandan Roy]

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
	 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

This Software runs NSGA-II procedure for different testfunctions
You are free to change, modify, share and distribute this software
Please Acknowledge the author for any use of this software
@author Proteek Chandan Roy, Department of CSE, Michigan State University, USA
email: royprote@egr.msu.edu


 Project Code : YPEA126
 Project Title : Non - dominated Sorting Genetic Algorithm III(NSGA - III)
 Publisher : Yarpiz(www.yarpiz.com)

 Implemented by : S.Mostapha Kalami Heris, PhD(member of Yarpiz Team)

 Contact Info : sm.kalami@gmail.com, info@yarpiz.com

 Base Reference Paper :
 K.Deb and H.Jain, "An Evolutionary Many-Objective Optimization Algorithm 
 Using Reference - Point - Based Nondominated Sorting Approach, Part I : Solving
 Problems With Box Constraints, "
 in IEEE Transactions on Evolutionary Computation,
vol. 18, no. 4, pp. 577 - 601, Aug. 2014.

 Reference Papaer URL : http://doi.org/10.1109/TEVC.2013.2281535
*/




#include "UNSGAIII.h"
#include "NSGAII.h"
#include "Define.h"
#include "Const.h"
#include "Random.h"
#include "Support_Functions.h"
#include <ctime>
#include "LUP.h"

//UNSGAIII have to give output as a vector of individuals or population (vector of individuals should be easier)

//!!!!! CHeck if changing saving to after elitism will give different results (the same with PF checking)
//!!!! write (if was erased) and copy original parent function



extern time_t selec_t;										//Time of Selection
extern time_t mut_t;					//Time of mutation
extern time_t col_t;					//time of collective
extern time_t elit_t;					//Time of UNSGAIII
//extern short MLSt;						//MLSt type

namespace Matrix
{
	double det(std::vector <std::vector<double>> & matrix, int n);
	double norm(std::vector <double> & matrix);
}
namespace UNSGAIII 
{
	void Dirs_Create(int nobj);
	void NSGAIII_Selection(std::vector<individual> & mix_pop, std::vector<individual> & parent_pop, std::vector<short>& fit_indexes);
	std::vector<short> bos(std::vector<individual> &mix_pop, std::vector<short> &feasible_index, std::vector<short>& fit_indexes);
	std::vector<individual> Select(std::vector<individual> & old_pop, std::vector<short> &fit_indexes);
	short Lex_Dominate(std::vector<double> &obj1, std::vector<double> &obj2);
	individual Tournament(individual & ind1, individual & ind2, std::vector<short> & fit_indexes);
	struct CV_Val_Index
	{
		short index;
		float CV_val;
	};
	

	int numdir;  //reference direction - Check generation of them in "Basic_param" - probably for whole population.
	std::vector<std::vector<double>> dirs;
	std::vector<bool> associationsready;
	short col_ix;

}
/*
*Calculate new collective using UNSGAIII algorithm*
@param Col - address of a given collective
*/
void UNSGAIII::UNSGAIII_Calc(collective & col)
{
	//get the number of objectives
	short nobj = col.FCode_Show()[0].Objs();
	col_ix = col.Index_Show() - 1;

	//Copy the parent vector
	std::vector<individual> parent_pop = col.Indiv_Show();		//parent population



	//Check if collective need initial crowding distance calculation and ranking
	if (col.Was_Erased())
	{
		time_t elite_t_temp = clock();				//starting time of whole operation

		//calculate crowding distance and ranking
		if (nobj == 2)
			NSGAII::Rank_Crowding_Distance_Assign(parent_pop, col.fit_index[0]);
		else
			associationsready[col_ix] = false;


		//set flag to has not been erased
		col.Clear_Erased();

		//calculate time
		elit_t += clock() - elite_t_temp;
	}

	std::vector<individual> off_pop;

	//Create the offspring vector
	if (nobj == 2)
		off_pop = col.Crossover_NSGAII(NSGAII::Select(parent_pop, col.fit_index[0])); //offspring population
	else
		off_pop = col.Crossover_NSGAII(Select(parent_pop, col.fit_index[0])); //offspring population

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
	if (nobj == 2) //Do NSGAII
		NSGAII::Nondominated_Sort_Fill(mix_pop, parent_pop, col.fit_index[0]);
	else if (nobj > 2) //Do NSGAIII
		NSGAIII_Selection(mix_pop, parent_pop, col.fit_index[0]);
	else
		abort();

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

void UNSGAIII::Init(short n_obj, short n_col, population & pop)
{
	numdir = pop.Size_Show();
	Dirs_Create(n_obj);
	associationsready = std::vector<bool>(n_col, false);
}
void UNSGAIII::NSGAIII_Selection(std::vector<individual> & mix_pop, std::vector<individual> & parent_pop, std::vector<short>& fit_indexes)
{
	short nobj = fit_indexes.size();
	int mix_pop_size = mix_pop.size();		//n
	int pop_size = parent_pop.size();

	//std::vector<float> opt_R

	std::vector<short> selected_pop_index;
	std::vector<short> feasible_index;
	std::vector<short> infeasible_index;
	for (int i = 0; i < mix_pop.size(); i++)
	{
		if (mix_pop[i].Cons_Viol_Show() == false)
			feasible_index.push_back(i);
		else
			infeasible_index.push_back(i);
	}


	std::vector<std::vector<short>> Fronts(mix_pop_size,std::vector<short>());

	if (!feasible_index.empty())
	{
		//calculate ranks
		std::vector<short> Ranks = bos(mix_pop, feasible_index, fit_indexes);

		for (int ix = 0; ix < feasible_index.size(); ix++)
		{
			Fronts[Ranks[ix]-1].push_back(feasible_index[ix]);
		}

	}

	std::vector<int> pop2Dir(pop_size, 0);
	std::vector<double> pop2DirDistances(pop_size, 0);
	int idx = 0;
	
	associationsready[col_ix] = false;

	if (feasible_index.size() < pop_size)
	{
		selected_pop_index = feasible_index;


		if (!infeasible_index.empty())
		{
			short ncons = mix_pop[0].Cons_Show().size();
			std::vector< CV_Val_Index> CV_index;
			for (int ix = 0; ix < infeasible_index.size(); ix++)
			{
				CV_Val_Index temp_CV;

				temp_CV.index = infeasible_index[ix];

				float temp_cv_val = 0;
				std::vector<double> indi_cons_val = mix_pop[temp_CV.index].Cons_Show();
				for (short ix_cons = 0; ix_cons < ncons; ix_cons++)
				{
					if (indi_cons_val[ix_cons] < 0.)
						temp_cv_val += indi_cons_val[ix_cons];
				}
				temp_CV.CV_val = temp_cv_val;
				CV_index.push_back(temp_CV);
			}

			std::sort(CV_index.begin(), CV_index.end(), [](const CV_Val_Index&a, const CV_Val_Index&b) {return a.CV_val > b.CV_val; });


			for (int ix = 0; ix < infeasible_index.size(); ix++)
			{
				selected_pop_index.push_back(CV_index[ix].index);
				if (selected_pop_index.size() >= pop_size)
					break;
			}
		}
	}
	else
	{
		associationsready[col_ix] = true;

		std::vector<short> cd(mix_pop_size,0);

		cd[0] = Fronts[0].size();
		for (int ix = 1; ix < mix_pop_size; ix++)
		{
			cd[ix] = cd[ix - 1] + Fronts[ix].size();
		}
		short lastfront = 0;
		for (int ix = 0; ix < cd.size(); ix++)
		{
			if (cd[ix] < pop_size)
			{
				continue;
			}
			else
			{
				lastfront = ix;
				break;
			}
		}

		for (int ix = 0; ix < lastfront; ix++)
			selected_pop_index.insert(selected_pop_index.end(), Fronts[ix].begin(), Fronts[ix].end());

		if (selected_pop_index.size() > pop_size)
			abort();

		std::vector<short> lastpopindex = Fronts[lastfront];
		
		std::vector<short> combinepopindex = selected_pop_index;
		combinepopindex.insert(combinepopindex.end(), lastpopindex.begin(), lastpopindex.end());

		short combined_size = combinepopindex.size();

		std::vector<std::vector<double>> combinedobj(nobj, std::vector<double>(combined_size, 0.));

		for (int ix = 0; ix < combined_size; ix++)
		{
			std::vector<double> fit_vect_temp = mix_pop[combinepopindex[ix]].Fitness_Show();
			std::vector<double> fit_vect;
			for (int i_obj = 0; i_obj < nobj; i_obj++)
			{
				fit_vect.push_back(fit_vect_temp[fit_indexes[i_obj] - 1]);
			}
			
			for (int i_obj = 0; i_obj < nobj; i_obj++)
			{
				combinedobj[i_obj][ix] = fit_vect[i_obj];
			}
		}
		std::vector<double> z;
		for (int i_obj = 0; i_obj < nobj; i_obj++)
			z.push_back(*std::min_element(combinedobj[i_obj].begin(), combinedobj[i_obj].end()));

		for (int ix = 0; ix < combined_size; ix++)
		{
			for (int i_obj = 0; i_obj < nobj; i_obj++)
			{
				combinedobj[i_obj][ix] -= z[i_obj];
			}
		}

		std::vector<std::vector<double>> S(nobj, std::vector<double>(nobj, 0));

		for (int i_obj = 0; i_obj < nobj; i_obj++)
		{
			std::vector<double> w(nobj, 1E-16);
			w[i_obj] = 1;

			double norm_w = Matrix::norm(w);

			for (int i = 0; i < nobj; i++)
				w[i] /= norm_w;



			//std::vector<std::vector<double>> tobj = combinedobj;
			int index_min = 0;
			double temp_min = INF;
			for (int i = 0; i < combined_size; i++)
			{
				double temp_max = combinedobj[0][i] / w[0];
				for (int j = 1; j < nobj; j++)
				{
					double temp_val = combinedobj[j][i] / w[j];
					if (temp_val > temp_max)
						temp_max = temp_val;

				}
				if (i == 0)
					temp_min = temp_max;
				else if (temp_max < temp_min)
				{
					temp_min = temp_max;
					index_min = i;
				}
			}
			for (int i = 0; i < nobj; i++)
				S[i_obj][i] = combinedobj[i][index_min];

		}

		std::vector<double> A(nobj, 0);

		if (Matrix::det(S,nobj) < 1E-16)
		{
			for (int i = 0; i < nobj; i++)
				A[i] = *std::max_element(combinedobj[i].begin(), combinedobj[i].end());
		}
		else
		{
			A = LU_Solve(S, std::vector<double>(nobj, 1));
			for (int i = 0; i < nobj; i++)
				A[i] = 1 / A[i];
		}

		//Normalize with intercept
		for (int ix = 0; ix < combined_size; ix++)
		{
			for (int i_obj = 0; i_obj < nobj; i_obj++)
			{
				combinedobj[i_obj][ix] /= A[i_obj];
			}
		}

		std::vector<int> Count(numdir, 0);

		/*for selected indices*/
		std::vector<int>tempPopDir(mix_pop_size, 0);
		std::vector<double>distance(mix_pop_size, INT_MAX);

		for (int i = 0; i < selected_pop_index.size(); i++)
		{
			std::vector<double> obj;
			int s = selected_pop_index[i];


			for (int j = i; j < combined_size; j++)
			{
				if (combinepopindex[j] == s)
				{
					for (int i_obj = 0; i_obj < nobj; i_obj++)
						obj.push_back(combinedobj[i_obj][j]);
					break;
				}
			}


			for (int j = 0; j < numdir; j++)
			{
				std::vector<double> w;
				for (int i_obj2 = 0; i_obj2 < nobj; i_obj2++)
					w.push_back(dirs[j][fit_indexes[i_obj2]-1]);
				double sum_obj = 0;
				double norm_w = Matrix::norm(w);
				std::vector<double> temp_obj(nobj, 0);
				for (int k = 0; k < nobj; k++)
				{
					sum_obj = w[k] * obj[k];
				}
				for (int k = 0; k < nobj; k++)
				{
					temp_obj[k] = obj[k] - ((sum_obj)*w[k]) / (norm_w*norm_w);
				}
				double d = Matrix::norm(temp_obj);


				if (d < distance[s])
				{
					tempPopDir[s] = j;
					distance[s] = d;
				}

			}

			Count[tempPopDir[s]]++;

			pop2Dir[idx] = tempPopDir[s];
			pop2DirDistances[idx] = distance[s];
			idx++;
		}

		//for last indices
		std::vector<std::vector<int>> DirLast(numdir, std::vector<int>());
		std::vector<std::vector<double>> DirLastDist(numdir, std::vector<double>());

		tempPopDir = std::vector<int>(mix_pop_size, 0);
		distance = std::vector<double>(mix_pop_size, INT_MAX);

		for (int i = 0; i < lastpopindex.size(); i++)
		{
			std::vector<double> obj;
			int s = lastpopindex[i];


			for (int j = i; j < combined_size; j++)
			{
				if (combinepopindex[j] == s)
				{
					for (int i_obj = 0; i_obj < nobj; i_obj++)
						obj.push_back(combinedobj[i_obj][j]);
					break;
				}
			}


			for (int j = 0; j < numdir; j++)
			{
				std::vector<double> w;
				for (int i_obj2 = 0; i_obj2 < nobj; i_obj2++)
					w.push_back(dirs[j][fit_indexes[i_obj2]-1]);
				double sum_obj = 0;
				double norm_w = Matrix::norm(w);
				std::vector<double> temp_obj(nobj, 0);
				for (int k = 0; k < nobj; k++)
				{
					sum_obj = w[k] * obj[k];
				}
				for (int k = 0; k < nobj; k++)
				{
					temp_obj[k] = obj[k] - ((sum_obj)*w[k]) / (norm_w*norm_w);
				}
				double d = Matrix::norm(temp_obj);


				if (d < distance[s])
				{
					tempPopDir[s] = j;
					distance[s] = d;
				}

			}

			DirLast[tempPopDir[s]].push_back(s);
			DirLastDist[tempPopDir[s]].push_back(distance[s]);

		}

		//niching preservation
		for (int i = 0; i < numdir; i++)
		{
			if (DirLast[i].size() > 1)
			{
				std::vector<double> DirLastDist_temp = DirLastDist[i];
				std::sort(DirLastDist[i].begin(), DirLastDist[i].end());
				//copy the DirLast
				std::vector<int> DirLast_temp = DirLast[i];
				for (int j = 0; j < DirLast_temp.size(); j++)
				{
					for (int k = 0; k < DirLast_temp.size(); k++)
					{
						if (DirLastDist[i][j] == DirLastDist_temp[k])
						{
							DirLast[i][j] = DirLast_temp[k];
						}
					}
				}
			}
		}

		while (selected_pop_index.size() < pop_size)
		{
			int j = std::distance(Count.begin(), std::min_element(Count.begin(), Count.end()));

			if (DirLast[j].empty())
				Count[j] = INT_MAX;
			else if (Count[j] == 0)
			{
				selected_pop_index.push_back(DirLast[j][0]);
				
			
				pop2Dir[idx] = DirLast[j][0];
				pop2DirDistances[idx] = DirLastDist[j][0];
				idx++;


				Count[j]++;
				DirLast[j].erase(DirLast[j].begin());

			}
			else
			{
				int index = Random_I(0, DirLast[j].size()-1);

				selected_pop_index.push_back(DirLast[j][index]);

				pop2Dir[idx] = DirLast[j][index];
				pop2DirDistances[idx] = DirLastDist[j][index];
				idx++;

				Count[j]++;
				DirLast[j].erase(DirLast[j].begin()+index);
			}
		}

	}
	if (selected_pop_index.size() != pop_size)
		abort();
	parent_pop.clear();
	for (int i = 0; i < pop_size; i++)
	{
		individual temp = mix_pop[selected_pop_index[i]];
		temp.Rank_Set(pop2Dir[i]);
		temp.Crowd_Dist_Set(pop2DirDistances[i]);
		parent_pop.push_back(temp);

	}
}
void UNSGAIII::Dirs_Create(int nobj)
{
	int dirCount = INT_MAX;
	int tempN = numdir;

	while (dirCount > numdir)
	{
		dirs = Uniform_Weights_Generate(tempN, nobj, 3);
		dirCount = dirs.size();
		tempN--;
	}
	numdir = dirs.size();
}

short UNSGAIII::Lex_Dominate(std::vector<double> &obj1, std::vector<double> &obj2)
{
	short equal = 1;
	short d = 1;
	int size = obj1.size();

	for (int i = 0; i < size; i++)
	{
		if (obj1[i] > obj2[i])
		{
			d = 0;
			break;
		}
		else if (equal == 1 && (obj1[i] < obj2[i]))
			equal = 0;
	}
	if (d == 1 && equal == 1)
		d = 0;
	return d;
}

std::vector<short> UNSGAIII::bos(std::vector<individual> &mix_pop, std::vector<short> & feasible_index, std::vector<short>& fit_indexes)
{
	int pop_size = feasible_index.size();
	std::vector<short> rank_vect(pop_size, SHRT_MAX);
	
	int flag;
	int i;
	int end;
	int front_size;
	int rank = 1;
	NSGAII::list *orig;
	NSGAII::list *cur;
	NSGAII::list *temp1, *temp2;
	orig = (NSGAII::list *)malloc(sizeof(NSGAII::list));
	cur = (NSGAII::list *)malloc(sizeof(NSGAII::list));
	front_size = 0;
	orig->index = -1;
	orig->parent = NULL;
	orig->child = NULL;
	cur->index = -1;
	cur->parent = NULL;
	cur->child = NULL;
	temp1 = orig;
	for (i = 0; i < pop_size; i++)
	{
		NSGAII::Insert(temp1, i);
		temp1 = temp1->child;
	}
	do
	{
		if (orig->child->child == NULL)
		{
			rank_vect[orig->child->index] = rank;
			break;
		}
		temp1 = orig->child;
		NSGAII::Insert(cur, temp1->index);
		front_size = 1;
		temp2 = cur->child;
		temp1 = NSGAII::Del(temp1);
		temp1 = temp1->child;
		do
		{
			temp2 = cur->child;
			do
			{
				end = 0;
				flag = Dominance_Check(mix_pop[feasible_index[temp1->index]], mix_pop[feasible_index[temp2->index]], fit_indexes);
				if (flag == 1)
				{
					NSGAII::Insert(orig, temp2->index);
					temp2 = NSGAII::Del(temp2);
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
				NSGAII::Insert(cur, temp1->index);
				front_size++;
				temp1 = NSGAII::Del(temp1);
			}
			temp1 = temp1->child;
		} while (temp1 != NULL);
		temp2 = cur->child;
		do
		{
			rank_vect[temp2->index] = rank;
			temp2 = temp2->child;
		} while (temp2 != NULL);
		temp2 = cur->child;
		do
		{
			temp2 = NSGAII::Del(temp2);
			temp2 = temp2->child;
		} while (cur->child != NULL);
		rank += 1;
	} while (orig->child != NULL);
	free(orig);
	free(cur);

	return rank_vect;

}

/*
*Select individuals for crossover*
@param old_pop - address of the old population
*/
std::vector<individual> UNSGAIII::Select(std::vector<individual> & old_pop, std::vector<short> &fit_indexes)
{
	time_t select_t_temp = clock();										//starting time of the selection
	std::vector<int> a1, a2;				//Value vector
	std::vector<individual> off_pop;			//offspring population - output
	int pop_size = old_pop.size();			//size of the old population


	//Assign the values to the value vector
	for (int i = 0; i < pop_size; i++)
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
individual UNSGAIII::Tournament(individual & ind1, individual & ind2, std::vector<short>& fit_indexes)
{
	int flag;
	//use check dominance specific to selected mode

	if (!ind1.Cons_Viol_Show() && !ind2.Cons_Viol_Show())
	{
		if (associationsready[col_ix] && (ind1.Rank_Show() == ind2.Rank_Show()))
		{
			std::vector<double> obj1 = ind1.Fitness_Show();
			std::vector<double> obj2 = ind2.Fitness_Show();

			flag = Dominance_Check_NCons(ind1, ind2, fit_indexes);

			if (flag == 1)
			{
				return (ind1);
			}
			else if (flag == -1)
			{
				return (ind2);
			}
			else
			{
				if (ind1.Crowd_Dist_Show()< ind2.Crowd_Dist_Show())
					return (ind1);
				else
					return (ind2);
			}

		}
		else
		{
			if ((Random()) <= 0.5)
			{
				return(ind1);
			}
			else
			{
				return(ind2);
			}
		}

		
	}
	else
	{
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
	
}








double Matrix::det(std::vector<std::vector<double>> & matrix, int n) 
{
	double det = 0;
	std::vector <std::vector<double>> submatrix(n - 1, std::vector<double>(n - 1, 0));
	if (n == 2)
		return ((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]));
	else 
	{
		for (int i_lvl = 0; i_lvl < n; i_lvl++) 
		{
			int subi = 0;
			for (int i = 1; i < n; i++) 
			{
				int subj = 0;
				for (int j = 0; j < n; j++) 
				{
					if (j == i_lvl)
						continue;
					submatrix[subi][subj] = matrix[i][j];
					subj++;
				}
				subi++;
			}
			det += (pow(-1, i_lvl) * matrix[0][i_lvl] * Matrix::det(submatrix, n - 1));
		}
	}
	return det;
}

double Matrix::norm(std::vector <double> &matrix)
{
	double norm = 0;

	for (int i = 0; i < matrix.size(); i++)
	{
		norm += pow(matrix[i], 2);
	}

	norm = sqrt(norm);

	return norm;
}