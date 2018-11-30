#include "DMOEADD.h"
#include "Class.h"
#include "Support_Functions.h"
#include "NSGAII.h"
#include <ctime>

extern time_t selec_t;										//Time of Selection
extern time_t mut_t;					//Time of mutation
extern time_t col_t;					//time of collective
extern time_t elit_t;					//Time of DMOEADD
extern time_t cross_t;					//Time of crossover
//extern short MLSt;						//MLSt type
extern int nfes;						//Number of evaluations

namespace DMOEADD {
	short ix_D;								//index of the collective
	int nvar_D;								//number of variables
	int nobj;								//number of objectives

	void Assign_DMOEADD_Param(std::vector<individual> &pop);
	void Assign_DMOEADD_Param2(std::vector<individual> &parent_pop, std::vector<individual> &off_pop);
	std::vector<individual> Select(std::vector<individual> & old_pop, int size);
	individual Tournament(individual & ind1, individual & ind2);
	int Tournament2(individual & ind1, individual & ind2);
	individual Diff_Evo_XoverB(individual &ind0, individual &ind1, individual &ind2, double rate, function & fcode);
	static bool Sort(const individual& c1, const individual& c2);
	static bool Sort2(const individual& c1, const individual& c2);
}



std::vector<individual> DMOEADD::DMOEADD_Calc(collective & col)
{
	//Copy the parent vector
	std::vector<individual> parent_pop = col.Indiv_Show();		//parent population

	std::vector<individual> P_front;		//Pareto front
	
	nobj = col.FCode_Show()[0].Objs();

	nvar_D = col.FCode_Show()[0].Vars();
	
	if (col.Was_Erased())
	{
		//calculate crowding distance and ranking
		Assign_DMOEADD_Param(parent_pop);

		col.Clear_Erased();
	}
	
	
	
	//copy the index of the collective
	ix_D = col.Index_Show();

	

	//Do crossover on specified amt of individuals
	std::vector<individual> off_pop;					//offspring population
	int off_pop_size = col.Size_Show()*DMOEADD_Offspring;	//size of the offspring population
	if (DMOEADD_NEW || !DMEADD_MUTATE_ALL)
	{
		off_pop_size *= 10;

		if (off_pop_size > col.Size_Show())
			off_pop_size = col.Size_Show();
	}
	//Check if offspring population is not too small
	if (off_pop_size < 10)
	{
		if (col.Size_Show() > 20)
			off_pop_size = 10;
		else
			off_pop_size = col.Size_Show() / 2;
	}

	//Select the individuals for crossover
	std::vector<individual> select_pop = Select(parent_pop, off_pop_size); //offspring population

	for (int i = 0; i < off_pop_size; i++)
	{
		double rate2 = 0.5; //rate + 0.25*(rnd_uni(&rnd_uni_init) - 0.5);

		off_pop.push_back(Diff_Evo_XoverB(select_pop[i * 3], select_pop[i * 3 + 1], select_pop[i * 3 + 2], rate2, col.FCode_Show()[0]));

		//Calculate the objectives of the offspring
		if (!DMOEADD_NEW)
		{
			off_pop[i].Fitness_Calc(col.FCode_Show()[0]);
			off_pop[i].utility = 0;
			nfes++;
		}
		
	}
	
	if (DMOEADD_NEW)
	{
		//here we sort all individuals according to fitness and remove the worst ones
		if (DMEADD_MUTATE_ALL)
		{
			//Mutate population
			col.Mutation_NSGAII(off_pop);
			nfes += off_pop.size();
			//Calculate parameters of offsprings
			Assign_DMOEADD_Param2(parent_pop, off_pop);

			//Sort the population
			std::sort(parent_pop.begin(), parent_pop.end(), Sort);

			//Add offsprings
			int parent_size = parent_pop.size();
			for (int i = 0; i < off_pop_size; i++)
			{
				for (int j = 1; j <= off_pop_size * 2; j++)
				{
					int flag = Tournament2(off_pop[i], parent_pop[parent_size - j]);
					int a = 0;
					if (flag == 1)
					{
						parent_pop[parent_size - j] = off_pop[i];
						//Sort the population
						std::sort(parent_pop.begin(), parent_pop.end(), Sort);				//Can change to puting better offsprings into vector (when are better that the best of the worst), and replace at the end
						break;
					}
				}
			}
		}
		else
		{
			abort();
		}
	}
	else
	{
		if (DMEADD_MUTATE_ALL)
		{
			//Calculate parameters of offsprings
			Assign_DMOEADD_Param2(parent_pop, off_pop);

			//Sort the population
			std::sort(parent_pop.begin(), parent_pop.end(), Sort);

			//Add offsprings
			int parent_size = parent_pop.size();
			for (int i = 0; i < off_pop_size; i++)
			{
				for (int j = 1; j <= off_pop_size * 2; j++)
				{
					int flag = Tournament2(off_pop[i], parent_pop[parent_size - j]);
					int a = 0;
					if (flag == 1)
					{
						parent_pop[parent_size - j] = off_pop[i];
						//Sort the population
						std::sort(parent_pop.begin(), parent_pop.end(), Sort);				//Can change to puting better offsprings into vector (when are better that the best of the worst), and replace at the end
						break;
					}
				}
			}

			//Mutate
			std::vector<individual> parent_copy = parent_pop;
			std::vector<individual> parent_copy2 = parent_pop;
			col.Mutation_NSGAII(parent_copy);
			nfes += parent_pop.size();

			//Calculate parameters of offsprings
			Assign_DMOEADD_Param2(parent_copy2, parent_copy);

			//Crate new population (from mutated and offsprings)
			for (int i = 0; i < parent_size; i++)
			{
				if (parent_copy[i].saved_fitness < parent_pop[i].saved_fitness)
					parent_pop[i] = parent_copy[i];
				//parent_pop[i] = Tournament(parent_copy[i], parent_pop[i]);
			}
		}
		else
		{

			//Mutate offsprings
			std::vector<individual> off_copy = off_pop;
			std::vector<individual> off_copy2 = off_pop;
			col.Mutation_NSGAII(off_copy);
			nfes += off_pop.size();
			int off_size =off_pop.size();
			//Calculate parameters of offsprings
			Assign_DMOEADD_Param2(off_copy2, off_copy);

			//Crate new population (from mutated and offsprings)
			for (int i = 0; i < off_size; i++)
			{
				if (off_copy[i].saved_fitness < off_pop[i].saved_fitness)
					off_pop[i] = off_copy[i];
				//off_pop[i] = Tournament(off_copy[i], off_pop[i]);
			}

			nfes += off_pop.size();			

			//Calculate parameters of offsprings
			Assign_DMOEADD_Param2(parent_pop, off_pop);

			//Sort the population
			std::sort(parent_pop.begin(), parent_pop.end(), Sort);

			//Add offsprings
			int parent_size = parent_pop.size();
			for (int i = 0; i < off_pop_size; i++)
			{
				for (int j = 1; j <= off_pop_size * 2; j++)
				{
					int flag = Tournament2(off_pop[i], parent_pop[parent_size - j]);
					if (flag == 1)
					{
						parent_pop[parent_size - j] = off_pop[i];
						//Sort the population
						std::sort(parent_pop.begin(), parent_pop.end(), Sort);				//Can change to puting better offsprings into vector (when are better that the best of the worst), and replace at the end
						break;
					}
				}
			}
		}
	}

	if (!DMOEADD_NEW)
	{
		//Sort the population
		//std::sort(parent_pop.begin(), parent_pop.end(), Sort);
	}
	//Save pareto front
	for (int i = 0; i < parent_pop.size(); i++)
	{
		if (parent_pop[i].Rank_Show() == 0)
			P_front.push_back(parent_pop[i]);
	}

	time_t col_t_temp = clock();				//Starting time of collective operations

	if (parent_pop.size() != col.Indiv_Show().size())
		abort();
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

	return P_front;
}

/*
*Assign DMOEADD parameters to given population*
@param pop - given population
*/
void DMOEADD::Assign_DMOEADD_Param(std::vector<individual> &pop)
{
	time_t elite_t_temp = clock();				//starting time of whole operation

	int psize = pop.size();					//Popualtion size
	std::vector<int> ranks(psize, 0);
	//Calculate and assign pareto-rank values and calculate Z values
	double Z = 0;
	for (int i = 0; i < psize; i++)
	{
		for (int j = i + 1; j < psize; j++)
		{
			int flag = 0;
			if (DMOEADD_OVERRIDE)
			{
				abort();
				//flag = Dominance_Check(pop[i], pop[j]);
				if (flag == -1)
				{
					ranks[i]++;
				}
				else if (flag == 1)
				{
					ranks[j]++;
				}
			}
			else
			{
				abort();
				//flag = Dominance_Check2(pop[i], pop[j],ix_D);
				if (flag == -1)
				{
					ranks[i]++;
				}
				else if (flag == 1)
				{
					ranks[j]++;
				}
			}
		}
		pop[i].Rank_Set(ranks[i]);

		Z += exp(-(double)ranks[i] / (double)DMOEADD_Temp);
	}

	//Sort the individuals according to rank
	std::sort(pop.begin(), pop.end(), Sort2);

	//Calculate S values
	std::vector<double> S;

	for (int i = 0; i < psize; i++)
	{
		double p = (1. / Z)*exp(-(double)pop[i].Rank_Show() / (double)DMOEADD_Temp);
		S.push_back(-p*log(p));
	}

	

	//Calculate and assign crowding distance
	int k = 0;			//how many individuals in current rank
	int j = 0;			//Shows which individual is beginning of rank
	for (int i = 0; i < psize-1; i++)
	{
		
		if (pop[i].Rank_Show() == pop[i + 1].Rank_Show())
		{
			k++;
			continue;
		}
		else
		{
			NSGAII::Crowding_Distance_Indices_Assign(pop, j, j + k);
			k = 0;
			j = i+1;
		}
	}

	//Calculate and assign fitness
	for (int i = 0; i < psize; i++)
	{
		pop[i].saved_fitness.clear();

		//check if crowding distance is not inf
		double d = pop[i].Crowd_Dist_Show();

		if (d > 1.5)
			d = 5;

		double temp_fit = (double)pop[i].Rank_Show() - (S[i]*(double)DMOEADD_Temp)- d;
		pop[i].saved_fitness.push_back(temp_fit);

	}

	//calculate time
	elit_t += clock() - elite_t_temp;
}

/*
*Assign DMOEADD parameters to given population*
@param parent_pop - given parent population
@param off_pop - given offspring population
*/
void DMOEADD::Assign_DMOEADD_Param2(std::vector<individual> &parent_pop, std::vector<individual> &off_pop)
{
	std::vector<individual> parent_pop_copy = parent_pop;

	//Assign utiltiy values
	for (int i = 0; i < off_pop.size(); i++)
	{
		off_pop[i].utility = (double)(i + 1);
	}

	//Add offsprings
	parent_pop_copy.insert(parent_pop_copy.end(), off_pop.begin(), off_pop.end());

	//Calculate values
	Assign_DMOEADD_Param(parent_pop_copy);

	//copy values to offspring and parent
	parent_pop.clear();

	for (int i = 0; i < parent_pop_copy.size(); i++)
	{
		int utility = parent_pop_copy[i].utility;
		if (utility == 0)
			parent_pop.push_back(parent_pop_copy[i]);
		else
			off_pop[utility - 1] = parent_pop_copy[i];
	}

	//reset utility value
	for (int i = 0; i < off_pop.size(); i++)
	{
		off_pop[i].utility = 0;
	}
}


/*
*Select individuals for crossover*
@param old_pop - address of the old population
*/
std::vector<individual> DMOEADD::Select(std::vector<individual> & old_pop, int size)
{
	time_t select_t_temp = clock();										//starting time of the selection
	std::vector<int> a1, a2;				//Value vector
	std::vector<individual> off_pop;			//offspring population - output
	int select_size = size * 3;			//size of the old population
	int p_size = old_pop.size();

											//Assign the values to the value vector
	for (int i = 0; i<select_size; i++)
	{
		a1.push_back(i%p_size);
		a2.push_back(i%p_size);
	}

	//increase the size of the value vector to the multiply of 4
	int mult_size = 0;			//How many number have been added
	while ((a1.size() % 4) != 0)
	{
		a1.push_back(Random_I(0, p_size - 1));
		a2.push_back(Random_I(0, p_size - 1));
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
	for (int i = 0; i < select_size; i += 4)
	{
		off_pop.push_back(Tournament(old_pop[a1[i]], old_pop[a1[i + 1]]));
		off_pop.push_back(Tournament(old_pop[a1[i + 2]], old_pop[a1[i + 3]]));
		off_pop.push_back(Tournament(old_pop[a2[i]], old_pop[a2[i + 1]]));
		off_pop.push_back(Tournament(old_pop[a2[i + 2]], old_pop[a2[i + 3]]));
	}

	//remove the additional individuals
	for (int i = 0; i < mult_size; i++)
		off_pop.pop_back();

	//calculate the time of the selection
	selec_t += clock() - select_t_temp;

	if (off_pop.size() != size*3)
		abort();
	return off_pop;
}

/*
*Do binary Tournament*
@param ind1 - address of 1st (new) individual
@param ind2 - address of 2nd (old) individual
*/
individual DMOEADD::Tournament(individual & ind1, individual & ind2)
{
	if (ind1.Cons_Viol_Show() && !ind2.Cons_Viol_Show())
	{
		return ind2;
	}
	else if (!ind1.Cons_Viol_Show() && ind2.Cons_Viol_Show())
	{
		return ind1;
	}

	if (ind1.Rank_Show() < ind2.Rank_Show())
	{
		return (ind1);
	}
	else if (ind1.Rank_Show() == ind2.Rank_Show() && ind1.Crowd_Dist_Show() > ind2.Crowd_Dist_Show())
	{
		return (ind1);
	}
	else if (exp((double)(ind1.Rank_Show() - ind2.Rank_Show())/(double)DMOEADD_Temp)>Random())
	{
		return (ind1);
	}
	else
	{
		return(ind2);
	}
}

/*
*Do binary Tournament*
@param ind1 - address of 1st (new) individual
@param ind2 - address of 2nd (old) individual
*/
int DMOEADD::Tournament2(individual & ind1, individual & ind2)
{
	if (ind1.Cons_Viol_Show() && !ind2.Cons_Viol_Show())
	{
		return 2;
	}
	else if (!ind1.Cons_Viol_Show() && ind2.Cons_Viol_Show())
	{
		return 1;
	}

	if (ind1.Rank_Show() < ind2.Rank_Show())
	{
		return 1;
	}
	else if (ind1.Rank_Show() == ind2.Rank_Show() && ind1.Crowd_Dist_Show() > ind2.Crowd_Dist_Show())
	{
		return 1;
	}
	else if (exp((double)(ind1.Rank_Show() - ind2.Rank_Show()) / (double)DMOEADD_Temp)>Random())
	{
		return 1;
	}
	else
	{
		return 2;
	}
}

individual DMOEADD::Diff_Evo_XoverB(individual &ind0, individual &ind1, individual &ind2, double rate, function & fcode)
{
	time_t cross_t_temp = clock();
	int idx_rnd = int(Random()*nvar_D);

	std::vector<double> child_code = std::vector<double>(nvar_D, 0.0);

	for (int n = 0; n<nvar_D; n++)
	{
		double upper_b = fcode.Bound(n, "upper");	//upper boundary
		double lower_b = fcode.Bound(n, "lower");	//lower boundary
													/*Selected Two Parents*/

													// strategy one 
													// child_code.x_var[n] = ind0.x_var[n] + rate*(ind2.x_var[n] - ind1.x_var[n]);

													//*
													// strategy two
		double rnd1 = Random();
		double CR = 1.0;
		if (rnd1<CR || n == idx_rnd)
			child_code[n] = ind0.Code_Show(n) + rate*(ind2.Code_Show(n) - ind1.Code_Show(n));
		else
			child_code[n] = ind0.Code_Show(n);
		//*/


		// handle the boundary voilation
		if (child_code[n]<lower_b)
		{
			double rnd = Random();
			child_code[n] = lower_b + rnd*(ind0.Code_Show(n) - lower_b);
		}
		if (child_code[n]>upper_b)
		{
			double rnd = Random();
			child_code[n] = upper_b - rnd*(upper_b - ind0.Code_Show(n));
		}
	}
	cross_t += clock() - cross_t_temp;
	return individual{ child_code,fcode };
}

//sort according to fitness
static bool DMOEADD::Sort(const individual& c1, const individual& c2)
{
	return c1.saved_fitness[0] < c2.saved_fitness[0];
}

//sort according to rank
static bool DMOEADD::Sort2(const individual& c1, const individual& c2)
{
	return c1.Rank_Show() < c2.Rank_Show();
}