#include "IBEA.h"
#include <time.h>
#include "Support_Functions.h"

std::vector<std::vector<individual>> IBEA_Archive;	//Archive of individuals for each collective
extern int nfes;
extern time_t elit_t;
extern time_t col_t;
std::vector<double> max_indicator_val;
std::vector<std::vector <std::vector<double>>> indicator_val;

namespace IBEA
{
	void Archive_Create(collective & col);
	void Remove_Dominated(std::vector<individual> & pop);
	void Fitness_Calc(std::vector<individual> & pop, std::vector<short> fit_indexes, short col_ix);
	void IndicatorHD_Calc(std::vector<individual> &pop, std::vector<double> &max_val, std::vector<double> &min_val, std::vector<short> fit_indexes, short col_ix);
	void Remove_Worst(std::vector<individual> &pop, short col_ix);
	void Ind_Fit_Calc(std::vector<individual> &pop, int i_ind, short col_ix);
	double HV_calc(std::vector<double> &ind_1, std::vector<double> &ind_2, short d, std::vector<double> &max_val, std::vector<double> &min_val);
	std::vector<individual> Binary_Tournament(std::vector<individual> & pop, short amt);
}

std::vector<individual> IBEA::IBEA_Calc(collective & col, int iGen)
{
	//get the col index and col size
	int col_size = col.Size_Show();
	int n_obj = col.FCode_Show()[0].Objs();
	short col_ix = col.Index_Show() - 1;

	if (col.Was_Erased())
	{
		IBEA_Pop_Init(n_obj, col);
	}

	time_t elite_t_temp = clock();				//starting time of IBEA

	//create the offspring population
	std::vector<individual> offsprings;

	while (offsprings.size() < col_size)
	{
		std::vector<individual> parents;

		//Get the parents
		parents = Binary_Tournament(IBEA_Archive[col_ix], 2);

		std::vector<individual> child = col.Crossover_NSGAII(parents);

		col.Mutation(child[0]);

		child[0].Fitness_Calc(col.FCode_Show()[0]);
		nfes++;

		offsprings.push_back(child[0]);
	}

	col.Indiv_Set() = offsprings;

	//Create archive
	Archive_Create(col);

	//calculate time
	elit_t += clock() - elite_t_temp;

	//save the inidividuals
	col.save();

	time_t col_t_temp = clock();				//Starting time of collective operations

												//Calculate fitness of the collective
	col.Fitness_Calc();

	//save col time
	col_t += clock() - col_t_temp;


	return IBEA_Archive[col_ix];
}

void IBEA::IBEA_Init(short n_col, population & pop)
{
	IBEA_Archive.clear();
	max_indicator_val.clear();
	indicator_val.clear();
	
	max_indicator_val = std::vector<double>(n_col, 0);
	indicator_val = std::vector<std::vector <std::vector<double>>>(n_col, std::vector <std::vector<double>>());
	IBEA_Archive = std::vector<std::vector<individual>>(n_col, std::vector<individual>());

	//Calculate fitness of the population
	pop.Fitness_Calc();

	nfes += pop.Size_Show();

}

void IBEA::IBEA_Pop_Init(short n_obj, collective & col)
{
	short col_ix = col.Index_Show() - 1;
	int col_size = col.Size_Show();

	IBEA_Archive[col_ix].clear();

	Archive_Create(col);
}

void IBEA::Archive_Create(collective & col)
{
	short col_ix = col.Index_Show() - 1;
	int col_size = col.Size_Show();

	//Create the union
	std::vector<individual> IBEA_Union = col.Indiv_Show();
	IBEA_Union.insert(IBEA_Union.end(), IBEA_Archive[col_ix].begin(), IBEA_Archive[col_ix].end());

	//Remove the dominated points
	Remove_Dominated(IBEA_Union);

	Fitness_Calc(IBEA_Union, col.fit_index[0], col_ix);
	//add the individuals to the archive
	IBEA_Archive[col_ix] = IBEA_Union;

	//refine the archive
	while (IBEA_Archive[col_ix].size() > col_size)
		Remove_Worst(IBEA_Archive[col_ix], col_ix);
}
void IBEA::Remove_Dominated(std::vector<individual> & pop)
{
	short n_obj = pop[0].Fitness_Show().size();
	//Find non-dominated points
	for (int i = 0; i < pop.size(); i++)		//loop for each point
	{
		//copy the fitness
		std::vector<double> ind_i = pop[i].Fitness_Show();
		//Find non-dominated points
		for (int j = i + 1; j < pop.size(); j++)		//loop for each point
		{
			//copy the fitness
			std::vector<double> ind_j = pop[j].Fitness_Show();


			int val = 0;			//temporary value for result saving
			short m_i = 0;			//temporary value for checking if i variable dominates
			short m_j = 0;			//temporary value for checking if j variable dominates

			//check, for every fitness, which point is dominated
			for (short k = 0; k < n_obj; k++)
			{
				//check if i dominates j
				if (ind_i[k] <= ind_j[k])
					m_i++;
				//j dominates i
				else
					m_j++;
			}
			//check if domination or duplication occures (occures if all fitnesses are dominated or all are the same)
			if (m_i == n_obj)
			{
				pop.erase(pop.begin() + j);
				j--;
			}
			else if (m_j == n_obj)
			{
				pop.erase(pop.begin() + i);
				i--;
				break;
			}
		}
	}
}

void IBEA::Fitness_Calc(std::vector<individual> & pop, std::vector<short> fit_indexes, short col_ix)
{
	short n_obj = pop[0].Fitness_Show().size();
	int pop_size = pop.size();

	std::vector<double> max_val(n_obj,-INF);
	std::vector<double> min_val(n_obj, INF);

	for (int i_ind = 0; i_ind < pop_size; i_ind++)
	{
		for (int i_obj = 0; i_obj < n_obj; i_obj++)
		{
			double val = pop[i_ind].Fitness_Show(i_obj);
			if (val > max_val[i_obj])
				max_val[i_obj] = val;
			if (val < min_val[i_obj])
				min_val[i_obj] = val;
		}
	}

	IndicatorHD_Calc(pop, max_val, min_val, fit_indexes, col_ix);

	for (int i_ind = 0; i_ind < pop_size; i_ind++)
	{
		Ind_Fit_Calc(pop, i_ind, col_ix);
	}
}

void IBEA::IndicatorHD_Calc(std::vector<individual> &pop, std::vector<double> &max_val, std::vector<double> &min_val, std::vector<short> fit_indexes, short col_ix)
{
	std::vector<double> ind_i, ind_j;
	short n_obj = pop[0].Fitness_Show().size();
	indicator_val[col_ix].clear();


	for (int i = 0; i < pop.size(); i++)		//loop for each point
	{
		//copy the fitness
		ind_i = pop[i].Fitness_Show();

		std::vector<double> aux;
		//Find non-dominated points
		for (int j = 0; j < pop.size(); j++)		//loop for each point
		{
			//copy the fitness
			ind_j = pop[j].Fitness_Show();

			int flag = Dominance_Check(pop[i], pop[j], fit_indexes);

			double val = 0.;
			if (flag == 1)
			{
				val = -HV_calc(ind_i, ind_j, n_obj, max_val, min_val);
			}
			else
			{
				val = HV_calc(ind_j, ind_i, n_obj, max_val, min_val);
			}

			if (abs(val) > max_indicator_val[col_ix])
				max_indicator_val[col_ix] = abs(val);

			aux.push_back(val);
		}
		indicator_val[col_ix].push_back(aux);

	}
}
void IBEA::Remove_Worst(std::vector<individual> &pop, short col_ix)
{
	double worst = pop[0].utility;
	int ix_worst = 0;
	double kappa = 0.05;

	for (int i = 1; i < pop.size(); i++)
	{
		if (pop[i].utility > worst)
		{
			worst = pop[i].utility;
			ix_worst = i;
		}
	}

	for (int i = 1; i < pop.size(); i++)
	{
		if (i != ix_worst)
		{
			double indicator_fit = pop[i].utility;
			indicator_fit -= exp((-indicator_val[col_ix][ix_worst][i] / max_indicator_val[col_ix]) / kappa);
			pop[i].utility = indicator_fit;
		}
	}

	indicator_val[col_ix].erase(indicator_val[col_ix].begin() + ix_worst);
	for (int ix = 0; ix < indicator_val[col_ix].size(); ix++)
		indicator_val[col_ix][ix].erase(indicator_val[col_ix][ix].begin() + ix_worst);

	pop.erase(pop.begin() + ix_worst);
}
void IBEA::Ind_Fit_Calc(std::vector<individual> &pop, int i_ind, short col_ix)
{
	double indicator_fit = 0;
	double kappa = 0.05;

	for (int i = 0; i < pop.size(); i++)
	{
		if (i != i_ind)
		{
			indicator_fit += exp((-indicator_val[col_ix][i][i_ind] / max_indicator_val[col_ix]) / kappa);
		}
	}
	pop[i_ind].utility = indicator_fit;
}

double IBEA::HV_calc(std::vector<double> &ind_1, std::vector<double> &ind_2, short d, std::vector<double> &max_val, std::vector<double> &min_val)
{
	double a, b, r, max;
	double vol = 0.;
	double rho = 2.0;

	r = rho * (max_val[d - 1] - min_val[d - 1]);
	max = min_val[d - 1] + r;

	a = ind_1[d - 1];
	if (ind_2.empty())
		b = max;
	else
		b = ind_2[d - 1];

	if (d == 1)
	{
		if (a < b)
			vol = (b - a) / r;
		else
			vol = 0;
	}
	else
	{
		if (a < b)
		{
			vol = HV_calc(ind_1, std::vector<double>(), d - 1, max_val, min_val) * (b - a) / r;
			vol += HV_calc(ind_1, ind_2, d - 1, max_val, min_val) * (max - b) / r;
		}
		else
		{
			vol = HV_calc(ind_1, ind_2, d - 1, max_val, min_val) * (max - b) / r;
		}
	}
	return vol;
}
std::vector<individual> IBEA::Binary_Tournament(std::vector<individual> & pop, short amt)
{
	std::vector<individual> parents;
	int pop_size = pop.size();

	for (int ix = 0; ix < amt; ix++)
	{
		int ix_parent1, ix_parent2;
		ix_parent1 = Random_I(0, pop_size - 1);
		ix_parent2 = Random_I(0, pop_size - 1);

		if (pop_size >= 2)
		{
			while (ix_parent2 == ix_parent1)
				ix_parent2 = Random_I(0, pop_size - 1);
		}
		double fit_parent1 = pop[ix_parent1].utility;
		double fit_parent2 = pop[ix_parent2].utility;

		if (fit_parent1 < fit_parent2)
			parents.push_back(pop[ix_parent1]);
		else if (fit_parent1 > fit_parent2)
			parents.push_back(pop[ix_parent2]);
		else
		{
			if (Random() < 0.5)
				parents.push_back(pop[ix_parent1]);
			else
				parents.push_back(pop[ix_parent2]);
		}
	}

	return parents;
}