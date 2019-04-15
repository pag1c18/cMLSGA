
/*This code is based on: L.-Y. Tseng and C. Chen, “Multiple Trajectory Search for Multiobjective Optimization,” in 2007 IEEE Congress on Evolutionary Computation (CEC 2007), 2007, pp. 3609–3616.
It has been modified by Przemyslaw Grudniewski for the purposed of MLSGA-framework, 2019*/

#define _CRT_SECURE_NO_WARNINGS

#include "MTS.h"
#include "Struct.h"
#include "Random.h"
#include <ctime>
#include "Support_Functions.h"



int LSCNT1 = 0;
int LSCNT2 = 0;
int LSCNT3 = 0;
//std::vector<std::vector<double>> RD;			///check!! - probably need to put in individuals
double min_search_range = 10e-8;
std::vector<STRUCTURES::boundaries> bound;
extern int nfes;											//number of function evaluations for MOEAD and MTS
pareto_front * PF_pointer;
std::vector<int> foreground_size;
function * fcode;
bool dynamic_on;
std::vector<std::vector<bool>> improve;
int ncons;
extern int pop_size;	//current population size
extern int dyn_tau;

extern time_t mut_t;
namespace MTS
{
	/*Return the better individual - true if ind1 is better*/
	bool Better1(individual & ind1, individual & ind2, std::vector<short> &fit_indexes);
	/*
	*return the index of the best value given among three*
	*/
	int Better3(int val1, int val2, int val3);
	/*
	*Choose the foreground solutions*
	@param col - given collective
	*/
	void ChooseSolution(collective & col);

	int LocalSearch1(individual & ind, int ix, std::vector<short> &fit_indexes);
	int LocalSearch2(individual & ind, int ix, std::vector<short> &fit_indexes);
	int LocalSearch3(individual & ind, int ix, std::vector<short> &fit_indexes);
	/*
	*Calculate the parameters of an object*
	@param ind - individual for which calculatio nwill be made
	@param save - if results should be saved
	*/
	int Object(individual & ind, bool save);
	/*
	*Fix the value according to boundaries*
	*/
	void Fix_Range(std::vector<double> & code);
	int UpdateSolution(individual & ind);
	/*Create sequence*/
	void CreateSequence(std::vector<int> & seq);

	std::vector<std::vector<unsigned int>> Factor_Create(unsigned int M, unsigned int N, unsigned int Q);
	void Rand_Create(std::vector<std::vector<unsigned int>> & factor, unsigned int M, unsigned int N, unsigned int Q);
	void Create_Array(std::vector<int> & cnt, unsigned int Q);

}


void MTS::MTS_Init(population & pop, pareto_front & PF, short ncol)
{
	PF_pointer = &PF;
	ncons = pop.FCode_Show()[0].Cons();
	//clear ald values
	foreground_size.clear();
	improve.clear();

	//intialise new values
	foreground_size = std::vector<int>(ncol, 0);
	improve = std::vector<std::vector<bool>>(ncol, std::vector<bool>(3, true));


	
	fcode = pop.FCode_Show();
	bound = fcode[0].Bound();
	dynamic_on = fcode[0].Time_Dep();
	int var_n = fcode[0].Vars();
	std::vector<std::vector<unsigned int>> factor;
	int psize = pop.Size_Show();
	//loop for MTS initialisation - do not make any change actually
	if (true)
	{
		factor = Factor_Create(psize, var_n, psize);
		for (int i = 0; i < psize; i++)
		{
			for (int j = 0; j < var_n; j++)
			{
				pop.Indiv_Set(i).Code_Set(j) = bound[j].lower + double(factor[i][j])*(bound[j].upper - bound[j].lower) / double(psize);
			}
		}
	}
	//Set initial values
	for (int i = 0; i < psize; i++)
	{
		//set lcount[i] = 0; zero_count(); - each individual should have one value
		pop.Indiv_Set(i).utility = 0;
		pop.Indiv_Set(i).table = std::vector<int>{ 1 };
		//calculate the fitness of the individuals
		Object(pop.Indiv_Set(i), false);
	}
}

void MTS::MTS_Init_Col(collective & col)
{
	//push the size of the foreground of given collective
	foreground_size[col.Index_Show()-1] = floor(col.Size_Show()*MTS_foreground_multp);
	if (foreground_size[foreground_size.size() - 1] < 1)
		foreground_size[foreground_size.size() - 1] = 1;
	col.save();
}

void MTS::MTS_Calc(collective & col)
{
	//set RD - its probably search range??
	if (col.Was_Erased())
	{
		for (int i = 0; i < col.Size_Show(); i++)
		{
			//set lcount[i] = 0; zero_count(); - each individual should have one value
			col.Indiv_Set(i).utility = 0;
			col.Indiv_Set(i).table = std::vector<int>{ 1 };
		}
		//set the marker to not erased
		col.Clear_Erased();
	}
	//or - check which one is better
	//ChooseSolution(col);

	int ix = col.Index_Show();	//index of the collective

	time_t mut_t_temp = clock();										//starting time of the mutation

	std::vector<short> fit_index = col.fit_index[0];

	//Evolve the population
	for (int i = 0; i < col.Size_Show(); i++)
	{
		if (Termination_Check(dynamic_on))
			break;
		if (col.Indiv_Show(i).table[0] == 1)
		{
			//best solution array?

			//Set local searches to 0;
			LSCNT1 = LSCNT2 = LSCNT3 = 0;

			//Test Local searches
			for (int j = 0; j < MTS_LocalSearch_test_amount; j++)
			{
				LSCNT1 += LocalSearch1(col.Indiv_Set(i), ix, fit_index);
				LSCNT2 += LocalSearch2(col.Indiv_Set(i), ix, fit_index);
				LSCNT3 += LocalSearch3(col.Indiv_Set(i), ix, fit_index);
			}
			//choose the best local search method
			int sw = Better3(LSCNT1, LSCNT2, LSCNT3);

			//perform local searches for given amount of iterations
			for (int j = 0; j < MTS_LocalSearch_amount; j++)
			{
				if (Termination_Check(dynamic_on))
					break;
				switch (sw)
				{
				case 1:
					col.Indiv_Set(i).utility += LocalSearch1(col.Indiv_Set(i), ix, fit_index);
					break;
				case 2:
					col.Indiv_Set(i).utility += LocalSearch1(col.Indiv_Set(i), ix, fit_index);
					break;
				case 3:
					col.Indiv_Set(i).utility += LocalSearch1(col.Indiv_Set(i), ix, fit_index);
					break;
				default:
					abort();
					break;
				}
			}

		}
	}
	//calculate the time of the mutation
	mut_t += clock() - mut_t_temp;

	//Choose the best solutions to the foreground
	ChooseSolution(col);
	
	
	col.save();
}


/*Return the better individual - true if ind1 is better*/
bool MTS::Better1(individual & ind1, individual & ind2, std::vector<short> &fit_indexes)
{
	int a = 0;
	int b = 0;
	


	for (int i = 0; i < fit_indexes.size(); i++)
	{
		short index = fit_indexes[i] - 1;
		if (TGM == false)
		{
			if (ind1.Fitness_Show(index) < ind2.Fitness_Show(index))
				a++;
			else if (ind1.Fitness_Show(index) > ind2.Fitness_Show(index))
				b++;
		}
		else
		{
			if (ind1.TGM_fitness[0][index] < ind2.TGM_fitness[0][index])
				a++;
			else if (ind1.TGM_fitness[0][index]  > ind2.TGM_fitness[0][index])
				b++;
		}
	}
	if (ind1.Cons_Viol_Show() || ind2.Cons_Viol_Show())
	{
		for (int i = 0; i < ncons; i++)
		{
			if (ind1.Cons_Show()[i] < cons_check_param)
				b++;
			if (ind2.Cons_Show()[i] < cons_check_param)
				a++;
		}
	}
	if (a > b)
		return true;
	else
		return false;
}

int MTS::Better3(int val1, int val2, int val3)
{
	if (val1 > val2)
	{
		if (val1 < val3)
			return 3;
		else
			return 1;
	}
	else
	{
		if (val2 < val3)
			return 3;
		else
			return 2;
	}
	
}

void MTS::ChooseSolution(collective & col)
{
	col.Sort_Individuals_MTS();
	for (int i = 0; i < col.Size_Show(); i++)
	{
		if (i < foreground_size[col.Index_Show()-1])
		{
			col.Indiv_Set(i).table[0] = 1;
			col.Indiv_Set(i).utility = 0;
		}
		else
			col.Indiv_Set(i).table[0] = 0;
	}
}

int MTS::LocalSearch1(individual & ind, int ix, std::vector<short> &fit_indexes)
{
	if (Termination_Check(dynamic_on))
		return 0;
	double of;
	bool AllMin;
	int l = 0;	//output
	int n_var = ind.Code_Show().size();			//amount of variables
	std::vector<double> s; // = ( LOWER_BOUND - UPPER_BOUND ) * 0.5;
	std::vector<double> s1; // = s * 0.5;
	std::vector<int> S(n_var,0);
	for (int i = 0; i < n_var; i++)
	{
		s.push_back((bound[i].lower - bound[i].upper) * 0.5);
		s1.push_back(s[i] * 0.5);
	}
	

	if (!improve[ix-1][0])
	{
		AllMin = false;
		for (int i = 0; i < n_var; i++)
		{
			s[i] = s1[i];
			if (std::abs(s[i]) >= min_search_range)
				AllMin = true;
		}
		if (AllMin)
		{
			for (int i = 0; i < n_var; i++)
			{
				s[i] = (bound[i].lower - bound[i].upper) * Random(0.4, 0.5);
				s1[i] = s[i] * 0.5;
			}
		}
		else
		{
			for (int i = 0; i < n_var; i++)
			{
				s1[i] *= 0.5;
			}
		}
	}

	improve[ix - 1][0] = false;

	individual BES = individual(ind);

	CreateSequence(S);
	for (int _i = 0; _i < n_var; _i++)
	{
		//stopping criterion
		if (nfes != 0)
		{
			if (Termination_Check(dynamic_on))
			{
				break;
			}
		}
		int i = S[_i];
		ind.Code_Set(i) += s[i];
		if (Object(ind, true))
			l += MTS_MPL1;
		if (Better1(ind, BES, fit_indexes))
		{
			BES = ind;
			l += MTS_MPL2;
			improve[ix - 1][0] = true;
		}
		else
		{
			//stopping criterion
			if (nfes != 0)
			{
				if (Termination_Check(dynamic_on))
				{
					break;
				}
			}
			if (!Better1(BES, ind, fit_indexes))
			{
				if (Random_I(0,49) == 0)
				{
					BES = ind;
					continue;
				}
			}
			
			ind = BES;
			ind.Code_Set(i) -= s1[i];
			if (Object(ind, true))
				l += MTS_MPL1;
			if (Better1(ind, BES, fit_indexes))
			{
				BES = ind;
				l += MTS_MPL2;
				improve[ix - 1][0] = true;
			}
			else
			{
				ind = BES;
			}
		}
	}
	return l;
}

int MTS::LocalSearch2(individual & ind, int ix, std::vector<short> &fit_indexes)
{
	if (Termination_Check(dynamic_on))
		return 0;
	int i, j, l = 0;
	bool AllMin;
	int n_var = ind.Code_Show().size();			//amount of variables
	std::vector<double> s; // = ( LOWER_BOUND - UPPER_BOUND ) * 0.5;
	std::vector<double> s1; // = s * 0.5;
	for (int i = 0; i < n_var; i++)
	{
		s.push_back((bound[i].lower - bound[i].upper) * 0.5);
		s1.push_back(s[i] * 0.5);
	}
	

	if (!improve[ix - 1][1])
	{
		//ÁY¤p·j´M½d³ò
		AllMin = false;
		for (i = 0; i < n_var; i++)
		{
			s[i] = s1[i];
			if (std::abs(s[i]) > min_search_range)
				AllMin = true;
		}
		if (AllMin)
		{
			for (i = 0; i < n_var; i++)
			{
				s[i] = (bound[i].lower - bound[i].upper) * Random(0.4, 0.5);
				s1[i] = s[i] * 0.5;
			}
		}
		else
		{
			for (i = 0; i < n_var; i++)
			{
				s1[i] *= 0.5;
			}
		}
	}

	bool * ch = new bool[n_var];
	double * D = new double[n_var];
	improve[ix - 1][1] = false;

	individual BES = individual(ind);

	for (i = 0; i < 10; i++)
	{
		//stopping criterion
		if (nfes != 0)
		{
			if (Termination_Check(dynamic_on))
			{
				break;
			}
		}
		//²£¥Í·j´M¤è¦V¥H¤Î­n·j´Mªºdimension
		for (j = 0; j < n_var; j++)
		{
			if (Random_I(0,3) == 0)
				ch[j] = true;
			else
				ch[j] = false;
			if (Random_I(0, 1) == 0)
				D[j] = 1;
			else
				D[j] = -1;
		}
		for (j = 0; j < n_var; j++)
		{
			if (ch[j])
				BES.Code_Set(j) = Random(ind.Code_Show(j) - s[j], ind.Code_Show(j) + s[j]);
			else
				BES.Code_Set(j) = ind.Code_Show(j);
		}
		if (Object(BES, true))
			l += MTS_MPL1;

		if (Better1(BES, ind, fit_indexes))
		{
			ind = BES;
			l += MTS_MPL2;
			improve[ix - 1][1] = true;
		}
		else; // ( _x[n_var] > x[n_var] )
		{
			//stopping criterion
			if (nfes != 0)
			{
				if (Termination_Check(dynamic_on))
				{
					break;
				}
			}
			if (!Better1(ind, BES, fit_indexes))
			{
				if (Random_I(0,9) == 0)
				{
					ind = BES;
					continue;
				}
			}
			for (j = 0; j < n_var; j++)
			{
				if (ch[j])
					BES.Code_Set(j) = Random(ind.Code_Show(j) - s1[j], ind.Code_Show(j) + s1[j]);
				else
					BES.Code_Set(j) = ind.Code_Show(j);
			}
			if (Object(BES, true))
				l += MTS_MPL1;
			if (Better1(BES, ind, fit_indexes))
			{
				ind = BES;
				l += MTS_MPL2;
				improve[ix - 1][1] = true;
			}
		}
	}
	delete ch;
	delete D;
	return l;

}

int MTS::LocalSearch3(individual & ind, int ix, std::vector<short> &fit_indexes)
{
	if (Termination_Check(dynamic_on))
		return 0;
	std::vector<double> U, L, dis;
	int n_var = ind.Code_Show().size();			//amount of variables
	individual BES = individual(ind);
	std::vector<int> S(n_var, 0);
	double M1;
	bool CT_SEARCH = true;
	int grade = 0;


	for (int i = 0; i < n_var; i++)
	{
		U.push_back(bound[i].upper);
		L.push_back(bound[i].lower);
		dis.push_back((bound[i].upper - bound[i].lower) / 10);
	}

	while (CT_SEARCH)
	{
		//stopping criterion
		if (nfes != 0)
		{
			if (Termination_Check(dynamic_on))
			{
				return 0;
			}
		}
		CreateSequence(S);
		for (int _i = 0; _i < n_var; _i++)
		{
			//stopping criterion
			if (nfes != 0)
			{
				if (Termination_Check(dynamic_on))
				{
					return 0;
				}
			}
			int i = S[_i];
			M1 = ind.Code_Show(i);
			while (M1 - dis[i] > L[i])
			{
				M1 = M1 - dis[i];
			}

			M1 = Random( L[i], M1);

			for (ind.Code_Set(i) = M1; ind.Code_Show(i) < U[i]; ind.Code_Set(i) += dis[i])
			{
				
				if (Object(ind, true))
				{
					grade += MTS_MPL1;
				}

				if (Better1(ind, BES, fit_indexes))
				{
					BES = ind;
					grade += MTS_MPL2;
				}
				//stopping criterion
				if (nfes != 0)
				{
					if (Termination_Check(dynamic_on))
					{
						return 0;
					}
				}
			}

			ind = BES;

			if (ind.Code_Set(i) - dis[i] - dis[i] < bound[i].lower)
			{
				L[i] = bound[i].lower;
			}
			else
			{
				L[i] = ind.Code_Show(i) - dis[i] - dis[i];
			}

			if (ind.Code_Set(i) + dis[i] + dis[i] > bound[i].upper)
			{
				U[i] = bound[i].upper;
			}
			else
			{
				U[i] = ind.Code_Show(i) + dis[i] + dis[i];
			}

			dis[i] = (U[i] - L[i]) / 10;
			if (dis[i] < 0.01)
				CT_SEARCH = false;
		}
	}
	return 0;
}


int MTS::Object(individual & ind, bool save)
{
	//Fix_Range(ind.Code_Set());
	if (nfes != 0)
	{
		if (Termination_Check(dynamic_on))
		{
			return 0;
		}
	}
	int isupdate;
	Fix_Range(ind.Code_Set());

	if (TGM == true)
	{
		ind.TGM_fitness.insert(ind.TGM_fitness.begin(), std::vector<double>(fcode[0].Objs(), 0));
		ind.TGM_fitness.pop_back();

		if (ind.TGM_fitness.size() != TGM_size + 1)
			abort();
	}

	ind.Fitness_Calc(fcode[0]);
	nfes++;
	if (save)
		ind.save();

	isupdate = UpdateSolution(ind);

	if (isupdate >= 0)
		return 1;
	else
		return 0;
}
void MTS::Fix_Range(std::vector<double> & code)
{
	for (int i = 0; i < code.size(); i++)
	{
		if (code[i] > bound[i].upper || code[i] < bound[i].lower)
			code[i] = Random(bound[i].lower, bound[i].upper);
	}
}

int MTS::UpdateSolution(individual & ind)
{
	return PF_pointer->Pareto_Search(ind);
}

/*Create sequence*/
void MTS::CreateSequence(std::vector<int> & seq)
{
	int r;
	int seq_size = seq.size();
	for (int i = 0; i < seq_size; i++)
	{
		seq[i] = 0;
	}
	int i = 1;
	while ( i < seq_size)
	{
		r = Random_I(0, seq_size - 1);
		if (seq[r] == 0)
		{
			seq[r] = i;
			i++;
		}
	}
}

std::vector<std::vector<unsigned int>> MTS::Factor_Create( unsigned int M, unsigned int N, unsigned int Q)
{
	std::vector<std::vector<unsigned int>> factor; //output
	for (int i = 0; i < M; i++)
	{
		factor.push_back(std::vector<unsigned int>(N, 0));
	}
	Rand_Create(factor, M, N, Q);

	return factor;
}

void MTS::Rand_Create(std::vector<std::vector<unsigned int>> & factor, unsigned int M, unsigned int N, unsigned int Q)
{
	std::vector<int> cnt = std::vector<int>(Q, 0);
	int k, r;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			if ((j % Q) == 0)
			{
				k = 0;
				Create_Array(cnt, Q);
			}
			do
			{
				r = Random_I(0, M - 1);
			} while (factor[r][i]);
			factor[r][i] = cnt[k++];
		}
	}
}

void MTS::Create_Array(std::vector<int> & cnt, unsigned int Q)
{
	int r;
	for (int i = 0; i < Q; i++)
	{
		cnt[i] = 0;
	}
	for (int i = 0; i < Q; i++)
	{
		do
		{
			r = Random_I(0, Q - 1);
		} while (cnt[r]);
		cnt[r] = i;
	}
}