/*==========================================================================
//  C++ Implementation of MOEA/D Based on Differential Evolution (DE) for Contest Multiobjective
//  Problems in CEC2009
//
//  Author: Hui Li
//
//  See the details of MOEA/D-DE and test problems in the following papers
//
//  1) H. Li and Q. Zhang, Comparison Between NSGA-II and MOEA/D on a Set of Multiobjective Optimization
//  Problems with Complicated Pareto Sets, Technical Report CES-476, Department of Computer Science,
//  University of Essex, 2007
//
//  2) H. Li and Q. Zhang, Multiobjective Optimization Problems with Complicated Pareto Sets, MOEA/D and NSGA-II,
//  IEEE Transaction on Evolutionary Computation, 2008, to appear.
//
//  If you have any questions about the codes, please contact
//  Dr. Hui Li       at   hzl@cs.nott.ac.uk   or
//  Dr. Qingfu Zhang at   qzhang@essex.ac.uk
//
//  Date: 14/NOV/2008
//
// ===========================================================================*/

/*Modified by Przemyslaw Grudniewski for the purposed of MLSGA-framework, 2019*/


#define _CRT_SECURE_NO_WARNINGS

#include "MOEAD.h"
#include "Define.h"
#include "Const.h"
#include "Clustering.h"
#include "Sobol.h"
#include "Support_Functions.h"
#include <ctime>

int nfes;						//number of function evaluations
int nobj;						//number of objectives
int nvar;						//number of variables
int nweightindex;				//count of wieght indexes used
//std::vector<individual> saved;	//vector of saved individuals

extern std::string MODE;											//current hybrid mode
extern time_t selec_t;										//Time of Selection
extern time_t mut_t;					//Time of mutation
extern time_t col_t;					//time of collective
extern time_t elit_t;					//Time of 
extern time_t cross_t;					//Time of crossover
extern int pop_size;					//population size

std::vector<double> idealpoint;			//Storage of ideal point
std::vector<double> nadirpoint;
std::vector<individual> PF_storage2;	//Storage of possible PF points
std::vector<std::vector<double>> weight_vector;		//Storage of weight vectors

namespace MOEAD
{
	double Dist_Vector(std::vector <double> &vec1, std::vector <double> &vec2);
	void Min_Fast_Sort(std::vector<double> &x, std::vector<int> &idx, int n, int m);
	void Evolve(collective & col, int iGen);
	void Tournament(int depth, std::vector<int> &selected, int csize, collective & col);
	void Mate_Selection(std::vector<int> &list, int cid, int size, int type, collective & col);

	void Permutation(std::vector<int> &perm);

}

/*
Calculate individuals using NSGAII algorithm
@param Col - address of a given collective
@param iGen - current generation
*/
std::vector<individual> MOEAD::MOEAD_Calc(collective & col, int iGen)
{
	
	//check if collective was erased or is new
	if (col.Was_Erased())
	{
		time_t elite_t_temp = clock();				//starting time of whole operation
		//create the neighbourhood
		MOEAD::Neighbourhood_Init(col);

		//set the marker to not erased
		col.Clear_Erased();

		//calculate time
		elit_t += clock() - elite_t_temp;
	}
	PF_storage2.clear();
	Evolve(col, iGen);

	
	//save the inidividuals
	col.save();

	time_t col_t_temp = clock();				//Starting time of collective operations

	//Calculate fitness of the collective
	col.Fitness_Calc();

	//save col time
	col_t += clock() - col_t_temp;
	return PF_storage2;
}

void MOEAD::MOEAD_Init(short n_obj, population & pop)
{
	time_t elite_t_temp = clock();				//starting time of whole operation

	nobj = n_obj;
	nvar = pop.Indiv_Show(1).Code_Show().size();

	//Calculate fitness of the population
	pop.Fitness_Calc();

	int psize = pop.Size_Show();
	nfes += psize;


	//clear the old weight_vector
	weight_vector.clear();
	if (MODE != "MOEADM2M")
	{
		//get namda for whole population
		if (MOEAD_weight_assign == "Pop" || MOEAD_weight_assign == "Col")
		{
			//Read the values from the file
			if (MOEAD_weight_generator == "File")
			{
				// Read weight vectors from a data file
				char filename[1024];
				sprintf(filename, "Input/MOEADWeight/W%dD_%d.dat", nobj, pop.Size_Show());
				std::ifstream readf(filename);


				for (int i = 0; i <  psize; i++)
				{
					std::vector<double> temp_namda = std::vector<double>(nobj, 0.);
					// Load weight vectors
					for (int j = 0; j < nobj; j++)
					{
						readf >> temp_namda[j];
						//printf("%f ", sub.namda[j]);
					}
					weight_vector.push_back(temp_namda);
					pop.Indiv_Set(i).namda = temp_namda;
				}
			}
			//Generate by uniform method
			else if (MOEAD_weight_generator == "Uniform")
			{
				weight_vector = Uniform_Weights_Generate(psize, nobj, 1);
			}
			//Generate by sobol
			else if (MOEAD_weight_generator == "Sobol")
			{
				abort();	// sobol not properly implemented
				//weight_vector = Sobol_Sequence(psize, nobj);
			}
			//wrong mode chosen
			else
				abort();

			while (weight_vector.size() < pop.Size_Show())
			{
				pop.Remove();
				psize--;
			}



			if (MOEAD_weight_assign == "Pop")
			{
				//Assign the weight vectors to whole population
				for (int i = 0; i <  psize; i++)
					pop.Indiv_Set(i).namda = weight_vector[i];
			}
			else // Collective assignment
				nweightindex = 0;
		}
		else if (MOEAD_weight_assign == "Sep")
		{
			//empty check - but abort if the file generation is chosen (cannot read from file for separate vectors)
			if (MOEAD_weight_generator == "File")
				abort();
		}
		else // Wrong mode chosen
			abort();
	}
	


	//calculate time
	elit_t += clock() - elite_t_temp;
}


void MOEAD::Population_Init(short n_obj, collective & col)
{
	time_t elite_t_temp = clock();				//starting time of whole operation

	//Get the size of the collective
	int csize = col.Size_Show();

	if (MODE != "MOEADM2M")
	{
		//Get namda for collective when calculated by collective
		if (MOEAD_weight_assign == "Col")
		{
			//get weights
			for (int i = 0; i < csize; i++)
			{
				col.Indiv_Set(i).namda = weight_vector[nweightindex];

				nweightindex++;

				//DEBUG
				if (nweightindex >= weight_vector.size())
					abort();
			}
		}
		//If population type the values are already assigned but do check anyway
		else if (MOEAD_weight_assign == "Pop")
			; //
		//else if each collective have separate weight vectors
		else if (MOEAD_weight_assign == "Sep")
		{

			if (MOEAD_weight_generator == "Uniform")
			{
				weight_vector = Uniform_Weights_Generate(csize, nobj , 1);
			}
			//Generate by sobol
			else if (MOEAD_weight_generator == "Sobol")
			{
				abort();	// sobol not properly implemented
							//weight_vector = Sobol_Sequence(csize, nobj);
			}
			//wrong mode chosen
			else
				abort();

			//get weights
			for (int i = 0; i < csize; i++)
			{
				col.Indiv_Set(i).namda = weight_vector[i];
			}
		}
		//wrong mode chosen
		else
			abort();
	}
	


	for (int i = 0; i < csize; i++)
	{
		//save the fitness of the previous individual
		if (TGM == false)
			col.Indiv_Set(i).saved_fitness = col.Indiv_Show(i).Fitness_Show();
		else
			col.Indiv_Set(i).saved_fitness = col.Indiv_Show(i).TGM_fitness[0];

		col.Indiv_Show(i).save();
	}
	
	//calculate time
	elit_t += clock() - elite_t_temp;
}

void MOEAD::Neighbourhood_Init(collective & col)
{
	//Get the size of the collective
	int csize = col.Size_Show();

	std::vector<double> dist = std::vector<double>(csize, 0);
	std::vector<int>    indx = std::vector<int>(csize, 0);

	//Calculate the sie of the neighbourhood
	int MOEAD_niche = csize * MOEAD_niche_multi;				//the neighbour size

	//Check if is not too small
	if (MOEAD_niche <= 3)
		MOEAD_niche = 4;

	for (int i = 0; i<csize; i++)
	{
		if (col.Indiv_Show(i).saved_fitness.size() != nobj)
		{
			//save the fitness of the previous individual
			if (TGM == false)
				col.Indiv_Set(i).saved_fitness = col.Indiv_Show(i).Fitness_Show();
			else
				col.Indiv_Set(i).saved_fitness = col.Indiv_Show(i).TGM_fitness[0];
		}

		//clear the old neighbourhood
		col.Indiv_Set(i).table.clear();
		// calculate the distances based on weight vectors
		for (int j = 0; j<csize; j++)
		{
			dist[j] = Dist_Vector(col.Indiv_Show(i).namda, col.Indiv_Show(j).namda);
			indx[j] = j;
		}

		// find 'niche' nearest neighboring subproblems
		Min_Fast_Sort(dist, indx, col.Size_Show(), MOEAD_niche);



		// save the indexes of the nearest 'niche' neighboring weight vectors
		for (int k = 0; k<MOEAD_niche; k++)
		{
			col.Indiv_Set(i).table.push_back(indx[k]);
		}
	}
	dist.clear();
	indx.clear();
}

void MOEAD::Utility_Comp(collective & col, int iGen)
{

	//Get the size of the collective
	int csize = col.Size_Show();

	double f1, f2, uti, delta;
	std::vector<short> fit_ix = col.fit_index[0];

	//calculate the utility
	for (int n = 0; n<csize; n++)
	{
		if (TGM == false)
			f1 = Fitness_Function(col.Indiv_Show(n).Fitness_Show(), col.Indiv_Show(n).namda, fit_ix, iGen);
		else
			f1 = Fitness_Function(col.Indiv_Show(n).TGM_fitness[0], col.Indiv_Show(n).namda, fit_ix, iGen);
		f2 = Fitness_Function(col.Indiv_Show(n).saved_fitness, col.Indiv_Show(n).namda, fit_ix, iGen);
		delta = (f2 - f1) / f2;
		//delta = (f2 - f1);
		if (delta>0.001)
		{
			col.Indiv_Set(n).utility = 1.0;
		}
		else
		{
			//uti        = 0.95*(1.0+delta/0.001)*utility[n];
			//utility[n] = uti<1.0?uti:1.0;
			col.Indiv_Set(n).utility *= (0.95 + 0.05*delta / 0.001);
		}
		if (TGM == false)
			col.Indiv_Set(n).saved_fitness = col.Indiv_Show(n).Fitness_Show();
		else
			col.Indiv_Set(n).saved_fitness = col.Indiv_Show(n).TGM_fitness[0];
	}
}

double MOEAD::Dist_Vector(std::vector <double> &vec1, std::vector <double> &vec2)
{
	int dim = vec1.size();
	double sum = 0;
	for (int n = 0; n<dim; n++)
		sum += (vec1[n] - vec2[n])*(vec1[n] - vec2[n]);
	return sqrt(sum);
}

void MOEAD::Min_Fast_Sort(std::vector<double> &x, std::vector<int> &idx, int n, int m)
{
	for (int i = 0; i<m; i++)
	{
		for (int j = i + 1; j<n; j++)
			if (x[i]>x[j])
			{
				double temp = x[i];
				x[i] = x[j];
				x[j] = temp;
				int id = idx[i];
				idx[i] = idx[j];
				idx[j] = id;
			}
	}
}

/*
*Evolve the current collective*
@param col - current collective
*/
void MOEAD::Evolve(collective & col, int iGen)
{
	std::vector<int> order;		//order of subproblems

	//Do tournament
	Tournament(10, order, col.Size_Show(), col);

	for (int sub = 0; sub < order.size(); sub++)
	{

		int c_sub = order[sub];    // random order

		int type;
		double rnd = Random();

		// mating selection based on probability
		if (rnd < MOEAD_mating_chance)
			type = 1;   // from neighborhood
		else            
			type = 2;   // from population

		// select the indexes of mating parents
		std::vector<int> plist;

		Mate_Selection(plist, c_sub, 2, type, col);  // neighborhood selection


		
		
		double rate2 = 0.5; //rate + 0.25*(rnd_uni(&rnd_uni_init) - 0.5);
		

		//create dummy child
		individual child = col.Indiv_Show(0);

		// produce a child solution
		if (MOEAD_evol_operation_type == 1)
			child = Diff_Evo_XoverB(col.Indiv_Show(c_sub), col.Indiv_Show(plist[0]), col.Indiv_Show(plist[1]), rate2, col.FCode_Show()[0], col.CCode_Show()[0], col.GAPara_Show()[0]);
		else if (MOEAD_evol_operation_type == 2)
			child = LL_Crossover(col.Indiv_Show(plist[0]), col.Indiv_Show(plist[1]), iGen, col.FCode_Show()[0]);
		else
			abort();
		plist.clear();

		time_t mut_t_temp = clock();
		// apply mutation
		if (MOEAD_evol_operation_type == 1)
			col.Mutation(child);
		else if (MOEAD_evol_operation_type == 2)
			LL_Mutation(child, iGen, col.GAPara_Show()[0], col.FCode_Show()[0]);
		else
			abort();
		

		mut_t += clock() - mut_t_temp;

		time_t elite_t_temp = clock();				//starting time of MOEAD

		// evaluate the child solution
		col.population::Fitness_Calc(child);
		//child.save();

		if (MODE == "MOEADMSF" || MODE == "MOEADPSF")
			Problem_Update_Global(child, col, c_sub, type, iGen);
		else
			Problem_Update(child, col, c_sub, type, iGen);

		nfes++;

		//calculate time
		elit_t += clock() - elite_t_temp;

		if (Termination_Check(col.FCode_Show()[0].Time_Dep()))
			break;
	}
}

/*
*Tournament selection*
@param depth - depth of selection
@param selected - vector of indexes of selected individuals
@param csize - size of the collective/population
@param cix - index of the collective
*/
void MOEAD::Tournament(int depth, std::vector<int> &selected, int csize, collective & col)
{
	time_t select_t_temp = clock();										//starting time of the selection
	// selection based on utility
	std::vector<int> candidate;
	for (int k = 0; k < nobj; k++)
		selected.push_back(k);   // select first m weights
	
	for (int n = nobj; n < csize; n++)    
		candidate.push_back(n);  // set of unselected weights

	while (selected.size() < int((float)csize / MOEAD_new_sol_multi))
	//while (selected.size() < int(csize))
	{
		int best_idd = int(Random()*candidate.size()), i2;
		int best_sub = candidate[best_idd], s2;
		for (int i = 1; i<depth; i++)
		{
			i2 = int(Random()*candidate.size());
			s2 = candidate[i2];
			if (col.Indiv_Show(s2).utility>col.Indiv_Show(best_sub).utility)
			{
				best_idd = i2;
				best_sub = s2;
			}
		}
		selected.push_back(best_sub);
		candidate.erase(candidate.begin() + best_idd);
	}
	//calculate the time of the selection
	selec_t += clock() - select_t_temp;
}

/*
*Parents selection*
@param list - the set of the indexes of selected mating parents
@param cid  - the id of current subproblem
@param size - the number of selected mating parents
@param type - 1 - neighborhood; otherwise - whole population
*/
void MOEAD::Mate_Selection(std::vector<int> &list, int cid, int size, int type, collective & col)
{
	time_t select_t_temp = clock();										//starting time of the selection
	int ss = col.Indiv_Show(cid).table.size(), id, parent;
	while (list.size()<size)
	{
		if (type == 1) 
		{
			id = int(ss*Random());
			parent = col.Indiv_Show(cid).table[id];
		}
		else
			parent = int(col.Size_Show()*Random());

		// avoid the repeated selection
		bool flag = true;
		for (int i = 0; i<list.size(); i++)
		{
			if (list[i] == parent) // parent is in the list
			{
				flag = false;
				break;
			}
		}

		if (flag) list.push_back(parent);
	}
	//calculate the time of the selection
	selec_t += clock() - select_t_temp;
}


individual MOEAD::Diff_Evo_XoverB(individual &ind0, individual &ind1, individual &ind2, double rate, function & fcode, crossover<individual> & ccode, GA_parameters & gapara)
{
	if (ENCODING == "Real")
	{
		nvar = ind0.Code_Show().size();
		time_t cross_t_temp = clock();
		int idx_rnd = int(Random()*nvar);



		std::vector<double> child_code = std::vector<double>(nvar, 0.0);

		for (int n = 0; n < nvar; n++)
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
			if (rnd1 < CR || n == idx_rnd)
				child_code[n] = ind0.Code_Show(n) + rate*(ind2.Code_Show(n) - ind1.Code_Show(n));
			else
				child_code[n] = ind0.Code_Show(n);
			//*/


			// handle the boundary violation
			if (child_code[n] < lower_b)
			{
				double rnd = Random();
				child_code[n] = lower_b + rnd*(ind0.Code_Show(n) - lower_b);
			}
			if (child_code[n] > upper_b)
			{
				double rnd = Random();
				child_code[n] = upper_b - rnd*(upper_b - ind0.Code_Show(n));
			}
		}
		cross_t += clock() - cross_t_temp;
		individual out_ind(child_code, fcode);			//output individual
		if (TGM == true)
		{
			out_ind.TGM_fitness = ind0.TGM_fitness;
			out_ind.TGM_fitness.insert(out_ind.TGM_fitness.begin(), std::vector<double>(fcode.Objs(), 0));
			out_ind.TGM_fitness.pop_back();

			if (out_ind.TGM_fitness.size() != TGM_size + 1)
				abort();
		}
		return out_ind;
	}
	else if (ENCODING == "Binary" || ENCODING == "Gray")
	{
		//Create the vector of selected individuals for the crossover
		std::vector<individual> selected;
		selected.push_back(ind0);
		if (Random() < 0.5)
			selected.push_back(ind1);
		else
			selected.push_back(ind2);

		std::vector<individual> new_pop;		//output from crossover
		new_pop = ccode.Crossover(selected, gapara, fcode);

		if (Random() < 0.5)
			return new_pop[0];
		else
			return new_pop[1];
	}
	else
		abort();
}

/*
*Update the problem*
@param indiv - child solution
@param col - given collective/population
@param id - the id of current subproblem
@param type - update solutions in - neighborhood (1) or whole population (otherwise)
*/
void MOEAD::Problem_Update(individual &indiv, collective & col, int &id, int &type, int iGen)
{
	
	int size, time = 0;
	int save_temp = 0;
	if (type == 1)	
		size = col.Indiv_Show(id).table.size();  // from neighborhood
	else        
		size = col.Size_Show();            // from whole population

													 // a random order to update
	std::vector<int> perm;
	for (int k = 0; k < size; k++)
		perm.push_back(k);

	//random_shuffle(perm.begin(), perm.end());
	// replaced by Aimin 2011.04.29
	Permutation(perm);

	//Copy NPCp for ease of operation
	std::vector<individual> temp_col = col.Indiv_Show();


	int MOEAD_limit = col.Size_Show() * MOEAD_limit_multi;				//maximal number of solutions replaced
	if (MOEAD_limit == 1)
		MOEAD_limit = 2;

	std::vector<short> fit_ix = col.fit_index[0];


	for (int i = 0; i<size; i++)
	{
		// Pick a subproblem to update
		int k;
		if (type == 1) 
			k = temp_col[id].table[perm[i]];
		else       
			k = perm[i];

		// calculate the values of objective function regarding the current subproblem
		double f1, f2;
		if (TGM == false)
		{
			f1 = Fitness_Function(temp_col[k].Fitness_Show(), temp_col[k].namda, fit_ix, iGen);
			f2 = Fitness_Function(indiv.Fitness_Show(), temp_col[k].namda, fit_ix, iGen);
		}
		else
		{
			f1 = Fitness_Function(temp_col[k].TGM_fitness[0], temp_col[k].namda, fit_ix, iGen);
			f2 = Fitness_Function(indiv.TGM_fitness[0], temp_col[k].namda, fit_ix, iGen);
		}
		if (f2<f1)
		{
			indiv.namda = temp_col[k].namda;
			indiv.table = temp_col[k].table;
			indiv.saved_fitness = temp_col[k].saved_fitness;
			temp_col[k] = indiv;
			col.Indiv_Set(k) = indiv;
			if (save_temp == 0 && (MODE == "MOEAD" || MODE == "MOEADMSF" || MODE == "MOEADPSF" || MODE == "MOEADM2M"))
			{
				PF_storage2.push_back(indiv);
				save_temp++;
			}
			time++;
		}
		// the maximal number of solutions updated is not allowed to exceed 'limit'
		if (time >= MOEAD_limit && type != 3)
		{
			perm.clear();
			return;
		}
	}
	perm.clear();
}

void MOEAD::Problem_Update_Global(individual &indiv, collective & col, int &id, int &type, int iGen)
{
	int pops = col.Size_Show();
	int save_temp = 0;
	std::vector<double> x(pops, 0.);
	std::vector<int> idx(pops, 0);

	int MOEAD_limit = col.Size_Show() * MOEAD_limit_multi;				//maximal number of solutions replaced

	int time = 0;

	int niche = col.Indiv_Show(0).table.size();
	//Copy col for ease of operation
	std::vector<individual> temp_col = col.Indiv_Show();
	
	std::vector<short> fit_ix = col.fit_index[0];

	
	for (int i = 0; i < pops; i++)
	{
		idx[i] = i;
		
		x[i] = Fitness_Function(indiv.Fitness_Show(), temp_col[i].namda, fit_ix, iGen);

	}
	Min_Fast_Sort(x, idx, pops, niche);
	//update_problem(indiv, idx[0]);

	int j;
	for (j = 0; j < niche; j++)
	{

		double fj;
		
		fj = Fitness_Function(temp_col[idx[j]].Fitness_Show(), temp_col[idx[j]].namda, fit_ix, iGen);
		
		if (fj>x[j])
		{
			indiv.namda = temp_col[idx[j]].namda;
			indiv.table = temp_col[idx[j]].table;
			indiv.saved_fitness = temp_col[idx[j]].saved_fitness;
			temp_col[idx[j]] = indiv;
			col.Indiv_Set(idx[j]) = indiv;
			if (save_temp == 0 && (MODE == "MOEAD" || MODE == "MOEADMSF" || MODE == "MOEADPSF" || MODE == "MOEADM2M"))
			{
				PF_storage2.push_back(indiv);
				save_temp++;
			}
			time++;
		}
		if (time >= MOEAD_limit && type != 3)
		{
			break;
		}
	}
}

individual MOEAD::LL_Crossover(individual &ind1, individual &ind2, int iGen, function & fcode)
{
	if (ENCODING != "Real")
		abort();

	double gen_val;

	if (T_con == "nfes")
	{
		gen_val = nfes / max_iterations;

	}
	else
	{
		gen_val = iGen / max_generations;
	}

	time_t cross_t_temp = clock();
	
	//temporary code of the child
	std::vector<double> child_code = std::vector<double>(nvar, 0.);


	//copy the boundary vector
	std::vector<STRUCTURES::boundaries> var_bound = fcode.Bound();

	//copy the parents codes
	std::vector<double> ind1_code = ind1.Code_Show();
	std::vector<double> ind2_code = ind2.Code_Show();

	int i;
	double a = pow(1.0 - 1.0*gen_val, 0.7);
	double rc = 2.0*(Random() - 0.5)*(1 - pow(Random(), -a));
	for (i = 0; i<nvar; i++) 
	{
		child_code[i] = ind1_code[i] + rc*(ind1_code[i] - ind2_code[i]);
		if (child_code[i]<var_bound[i].lower)
			child_code[i] = var_bound[i].lower + 0.5*Random()*(ind1_code[i] - var_bound[i].lower);
		if (child_code[i]>var_bound[i].upper)
			child_code[i] = var_bound[i].upper - 0.5*Random()*(var_bound[i].upper - ind1_code[i]);
	}




	cross_t += clock() - cross_t_temp;
	individual out_ind(child_code, fcode);			//output individual
	if (TGM == true)
	{
		out_ind.TGM_fitness = ind1.TGM_fitness;
		out_ind.TGM_fitness.insert(out_ind.TGM_fitness.begin(), std::vector<double>(fcode.Objs(), 0));
		out_ind.TGM_fitness.pop_back();

		if (out_ind.TGM_fitness.size() != TGM_size + 1)
			abort();
	}
	return out_ind;
}
void MOEAD::LL_Mutation(individual &child, int iGen, GA_parameters & gapara, function & fcode)
{
	if (ENCODING != "Real")
		abort();

	double gen_val;

	if (T_con == "nfes")
	{
		gen_val = nfes / max_iterations;

	}
	else
	{
		gen_val = iGen/max_generations;
	}

	double rate = gapara.Mut_Prob();
	//copy the boundary vector
	std::vector<STRUCTURES::boundaries> var_bound = fcode.Bound();


	double prm = 0.05;
	int i, mark = 1;
	double a = pow(1.0 - 1.0 * gen_val, 0.7);
	double rm, y;
	for (i = 0; i < nvar; i++) {
		if (Random() <= rate) {
			mark = 0; // at least one component of variable does mutation;
			rm = prm* (Random() - 0.5)*(1 - pow(Random(), -a));
			y = child.Code_Show(i) + rm * (var_bound[i].upper - var_bound[i].lower);

			if (y < var_bound[i].lower)
				y = var_bound[i].lower + 0.5 * Random()*(child.Code_Show(i) - var_bound[i].lower);
			if (y > var_bound[i].upper)
				y = var_bound[i].upper - 0.5 * Random()*(var_bound[i].upper - child.Code_Show(i));

			child.Code_Set(i) = y;
		}
	}
	if (mark) {
		i = int(nvar*Random());
		rm = prm * (Random() - 0.5)*(1 - pow(Random(), -a));
		y = child.Code_Show(i) + rm * (var_bound[i].upper - var_bound[i].lower);

		if (y < var_bound[i].lower)
			y = var_bound[i].lower + 0.5 * Random()*(child.Code_Show(i) - var_bound[i].lower);
		if (y > var_bound[i].upper)
			y = var_bound[i].upper - 0.5 * Random()*(var_bound[i].upper - child.Code_Show(i));

		child.Code_Set(i) = y;
	}
}

void MOEAD::Permutation(std::vector<int> &perm)
{
	unsigned int i, j, size = perm.size();
	std::vector<int> index(size);
	for (i = 0; i<size; i++) index[i] = i;
	for (i = 0; i<size; i++)
	{
		j = int(size*Random());
		while (1)
		{
			if (index[j] >= 0)
			{
				perm[i] = index[j];
				index[j] = -1;
				break;
			}
			else if (j == size - 1)
			{
				j = 0;
			}
			else
			{
				j++;
			}
		}
	}
}



