#include "UNSGAIII.h"
#include "Define.h"
#include "Const.h"
#include "Random.h"
#include "Support_Functions.h"
#include <ctime>

//UNSGAIII have to give output as a vector of individuals or population (vector of individuals should be easier)

//!!!!! CHeck if changing saving to after elitism will give different results (the same with PF checking)
//!!!! write (if was erased) and copy original parent function



extern time_t selec_t;										//Time of Selection
extern time_t mut_t;					//Time of mutation
extern time_t col_t;					//time of collective
extern time_t elit_t;					//Time of UNSGAIII
//extern short MLSt;						//MLSt type

namespace UNSGAIII 
{
	void NSGAII_Selection(std::vector<individual> & mix_pop, std::vector<individual> & parent_pop);
	void NSGAIII_Selection(std::vector<individual> & mix_pop, std::vector<individual> & parent_pop);
	std::vector<short> bos(std::vector<individual> &mix_pop, std::vector<short> feasible_index);

	struct CV_Val_Index
	{
		short index;
		float CV_val;
	};
}
/*
*Calculate new collective using UNSGAIII algorithm*
@param Col - address of a given collective
*/
void UNSGAIII::UNSGAIII_Calc(collective & col)
{
	//get the number of objectives
	short nobj = col.FCode_Show()[0].Objs();


	//Copy the parent vector
	std::vector<individual> parent_pop = col.Indiv_Show();		//parent population



	//Check if collective need initial crowding distance calculation and ranking
	if (col.Was_Erased())
	{
		time_t elite_t_temp = clock();				//starting time of whole operation

		//calculate crowding distance and ranking
		if (nobj == 2) //Do NSGAII
			Rank_Crowding_Distance_Assign(parent_pop, col.fit_index[0]);
		else //Do NSGAIII
			;

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
	if (nobj == 2) //Do NSGAII
		Nondominated_Sort_Fill(mix_pop, parent_pop, col.fit_index[0]);
	else //Do NSGAIII
		;

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

void UNSGAIII::NSGAII_Selection(std::vector<individual> & mix_pop, std::vector<individual> & parent_pop)
{

}
void UNSGAIII::NSGAIII_Selection(std::vector<individual> & mix_pop, std::vector<individual> & parent_pop)
{
	short nobj = mix_pop[0].Fitness_Show().size();
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
		std::vector<short> Ranks = bos(mix_pop, feasible_index);

		for (int ix = 0; ix < feasible_index.size(); ix++)
		{
			Fronts[Ranks[ix]-1].push_back(feasible_index[ix]);
		}

	}

	std::vector<float> pop2Dir(pop_size, 0);
	std::vector<float> pop2DirDistances(pop_size, 0);
	bool associationsready = false;

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
		associationsready = true;

		std::vector<short> cd(mix_pop_size,0);

		cd[0] = Fronts[0].size();
		for (int ix = 1; ix < mix_pop_size; ix++)
		{
			cd[ix] = cd[ix - 1] + Fronts[ix].size();
		}
		short lastfront = 0;
		for (int ix = 1; ix < cd.size(); ix++)
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

		std::vector<std::vector<double>> combinedobj(nobj, std::vector<double>(combinepopindex.size(), 0.));

		for (int ix = 0; ix < combined_size; ix++)
		{
			std::vector<double> fit_vect = mix_pop[combinepopindex[ix]].Fitness_Show();
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



	}
}