#define _CRT_SECURE_NO_WARNINGS

#include "PAES.h"
#include "Sobol.h"
#include "Support_Functions.h"
#include "Selection.h"
#include <ctime>

extern int nfes;						//number of function evaluations
extern int nweightindex;				//count of wieght indexes used
								//std::vector<individual> saved;	//vector of saved individuals


extern time_t selec_t;										//Time of Selection
extern time_t mut_t;					//Time of mutation
extern time_t col_t;					//time of collective
extern time_t elit_t;					//Time of 
extern time_t cross_t;					//Time of crossover

extern std::vector<std::vector<double>> weight_vector;		//Storage of weight vectors
extern std::vector<double> idealpoint;			//Storage of ideal point


/*
Calculate individuals using PAES algorithm
@param Col - address of a given collective
@param mode - mode which will be used: 1 - choose the best individual and replace the worst indiviiduals; 2 - for each individual; 3 - the same as 1 but with weight vectors; 4 - the same as 2 but with weight vectors; 5 - as 1 but with MLSt selection; 6 - as 2 but with MLSt selection; 7 - with weight vectors and neighbourhood by random; 8 - as 7 but by selection
*/
std::vector<individual> PAES_Calc(collective & col, int mode, selection<individual> & scode)
{
	//Clear the old PF storage
	std::vector<individual> PF_storage;			//output

	//Copy the parent vector
	std::vector<individual> parent_pop = col.Indiv_Show();		//parent population

	//Calculate PAES limit
	int PAES_limit = col.Size_Show() * MOEAD_niche_multi;
	
	//Check if is not too small
	if (PAES_limit < 4)
		PAES_limit = 4;

	//Calculate PAES replace limit
	int PAES_replace_limit = PAES_limit / 4;
	if (PAES_replace_limit <2)
		PAES_replace_limit = 2;

	//Get the size of the collective
	int csize = col.Size_Show();

	//Do PAES accordding to current mode
	if (mode == 1)
	{
		time_t select_t_temp = clock();
		//Sort the individuals
		col.Sort_Individuals();

		//Select the individuals
		std::vector<individual> selected;		//vector of selected indivudals
		
		selected = scode.Select(col.Indiv_Show(), PAES_limit, col.FCode_Show()[0].Cons(), col.fit_index[0]);
		selec_t += clock() - select_t_temp;

		//Mutate individuals
		for (int i = 0; i < PAES_limit; i++)
		{
			time_t mut_t_temp = clock();

			//Mutate the individual
			col.Mutation(selected[i]);

			mut_t += clock() - mut_t_temp;

			time_t col_t_temp = clock();
			//Calculate the fitness of the mutated individual
			selected[i].Fitness_Calc(col.FCode_Show()[0]);

			col_t += clock() - col_t_temp;

			nfes++;

			time_t elit_t_temp = clock();

			int PAES_replace_limit_temp = 0;			//Count of how many solutions have been replaced
														//Replace the wost solutions
			for (int j = 1; j <= csize/2; j++)
			{

				if (PAES_replace_limit_temp >= PAES_replace_limit)
					break;

				//Check dominance
				int flag = Dominance_Check(col.Indiv_Show(csize - j), selected[i], col.fit_index[0]);

				//The new solution is better
				if (flag = -1)
				{
					col.Indiv_Set(csize - j) = selected[i];
					PAES_replace_limit_temp++;
				}
				//Are the same
				else if (flag == 0)
				{
					//Choose solution by random
					float rand = Random_F();
					if (rand <= 0.5f)
					{
						col.Indiv_Set(csize - j) = selected[i];
						PAES_replace_limit_temp++;
					}
				}
			}

			elit_t += clock() - elit_t_temp;

		}
		PF_storage = selected;

	}
	else if (mode == 2)
	{
		//Get the size of the collective
		int csize = col.Size_Show();

		time_t mut_t_temp = clock();
		//Mutate the whole population
		col.Mutation();

		//Calcualte the fitness of the new population
		col.population::Fitness_Calc();

		PF_storage = col.Indiv_Show();

		mut_t += clock() - mut_t_temp;

		time_t elit_t_temp = clock();

		for (int i = 0; i < csize; i++)
		{

			//Check dominance
			int flag = Dominance_Check(col.Indiv_Show(i), parent_pop[i], col.fit_index[0]);
			
			//The new solution is worse
			if (flag == -1)
				//Update
				col.Indiv_Set(i) = parent_pop[i];
			//Are the same
			else if (flag == 0)
			{
				//Choose solution by random
				float rand = Random_F();
				if (rand <=0.5f)
					col.Indiv_Set(i) = parent_pop[i];
			}
			nfes++;
		}
		elit_t += clock() - elit_t_temp;
	}
	else if (mode == 3)
	{
		time_t select_t_temp = clock();
		//Sort the individuals
		col.Sort_Individuals();

		//Select the individuals
		std::vector<individual> selected;		//vector of selected indivudals
		selected = scode.Select(col.Indiv_Show(), PAES_limit, col.FCode_Show()[0].Cons(), col.fit_index[0]);

		selec_t += clock() - select_t_temp;
		//Mutate individuals
		for (int i = 0; i < PAES_limit; i++)
		{
			time_t mut_t_temp = clock();

			//Mutate the individual
			col.Mutation(selected[i]);

			mut_t += clock() - mut_t_temp;

			time_t col_t_temp = clock();
			//Calculate the fitness of the mutated individual
			selected[i].Fitness_Calc(col.FCode_Show()[0]);

			col_t += clock() - col_t_temp;

			nfes++;

			time_t elit_t_temp = clock();

			int PAES_replace_limit_temp = 0;			//Count of how many solutions have been replaced
														//Replace the wost solutions
			for (int j = 1; j <= csize / 2; j++)
			{

				if (PAES_replace_limit_temp >= PAES_replace_limit)
					break;

				//Check dominance
				int flag = Dominance_Check3(col.Indiv_Show(csize - j), selected[i]);

				//The new solution is better
				if (flag = -1)
				{
					col.Indiv_Set(csize - j) = selected[i];
					PAES_replace_limit_temp++;
				}
				//Are the same
				else if (flag == 0)
				{
					//Choose solution by random
					float rand = Random_F();
					if (rand <= 0.5f)
					{
						col.Indiv_Set(csize - j) = selected[i];
						PAES_replace_limit_temp++;
					}
				}
			}
			elit_t += clock() - elit_t_temp;

		}
		PF_storage = selected;
	}
	else if (mode == 4)
	{
		//Get the size of the collective
		int csize = col.Size_Show();

		time_t mut_t_temp = clock();
		//Mutate the whole population
		col.Mutation();

		//Calcualte the fitness of the new population
		col.population::Fitness_Calc();

		PF_storage = col.Indiv_Show();

		mut_t += clock() - mut_t_temp;

		time_t elit_t_temp = clock();

		for (int i = 0; i < csize; i++)
		{

			//Check dominance
			int flag = Dominance_Check3(col.Indiv_Show(i), parent_pop[i]);

			//The new solution is worse
			if (flag == -1)
				//Update
				col.Indiv_Set(i) = parent_pop[i];
			//Are the same
			else if (flag == 0)
			{
				//Choose solution by random
				float rand = Random_F();
				if (rand <= 0.5f)
					col.Indiv_Set(i) = parent_pop[i];
			}
			nfes++;
		}
		elit_t += clock() - elit_t_temp;
	}
	else if (mode == 5)
	{
		time_t select_t_temp = clock();
		//Sort the individuals
		col.Sort_Individuals();

		//Select the individuals
		std::vector<individual> selected;		//vector of selected indivudals
		selected = scode.Select(col.Indiv_Show(), PAES_limit, col.FCode_Show()[0].Cons(), col.fit_index[0]);

		selec_t += clock() - select_t_temp;
		//Mutate individuals
		for (int i = 0; i < PAES_limit; i++)
		{
			time_t mut_t_temp = clock();

			//Mutate the individual
			col.Mutation(selected[i]);

			mut_t += clock() - mut_t_temp;

			time_t col_t_temp = clock();
			//Calculate the fitness of the mutated individual
			selected[i].Fitness_Calc(col.FCode_Show()[0]);

			col_t += clock() - col_t_temp;

			nfes++;

			time_t elit_t_temp = clock();

			int PAES_replace_limit_temp = 0;			//Count of how many solutions have been replaced
														//Replace the wost solutions
			for (int j = 1; j <= csize / 2; j++)
			{

				if (PAES_replace_limit_temp >= PAES_replace_limit)
					break;

				//Check dominance
				int flag = Dominance_Check2(col.Indiv_Show(csize - j), selected[i], col.Index_Show(), col.fit_index[0]);

				//The new solution is better
				if (flag = -1)
				{
					col.Indiv_Set(csize - j) = selected[i];
					PAES_replace_limit_temp++;
				}
				//Are the same
				else if (flag == 0)
				{
					//Choose solution by random
					float rand = Random_F();
					if (rand <= 0.5f)
					{
						col.Indiv_Set(csize - j) = selected[i];
						PAES_replace_limit_temp++;
					}
				}
			}
			elit_t += clock() - elit_t_temp;
		}
		PF_storage = selected;
	}
	else if (mode == 6)
	{
		//Get the size of the collective
		int csize = col.Size_Show();

		time_t mut_t_temp = clock();
		//Mutate the whole population
		col.Mutation();

		//Calcualte the fitness of the new population
		col.population::Fitness_Calc();

		PF_storage = col.Indiv_Show();

		mut_t += clock() - mut_t_temp;

		time_t elit_t_temp = clock();

		for (int i = 0; i < csize; i++)
		{

			//Check dominance
			int flag = Dominance_Check2(col.Indiv_Show(i), parent_pop[i], col.Index_Show(), col.fit_index[0]);

			//The new solution is worse
			if (flag == -1)
				//Update
				col.Indiv_Set(i) = parent_pop[i];
			//Are the same
			else if (flag == 0)
			{
				//Choose solution by random
				float rand = Random_F();
				if (rand <= 0.5f)
					col.Indiv_Set(i) = parent_pop[i];
			}
			nfes++;
		}
		elit_t += clock() - elit_t_temp;
	}
	else if (mode == 7)
	{
		time_t select_t_temp = clock();
		//Sort the individuals
		col.Sort_Individuals_NAMDA();

		//Select individuals by random
		std::vector<individual> selected;		//vector of selected indivudals
		std::vector<int> selected_ix;		//vector of indexes of selected individuals
		for (int i = 0; i < PAES_limit; i++)
		{
			selected_ix.push_back(Random_I(0, csize - 1));
			selected.push_back(col.Indiv_Show(selected_ix[i]));
		}

		for (int i = 0; i < PAES_limit; i++)
		{
			
			time_t mut_t_temp = clock();

			//Mutate the individual
			col.Mutation(selected[i]);

			mut_t += clock() - mut_t_temp;

			time_t col_t_temp = clock();
			//Calculate the fitness of the mutated individual
			selected[i].Fitness_Calc(col.FCode_Show()[0]);

			col_t += clock() - col_t_temp;

			nfes++;

			time_t elit_t_temp = clock();

			//Replace the wost solutions
			for (int j = selected_ix[i] - (PAES_limit/2); j < selected_ix[i] + (PAES_limit / 2); j++)
			{
				//Check if first index is not too small or the last one is not too big
				if (j < 0)
					j = 0;
				else if (j == csize)
					break;
				else if (j == selected_ix[i])
					continue;


				//Check dominance
				int flag = Dominance_Check3(col.Indiv_Show(j), selected[i]);

				//The new solution is better
				if (flag = -1)
				{
					col.Indiv_Set(j) = selected[i];
				}
				//Are the same
				else if (flag == 0)
				{
					//Choose solution by random
					float rand = Random_F();
					if (rand <= 0.5f)
					{
						col.Indiv_Set(j) = selected[i];
					}
				}
			}
			elit_t += clock() - elit_t_temp;
		}
		PF_storage = selected;
	}
	else if (mode == 8)
	{
		time_t select_t_temp = clock();
		//Sort the individuals
		col.Sort_Individuals_NAMDA();

		//Select individuals by random
		std::vector<individual> selected;		//vector of selected indivudals
		std::vector<int> selected_ix;		//vector of indexes of selected individuals
		selected = scode.Select(col.Indiv_Show(), PAES_limit, col.FCode_Show()[0].Cons(), col.fit_index[0]);

		//Copy the indexes
		for (int i = 0; i < PAES_limit; i++)
		{
			for (int j = 0; j < csize; j++)
			{
				if (selected[i].namda == col.Indiv_Show(j).namda)
					selected_ix.push_back(j);
			}
		}

		for (int i = 0; i < PAES_limit; i++)
		{

			time_t mut_t_temp = clock();

			//Mutate the individual
			col.Mutation(selected[i]);

			mut_t += clock() - mut_t_temp;

			time_t col_t_temp = clock();
			//Calculate the fitness of the mutated individual
			selected[i].Fitness_Calc(col.FCode_Show()[0]);

			col_t += clock() - col_t_temp;

			nfes++;

			time_t elit_t_temp = clock();

			//Replace the wost solutions
			for (int j = selected_ix[i] - (PAES_limit / 2); j < selected_ix[i] + (PAES_limit / 2); j++)
			{
				//Check if first index is not too small or the last one is not too big
				if (j < 0)
					j = 0;
				else if (j == csize)
					break;
				else if (j == selected_ix[i])
					continue;


				//Check dominance
				int flag = Dominance_Check3(col.Indiv_Show(j), selected[i]);

				//The new solution is better
				if (flag = -1)
				{
					col.Indiv_Set(j) = selected[i];
				}
				//Are the same
				else if (flag == 0)
				{
					//Choose solution by random
					float rand = Random_F();
					if (rand <= 0.5f)
					{
						col.Indiv_Set(j) = selected[i];
					}
				}
			}
			elit_t += clock() - elit_t_temp;
		}
		PF_storage = selected;
	}
	else
	{
		abort();
	}

	//save the collective
	col.save();

	time_t col_t_temp = clock();				//Starting time of collective operations
	//calculate the fitness of the collective
	col.Fitness_Calc();

	col_t += clock() - col_t_temp;



	return PF_storage;
}

/*
*Initialise the PAES parameters*
@param n_obj - number of objectives
@param pop - population used for initialisation
*/
void PAES::PAES_Init(short n_obj, population & pop)
{
	time_t elite_t_temp = clock();				//starting time of whole operation

	//Erase the number of evaluations
	nfes = 0;
	idealpoint.clear();
	idealpoint = std::vector<double>(n_obj, 1.0e+30);


	//Calculate fitness
	pop.Fitness_Calc();

	//Add the fitness evaluations
	nfes += pop.Size_Show();

	//get namda for whole population only when values are from file
	abort(); //PAES is not working
/*	if (MOEAD_SOBOL == false && MOEAD_SEPARATE_WEIGHTS == false)
	{
		// Read weight vectors from a data file
		char filename[1024];
		sprintf(filename, "Input/MOEADWeight/W%dD_%d.dat", n_obj, pop.Size_Show());
		std::ifstream readf(filename);


		for (int i = 0; i < pop_size; i++)
		{
			std::vector<double> temp_namda = std::vector<double>(n_obj, 0.);
			pop.Indiv_Set(i).namda = temp_namda;
			// Load weight vectors
			for (int j = 0; j < n_obj; j++)
			{
				double namda_temp;
				readf >> pop.Indiv_Set(i).namda[j];
				//printf("%f ", sub.namda[j]);
			}
		}
	}
	else if (n_obj != 2)
	{
		abort();
	}
	else if (MOEAD_SEPARATE_WEIGHTS == true)
	{
		Separate_Weights_Generate(pop.Size_Show());
		nweightindex = 0;
	}*/
	//calculate time
	elit_t += clock() - elite_t_temp;
}

/*
*Initialise the PAES population*
@param n_obj - number of objectives
@param col - collective used for initialisation
*/
void PAES::PAES_Population_Init(short n_obj, collective & col)
{
	time_t elite_t_temp = clock();				//starting time of whole operation

	//Get the size of the collective
	int csize = col.Size_Show();

	//Get namda for collective when calculated by sobol
	abort(); //PAES is not working
	/*if (MOEAD_SOBOL == true)
	{
		//only for 2 obj - for now
		if (n_obj != 2)
			abort();

		//get sobol sequence vector
		weight_vector = Sobol_Sequence(csize, n_obj);

		//get weights
		for (int i = 0; i < csize; i++)
		{
			col.Indiv_Set(i).namda = std::vector<double>(n_obj, 0.);
			col.Indiv_Set(i).namda[0] = weight_vector[0][i];
			col.Indiv_Set(i).namda[1] = 1. - weight_vector[0][i];
		}
	}
	//Get namda for collective when calculated by separation
	else if (MOEAD_SEPARATE_WEIGHTS == true)
	{
		//only for 2 obj - for now
		if (n_obj != 2)
			abort();
		//get weights
		for (int i = 0; i < csize; i++)
		{
			col.Indiv_Set(i).namda = weight_vector[nweightindex];

			nweightindex++;
		}
	}*/

	for (int i = 0; i < csize; i++)
	{

		//Save the individual
		col.Indiv_Show(i).save();
	}

	//calculate time
	elit_t += clock() - elite_t_temp;
}




