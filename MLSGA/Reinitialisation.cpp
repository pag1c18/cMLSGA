#include "Reinitialisation.h"
#include "SVM.h"
#include "Clustering.h"
#include "Random_Label.h"
#include "Define.h"
#include "MTS.h"
#include "MOEAD.h"
#include "BCE.h"
#include "HEIA.h"
#include <ctime>
#include "Support_Functions.h"

extern short n_col;
extern time_t SVM_t;										//Time of SVM
extern int nfes;

void Reinit(std::vector<collective> & col, pareto_front & PF, std::string & re_mode, std::string & mode, std::vector<std::ofstream> & graph_v_vector, std::vector<std::vector<collective>> & col_memory)
{

	if (ONE_OBJ_OVERRIDE != true)
		//calculate new fitness of the old PF
		PF.Fitness_Calc();

	//Reinitialize according to mode chosen
	if (re_mode == "None")
	{
		for (short iCol = 0; iCol < col.size(); iCol++)
		{
			//Get the index of the current collective
			short col_index = col[iCol].Index_Show() - 1;

			if (VIDEO == true && FITNESS_ALL == true)
			{
				extern std::ofstream *graph_v;					//extern file outputs - MLSGA.cpp
																//Set the pointer to the current video data file
				graph_v = &graph_v_vector[col_index];
			}
			//recalculate the fitness of the individuals in the collective
			col[iCol].population::Fitness_Calc();

			//recalculate the fitness of the collective
			col[iCol].Fitness_Calc();

			//Check new PF
			//Check if multi objective optimisation
			if (ONE_OBJ_OVERRIDE != true)
			{
				//Search for pareto front
				PF.Pareto_Search(col[iCol]);		//PF search
			}

			//For MLSGA create new elite
			if (mode == "Normal")
				col[iCol].Elite_Create();
			else
				col[iCol].save();

			if (mode == "BCE")
				BCE::BCE_Time_Update(col[0].FCode_Show()[0]);
			else if (mode == "HEIA")
				HEIA::HEIA_Time_Update(col[0].FCode_Show()[0]);

			nfes += col[iCol].Size_Show();
		}
	}
	else if (re_mode == "Random")
	{
		abort();
		/***** GA Initialisation *****/
		/*int psize = col[0].GAPara_Show()[0].Pop_Size();
		//Initialise the population
		population p1{ col[0].FCode_Show()[0],col[0].GAPara_Show()[0],col[0].MCode_Show()[0],col[0].SCode_Show()[0],col[0].CCode_Show()[0], psize };
		//Clustering(p1);

		if (mode == "MOEAD" || mode == "MOEADMSF" || mode == "MOEADPSF" || mode == "MOEADM2M")
		{
			//initialise the MOEAD
			MOEAD::MOEAD_Init(col[0].FCode_Show()[0].Objs(), p1);
		}
		else if (mode == "MTS")
		{
			//Initialise MTS
			MTS::MTS_Init(p1, PF);
		}
		else if (mode == "BCE")
		{
			BCE::BCE_Init(col[0].FCode_Show()[0].Objs(), n_col, p1);
		}
		else if (mode == "HEIA")
		{
			HEIA::HEIA_Init(n_col, p1);
		}

		//Initialise the vector for storaging real labels
		//std::vector<int> real_label(GA_standard.Pop_Size(), 0);										//Real labels for population separation
		std::vector<short> real_label;
		//Get the starting time of the SVM
		time_t SVM_t_temp = clock();

		//Do SVM
		//If only 1 collective SVM is pointless
		if (n_col > 1)
		{
			if (GROUPING == "SVM")
				real_label = SVM(p1);
			else if (GROUPING == "Cluster")
				real_label = Clustering(p1);
			else if (GROUPING == "Random")
				real_label = Random_Label(p1);
			else
				abort();
			//real_label = SVM_LABEL(p1);
			//real_label = Random_Label(p1);
		}
		else if (n_col < 0)
		{
			std::cout << "**********************************\n"
				<< "Wrong paramters chosen. Number of collective cannot be negative. Check const.h\nProgram will terminate"
				<< "\n**********************************\n";
			system("pause");
			abort();
		}
		else if (n_col == 1)
		{
			real_label = std::vector<short>(psize, 1);
		}
		else if (n_col == 0 && NSGAII_OVERRIDE == false)
		{
			std::cout << "No collectives. Do You want to proceed?";
			system("pause");
		}

		//Calculate the time of SVM
		SVM_t += clock() - SVM_t_temp;
		*/
		/***** Collevtive generation *****/
		/*
		//Initialise the vector of collectives
		std::vector<collective> col_temp;																	//vector for storage of collectives
																											//int sum_c = 0;


		time_t col_t_temp = clock();


		std::vector<short> biggest_col{ 0,0,0,0 };														//Vector for storing which collectives are the biggest ones
																										//graph_v_vector.clear();
		int pop_size_temp = 0;
		for (short iCol = 0; iCol < n_col; iCol++)														//Loop for collective generation
		{
			//Assign the individuals to the collective
			col_temp.push_back(collective::collective(p1, real_label, iCol + 1));

			//if (VIDEO == true && FITNESS_ALL == true)
			//Open the data file for the first frame of the video
			//graph_v_vector.push_back(ofstream("Temp/" + Name_Get("graph_v_" + std::to_string(index_r) + "_c" + std::to_string(iCol + 1), 1)));

			//sum_c += col_temp[iCol].Size_Show();

			//Find the biggest collective
			if (n_col != 1)
			{
				if (col_temp[iCol].Size_Show() > biggest_col[1])
				{
					biggest_col[2] = biggest_col[0];
					biggest_col[3] = biggest_col[1];
					biggest_col[0] = iCol;
					biggest_col[1] = col_temp[iCol].Size_Show();
				}
			}
		}
		for (short iCol = 0; iCol < n_col; iCol++)
		{
			//Get the size of the collective
			int temp_i = col_temp[iCol].Size_Show();

			if (n_col != 1)
			{
				//Collective cannot be empty or too small
				//If is too small

				if (temp_i <= (psize / min_col_size_multi) || temp_i <= 2)
				{
					//Check which collective is the biggest
					if (biggest_col[1] >= biggest_col[3])
					{
						//Cut the individuals from the biggest collective to the new one
						for (int i = temp_i; i < (int)(psize / new_col_size_multi); i++)
						{
							//Get the individual to be cut
							individual temp_indi = col_temp[biggest_col[0]].Indiv_Show(0);		//Inidividual which will be copied

																								//Add the individual to the new collective
							col_temp[iCol].Add(temp_indi);

							//Remove it from the old one
							col_temp[biggest_col[0]].Remove();
						}
						//Correct the size of the biggest collective
						biggest_col[1] = col_temp[biggest_col[0]].Size_Show();
					}
					else
					{

						for (int i = temp_i; i < (int)(psize / new_col_size_multi); i++)
						{
							individual temp_indi = col_temp[biggest_col[2]].Indiv_Show(0);
							col_temp[iCol].Add(temp_indi);
							col_temp[biggest_col[2]].Remove();
						}
						biggest_col[3] = col_temp[biggest_col[2]].Size_Show();
					}
				}
			};
		}
		col = col_temp;
		for (short iCol = 0; iCol < n_col; iCol++)
		{
			if (VIDEO == true && FITNESS_ALL == true)
			{
				extern std::ofstream *graph_v;					//extern file outputs - MLSGA.cpp
				graph_v = &graph_v_vector[iCol];
			}

			if (mode != "MOEAD" && mode != "PAES" && mode != "MTS" || mode != "MOEADMSF" || mode != "MOEADPSF" || mode != "MOEADM2M")
			{
				//Calculate the fitness of individuals in the collective
				col[iCol].population::Fitness_Calc();
				if (mode == "NSGAII" || mode == "DMOEADD")
					col[iCol].save();
			}
			else if (mode == "MOEAD" || mode == "MOEADMSF" || mode == "MOEADPSF" || mode == "MOEADM2M")
				//initialise the population for MOEAD
				MOEAD::Population_Init(col[0].FCode_Show()[0].Objs(), col[iCol]);
			else if (mode == "MTS")
				//initialise the collective parameters for MTS
			{
			MTS::MTS_Init_Col(col[iCol]);
			col[iCol].population::Fitness_Calc();
			nfes += col[iCol].Size_Show();
			}
			else if (mode == "BCE")
			{
				BCE::BCE_Pop_Init(col[0].FCode_Show()[0].Objs(), col[iCol]);
			}
			else if (mode == "HEIA")
			{
				HEIA::HEIA_Pop_Init(col[iCol]);
			}


			if (ONE_OBJ_OVERRIDE != true)
			{
				//Find the PF
				PF.Pareto_Search(col[iCol]);			//PF search
			}
		}
		//Check if population size match
		for (short iCol = 0; iCol < n_col; iCol++)
			pop_size_temp += col[iCol].Size_Show();
		if (pop_size_temp != psize)
			abort();
				*/
	}
	else
	{
		//Copy the amount of variables and obj
		int var_size = col[0].FCode_Show()[0].Vars();		//amount of variables
		int obj_size = col[0].FCode_Show()[0].Objs();		//amount of objectives

		//Copy the boundaries vector
		std::vector<double> bound_up(var_size,0);		//vector of upper boundaries
		std::vector<double> bound_l(var_size, 0);		//vector of lower boundaries
		for (int iVar = 0; iVar < var_size; iVar++)
		{
			bound_up[iVar] = col[0].FCode_Show()[0].Bound(iVar, "upper");
			bound_l[iVar] = col[0].FCode_Show()[0].Bound(iVar, "lower");
		}
		for (int iCol = 0; iCol < col.size(); iCol++)
		{
			//get the size of the current collective
			int csize = col[iCol].Size_Show();		//current collective size
			
			//Find the indexes of matching collectives in the memory
			std::vector<int> index_col_gen(gen_memory_size, 0);
			if (col.size() > 1)
			{
				
				for (int iMem = 0; iMem < gen_memory_size; iMem++)
				{
					for (int iCol_mem = 0; iCol_mem < col.size(); iCol_mem++)
					{
						if (col[iCol].Mode_Show() == col_memory[iMem][iCol_mem].Mode_Show())
						{
							if (col[iCol].Index_Show() == col_memory[iMem][iCol_mem].Index_Show())
							{
								index_col_gen[iMem] = iCol_mem;
								break;
							}
						}
					}
				}
			}

			//Reinitialise every individual
			for (int iInd = 0; iInd < csize; iInd++)
			{
				//copy the variables and objectives of thecurrent individual
				std::vector<double> ind_obj = col[iCol].Indiv_Show(iInd).Fitness_Show();	//vector with fitness values of current individual
				std::vector<double> ind_var = col[iCol].Indiv_Show(iInd).Code_Show(); //vector with variables of current individual
				
				std::vector<int> c_index(gen_memory_size, 0);		//index of the clostest solutions for each individual in the population in memory slots
				for (int iMem = 0; iMem < gen_memory_size; iMem++)
				{
					double dist_temp, dist_min;			//values of current distance and the current min one
					//Calculate the distance and find the min one
					for (int iInd2 = 0; iInd2 < csize; iInd2++)
					{
						//According to the chosen method
						if (re_mode == "BR" || re_mode == "VP") //based on PS distance
						{
							//calculate the distance
							if (iMem == 0)
								dist_temp = Distance<double>(ind_var, col_memory[iMem][index_col_gen[iMem]].Indiv_Show(iInd2).Code_Show());
							else
								dist_temp = Distance<double>(col_memory[iMem - 1][index_col_gen[iMem-1]].Indiv_Show(c_index[iMem-1]).Code_Show(), col_memory[iMem][index_col_gen[iMem]].Indiv_Show(iInd2).Code_Show());
						}
						else if (re_mode == "CER")	//based on PF distance
						{
							//calculate the distance
							if (iMem == 0)
								dist_temp = Distance<double>(ind_obj, col_memory[iMem][index_col_gen[iMem]].Indiv_Show(iInd2).Fitness_Show());
							else
								dist_temp = Distance<double>(col_memory[iMem - 1][index_col_gen[iMem-1]].Indiv_Show(c_index[iMem-1]).Fitness_Show(), col_memory[iMem][index_col_gen[iMem]].Indiv_Show(iInd2).Fitness_Show());
						}
						//update the min value and index of min individual
						if (iInd2 == 0)
							dist_min = dist_temp;
						else if (dist_min > dist_temp)
						{
							dist_min = dist_temp;
							c_index[iMem] = iInd2;
						}
					}
				}

				std::vector<double> var_new;		//vector for storage of the new variables/code
			
				//store the variables of the clostest memory individuals
				std::vector <std::vector<double>> ind_past_var;
				for (int iMem = 0; iMem < gen_memory_size; iMem++)
				{
					ind_past_var.push_back(col_memory[iMem][index_col_gen[iMem]].Indiv_Show(c_index[iMem]).Code_Show());
				}

				//Reinitialise the individual based on the chosen method
				if (re_mode == "BR" || re_mode == "CER")
				{
					//Calculate the new variables for the individual
					for (int iVar = 0; iVar < var_size; iVar++)
					{
						//initialise the temp values
						double m, P;
						double K = 0;
						//Calculate the temp values						
						for (int iMem = 0; iMem < gen_memory_size; iMem++)
						{
							if (iMem == 0)
							{
								P = ind_var[iVar] - ind_past_var[iMem][iVar];
								K = (gen_memory_size - 1)*abs(P);
							}
							else
							{
								P = ind_past_var[iMem - 1][iVar] - ind_past_var[iMem][iVar];
								K -= abs(P);
							}
						}
						//Caclulate the offset
						if (Random() < 0.5)
							m = 1 + tanh(K);
						else
							//m = 1 + sgn(K)*N_Distribution(1, abs(K));
							m = 1 + K*abs(N_Distribution(0, 1));

						//Calculate the new variable and push it to the vector
						double var_temp = ind_var[iVar] + m*(ind_var[iVar] - ind_past_var[0][iVar]);
						
						
						//Check if the new variable is in the boundaries
						if (var_temp > bound_up[iVar] || var_temp < bound_l[iVar])
							//if not push the old value
							var_new.push_back(ind_var[iVar]);
						else
							//push the new one
							var_new.push_back(var_temp);
					}

				}
				else if (re_mode == "VP")
				{
					//calculate the standard deviation for predicted noise
					double std_dev;

					std_dev = sqrt(1. / (4 * var_size)*pow(Distance(ind_var, ind_past_var[0]),2));

					if (std_dev == 0)
						std_dev = 0.1;
					//Choose the method for the reinitialisation by random, 0- variation, 1 - prediction
					int VP_method = 0;
					if (Random() < 0.5)
						VP_method = 1;

					//Calculate the new variables for the individual
					for (int iVar = 0; iVar < var_size; iVar++)
					{
						//Calculate the noise
						double eps = N_Distribution(0, std_dev);
						
						//Calculate the new variable and push it to the vector
						double var_temp;		//temporary new variables
					
						if (VP_method == 0)  //Variation
							var_temp = ind_var[iVar] + eps;
						else //Prediction
							var_temp = ind_var[iVar] + (ind_var[iVar] - ind_past_var[0][iVar]) + eps;

						//Check if the new variable is in the boundaries
						if (var_temp > bound_up[iVar] || var_temp < bound_l[iVar])
							//if not push the old value
							var_new.push_back(ind_var[iVar]);
						else
							//push the new one
							var_new.push_back(var_temp);
					}
				}

				//Check the size of the new code/variables
				if (var_new.size() != var_size)
					abort();
				else
					//update the individual code
					col[iCol].Indiv_Set(iInd).Code_Set() = var_new;
			}
			//Get the index of the current collective
			short col_index = col[iCol].Index_Show() - 1;

			if (VIDEO == true && FITNESS_ALL == true)
			{
				extern std::ofstream *graph_v;					//extern file outputs - MLSGA.cpp
																//Set the pointer to the current video data file
				graph_v = &graph_v_vector[col_index];
			}
			//recalculate the fitness of the individuals in the collective
			col[iCol].population::Fitness_Calc();

			//recalculate the fitness of the collective
			col[iCol].Fitness_Calc();

			nfes += col[iCol].Size_Show();

			//Check new PF
			//Check if multi objective optimisation
			if (ONE_OBJ_OVERRIDE != true)
			{
				//Search for pareto front
				PF.Pareto_Search(col[iCol]);		//PF search
			}

			//For MLSGA create new elite
			if (mode == "Normal")
				col[iCol].Elite_Create();
			else
				col[iCol].save();

						
			if (mode == "BCE")
				BCE::BCE_Time_Update(col[0].FCode_Show()[0]);
			else if (mode == "HEIA")
				HEIA::HEIA_Time_Update(col[0].FCode_Show()[0]);
			
		}
		
		//update the memory of the past collectives
		for (int i = gen_memory_size - 1; i > 0; i--)
		{
			col_memory[i] = col_memory[i - 1];
		}
		col_memory[0] = col;

	}


}