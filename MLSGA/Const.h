/**
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
along with this program.If not, see < https:///<www.gnu.org/licenses/>. 

	 CONST header
Defining algorithm parameters and constants
	behaviour is set in Define.h

*/

#pragma once

#ifndef CONST_H
#define CONST_H

const double long pi = 3.141592653589793238462643383279502884197l;
const double long e = 2.71828182845904523536028747135266249775724709369995l;

const float Ro_a = 1.2466f;	///<air density for 10C
const float Ro_water = 1025; ///<water density
const double visc_water = 0.0000012; ///<water viscosity

#include <string>
#include <vector>

/*General parameters*
**List of implemented Functions**
*Unconstrained*
1-6 ZDT 1-6
7-13 DTLZ 1-7
14-20 MOP 1-7
21-30 UF 1-10
31-39 WFG 1-9
40-49 IMB 1-10
*Constrained*
50-59 CF 1-10
60-61 DTLZ 8-9
62-65 IMB 11-14
66-74 DAS-CMOP 1-9
*Dynamic Unconstrained*
75-82 UDF 1-9
83-92 JY 1-10
93-97 FDA 1-5
*Dynamic Constrained*
98-112 CDF1-15


*Not Implemented*
5 - ZDT5
JY4
JY9
JY10
FDA4
FDA5
*/




/**MLSGA setting*/
const std::vector<std::vector<std::string>> GA_mode{ {"UNSGAIII"}, {"MOEAD"}, {"MOEADMSF"}, {"MOEADPSF"}, {"BCE"}, {"IBEA"},{"MTS"}, {"HEIA"}, {"MOEADM2M"} }; ///< {"Normal"}, {"UNSGAIII"}, {"MOEAD"}, {"MOEADMSF"}, {"MOEADPSF"}, {"BCE"}, {"IBEA"},{"MTS"}, {"HEIA"}, {"MOEADM2M"}
const bool MLSGA_Hybrid = true;										///<defines if the hybid algorithm is used or the original one. Default false
const bool MLSGA_norm_obj = false;									///<Define if the objectives for fitness evaluations should be normalised (for MSF and PSF, but true makes PSF inoperative). Default 8
const unsigned short MLSGA_n_col_b = 8;								///<Numbers of collectives - begin ***with a choice of 4, 6, 8 (is overriden to 1 if MLSGA_Hybrid == false). Default 8
const unsigned short MLSGA_n_col_e = 8;								///<Numbers of collectives - end - step is always 2 (is overriden to 1 if MLSGA_Hybrid == false). Default 8
const unsigned short MLSGA_n_col_elim = 1;							///<Number of collectives elimiated ***with a choice of 1, 2  (is overriden to 0 if MLSGA_Hybrid == false). Default 1
const unsigned short MLSGA_col_elim_limit = 10;						///<Define how many genrations have to pass before col elimination occurs - if 1 elimiantion will occur every generation - MTS have override that it change to 1 - line 494 in MLSGA.cpp. Default 10
const float MLSGA_elite_size = 0.1f; 								///<How many individuals become elite <0.0-1.0>. Default 0.1
const unsigned short MLSGA_n_MLS_b = 7; 							///<Index of the first MLS type used. Default 7
const unsigned short MLSGA_n_MLS_e = 7; 							///<Index of the last MLS type used. Default 7


/**Parameter setting*/
const int Pop_size_b = 600;											///<Population size - beginning. Default 1000
const int Pop_size_e = 600;											///<Population size - end. Default 1000
const int Pop_size_step = 300;										///<Step for population size change. Default 200
const unsigned short n_func_b = 1;								///<Index of the first function used
const unsigned short n_func_e = 30; 								///<Index of the last function used
const unsigned short n_func_obj = 3;								///<Number of objectives for the scalable problems - DTLZ, WFG
const std::vector<std::string> func_skip{  "Three_obj", "Many_obj" };///<Decides if specific functions will be skipped "Two_obj", "Three_obj", "Many_obj" - DTLZ, WFG
const float mut_prob_min = 0.006f; 									///<Minimal mutation probability - for real coded it is the absolute value, for binary it is divided by number of variables
const float mut_prob_max = 0.006f; 									///<Maximal mutation probability
const float mut_prob_step = 0.02f; 									///<Mutation probability step
const float cross_prob_min = 1;										///<Minimal crossover probability. Default 1
const float cross_prob_max = 1;										///<Maximal crossover probability. Default 1
const float cross_prob_step = 0.10f; 								///<Crossover probability step. Default 0.1


/**Termination conditions - only one is used depending on T_con*/
const unsigned short n_runs = 1;									///<Number of runs. Default 30
const std::string T_con{ "nfes" };									///<Which termination condition is used: "nfes" - iterations number (int max_iterations line 80); "ngen" - number of generations (int max_generations line 81); "ntime" - calculaion time in seconds (int max_time line 82)
const int max_iterations = 300000;//5760000 - 65 - 500;				///<max iterations. Default 300000
const int max_generations = 300;									///<max genenertaions
const int max_time = 120;											///<max calcualtion time

/**Performance evaluation*/
const std::vector<std::string> perf_eval_type{ "IGD", "HV" };		///<Which performance evaluation methods are used "IGD" - inverted generational distance, "HV" - hypervolume

/**Additional benchmarking parameters*/
const unsigned short mut_ix = 1; 									///<Index of the mutation: 1 - polynomial mutation. Default 1
const unsigned short cross_ix = 1; 									///<Index of the crossover: 1 - SBX, 2 - CMRDX, 3 - DE. Default 1
const unsigned short select_ix = 1;									///<Index of the selection (applies only to "normal" MSLGA): 1 - Roulette wheel. Default 1



/**Dynamic functions parameter*/
const int n_steps_dyn = 1;											///<Number of distinct steps; represents the severity of change. Default 5
const int T_dyn = 50;												///<Window where the dynamic problem remains cosntant. Default 5
const std::vector<std::string> Reinit_Mode{ "None" };				///<list of reinitialisation modes that have to be used "None", "BR", "CER", "VP"
const int gen_memory_size = 2;										///<How many generations are saved and used for reinitialisation. Default 2
const int JY_const_window = 400;									///<Window where JY functions remain constatant. Default 300


/**Additional GA parameters*/
const unsigned short Di_c = 20;										///<Distribution index of crossover. Default 20
const unsigned short Di_m = 20;										///<Distribution index of mutation. Default 20
const short fit_index_sel = 2;  									///<Which fitness is used for indiv selection. Default 2
const short fit_index_col_sel = 1; 									///<Which fitness is used for col selection. Default 1
const int Binary_string_size = 1;									///<Lenght of each variable in binary string. Default 1
const int MPCross_size = 5;											///<Defines how many separation points are in multi-point crossover (1 is two point crossver). Default 2

/**GA parameters - MOEAD*/
const float MOEAD_niche_multi = 0.1f;								///<multiplier of the neighbour size (* pop size). Default 0.1
const float MOEAD_limit_multi = 0.01f;								///<multiplier of the maximal number of solutions replaced(* pop size). Default 0.01
const float MOEAD_mating_chance = 0.9f;								///<probability chance for mating. Default 0.9
const int MOEAD_new_sol_multi = 5;									///<How many new individuals are generated each generation = pop_size/MOEAD_new_sol_multi. Default 5
const int MOEAD_evol_operation_type = 1;							///<Define type of evol. operations; 1 - DE cross. and poly. mut.; 2 - LL cross. and mut.. Default 1
const std::string MOEAD_weight_generator{ "Uniform" };				///<Define how the MOEA/D generate weight vectors; "File" - get from file in \Input\MOEADWeight; "Sobol" - generated by sobol - not implemented properly; "Uniform" - uniform method from MSF/PSF and mine. Default "Uniform"			
const std::string MOEAD_weight_assign{ "Pop" };						///<Define how the weight wectors are assigned to the population - in case of MLS; "Pop" -Assigned to population before the collective split; "Col" - Assigned to separately to each collective; "Sep" - Each collective have separate weight vector from 0. to 1. (i.e. each covers the whole search space). Default "Pop".

/**GA parameters - MTS*/
const int MTS_LocalSearch_test_amount = 1;							///<How many times each local search is tested each generation. Default 1
const int MTS_LocalSearch_amount = 1;								///<How many local searches are made each generation. Default 1
const float MTS_foreground_multp = 0.125;							///<Size of the foreground. Default 0.125
const int MTS_MPL1 = 9;												///<Bonus 1. Default 9
const int MTS_MPL2 = 2;												///<Bonus 2. Default 2
const int MTS_PF_refine_size = 1000;								///<Size of the PF after normal refining, when PF is too big (not the final one) for MTS only. Default 1000

/**GA parameters - BCE*/
const short BCE_mode = 2;											///<setting variation in the NPC and PC evolution;  1: SBX crossover, 2: DE. Default 2

/*GA parameters - HEIA*/
const float HEIA_crossover_prob = 0.5f;	///<Default 0.5f			///<Crossover probabilty for HEIA (between SBX and DE). Default 0.5

/**GA parameters - M2M*/
const short M2M_weight_method = 2;									///<Which weight method is used for weight creation (also used in other methods when nobj >= 3). 1 - reference point (not implemented); 2 - recusion point. Default 2
const short M2M_class_num = 10;										///<How many classes are defined. Default 10
const std::string M2M_method{ "MOEAD" }; 							///<Which evolutionary method is used from: "M2M", "UNSGAIII", "MOEAD", "MTS";. Default "MOEAD"
const float M2M_select_pro = 0.9f;									///<Default 0.9f
const float M2M_global_pro = 0.7f;									///<Default 0.7f

/**GA parameters - Transgenerational memory - not working*/
const int TGM_size = 2;												///<How many generations are considered (1 - maternal, 2 - grandmaternal). Default 2
const float TGM_m_min = -0.f;										///<Minimal maternal effect
const float TGM_m_max = 0.f;										///<Maximal maternal effect
const float TGM_m_step = 0.1f;										///<Change step of maternal effect
const float TGM_g_min = -0.f;										///<Minimal grandmaternal effect
const float TGM_g_max = 0.f;										///<Maximal grandmaternal effect
const float TGM_g_step = 0.1f;										///<Change step of grandmaternal effect

/**Precision*/
const unsigned short prec = 6;										///<Precision of results showing. Default 6
const unsigned short PF_prec = 16;									///<Precision of PF. Default 16
const int PF_real_size = 10000; 									///<Size of the real PF - for the IGD calculation. Default 10000
const int PF_size = 1000;  											///<Output size of the calculated PF. Default 1000
const int PF_refine_size = 1000;									///<Size of the PF after normal refining, when PF is too big (not the final one). Default 1000
const unsigned short PF_res = 6; 									///<Pareto front search resolution. Default 6
const int c_plot_res = 500; 										///<Resolution of contour plot generation (higher -> better, but slower), Define size of the net. Default 500
const unsigned short c_plot_deg = 10;								///<Degree of contour plot - increase difference between best and worst points on contour plot, but high values may cover some of results. Default 10
const unsigned short min_col_size_multi = 30;						///<Min size multiplier for collective size; col_min_size = pop_size / min_pop_size_multi. Default 30
const unsigned short new_col_size_multi = 15;						///<Min size multiplier for collective size; col_min_size = pop_size / min_pop_size_multi. Default 15
const unsigned short pow_eq_zero = 14;								///<which decimal will decide if two numbers are equal. Default 14
const float cons_check_param = -0.000001;							///<When constrain is violated for IGD check. Default -0.000001


/**Random*/
const int rnd_seed = 1;												///<Seed for the pseudo random generator. Default 1

/**Output*/
const std::string real_PF_out = "PF_Real";							///<File name for real Pareto Front
const int temp_data_size = 100;										///<Define amount of runs after which the temp storage will be cleaned. Default 100


/**Input*/
const std::string input_name = "Input/LHCGA.csv";					///<File name and path for distribution input

/**Video*/
const unsigned short frame_h = 720;									///<Frame height for video output. Default 720
const unsigned short frame_w = 1280; 								///<Frame width for video output. Default 1280
const unsigned short FPS = 10;										///<Frame per second. Default 10
const unsigned short frame_img_multi = 1;							///<How many times one frame is added - increase if want to make video longer and find a suitable frame. Default 1
const unsigned short frame_skip = 1;								///<How many frames are skipepd - for PAES and MOEAD. Default 1

/**Parameters for DAS-CMOP*/
const float CMOP_eta = 0.f;
const float CMOP_zeta = 0.f;
const float CMOP_gamma = 0.5f;

#endif ///< !CONST_H
