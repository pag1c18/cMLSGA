#pragma once
//****************************************
//		 CONST header
//Defining algorithm parameters and constants
//	behaviour is set in Define.h
//*****************************************
#ifndef CONST_H
#define CONST_H

const double long pi = 3.141592653589793238462643383279502884197l;
const double long e = 2.71828182845904523536028747135266249775724709369995l;

#include <string>
#include <vector>

/*General parameters*/
/**GA_mode*
Normal - Standard MLSGA
NSGAII - MLSGA augmented with NSGAII
MOEAD - MLSGA augmented with MOEAD
PAES - MLSGA augmented with PAES
MTS - MLSGA augmented with MTS
**List of Functions**
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
*Dynamic Unconstrained*
66-73 UDF 1-9
74-83 JY 1-10
84-88 FDA 1-5
*Dynamic Constrained*
89-98 CDF1-10

*Applications*
99 - Genetics

*Not Implemented*
5 - ZDT5
JY4
JY9
JY10
FDA4
FDA5
*/
/*MLSGA setting*/
const std::vector<std::vector<std::string>> GA_mode{ {"HEIA"} };//{ "MOEADMSF" , "MOEADM2M" },{ "MOEADMSF" , "MOEAD" },{ "MOEADMSF" , "NSGAII" },{ "MOEADMSF" , "Normal" }, { "MOEADPSF" , "MOEADM2M" },{ "MOEADPSF" , "MOEAD" },{ "MOEADPSF" , "NSGAII" },{ "MOEADPSF" , "Normal" },{ "MOEADPSF" , "HEIA" },{ "MOEADPSF" , "MOEADMSF" } };// , { "Normal" }, { "NSGAII" }, { "MOEAD" }, { "MOEADPSF" }, { "MOEADMSF" }, { "MOEADM2M" }, { "MTS" }, { "BCE" }, { "HEIA" } };// , { "MOEAD", "Normal" }, { "MOEAD", "NSGAII" }};//{ "MTS", "Normal" },{ "MTS", "NSGAII" },{ "MTS", "MOEAD" },{ "MTS", "MOEADMSF" },{ "MTS", "MOEADPSF" },{ "MTS", "MOEADM2M" }, { "MTS", "BCE" },{ "MTS", "BCE" } };//, { "MOEAD" , "BCE" }, { "MOEAD" , "HEIA" }, { "BCE", "HEIA" } }; //{ "Normal", "HEIA" },{ "Normal", "BCE" },{ "Normal", "MOEADM2M" },{ "Normal", "NSGAII" },{ "Normal", "MOEAD" },{ "Normal", "MTS" },{ "Normal", "MOEADMSF" },{ "Normal", "MOEADPSF" }, { "NSGAII", "MOEAD" },{ "NSGAII", "MOEADM2M" },{ "NSGAII", "MTS" },{ "NSGAII", "BCE" },{ "MOEAD", "MOEADM2M" },{ "MOEAD", "MTS" },{ "MOEAD", "BCE" },{ "MOEADM2M" , "MTS" },{ "MOEADM2M" , "BCE" },{ "MTS" , "BCE" },{ "HEIA", "MOEAD" },{ "HEIA", "MOEADM2M" },{ "HEIA", "MTS" },{ "HEIA", "BCE" }, {"NSGAII", "HEIA"} };// { "Normal", "HEIA" }, { "Normal", "BCE" }, { "Normal", "MOEADM2M" }, { "Normal", "NSGAII" }, { "Normal", "MOEAD" }, { "Normal", "MTS" }, { "Normal", "MOEADMSF" }, { "Normal", "MOEADPSF" } };// , { "NSGAII", "MOEAD" }, { "NSGAII", "MOEADM2M" }, { "NSGAII", "MTS" }, { "NSGAII", "BCE" }, { "MOEAD", "MOEADM2M" }, { "MOEAD", "MTS" }, { "MOEAD", "BCE" }, { "MOEADM2M" , "MTS" }, { "MOEADM2M" , "BCE" }, { "MTS" , "BCE" } }; //list of modes that have to be used: {"Normal"}, {"NSGAII"}, {"MOEAD"}, {"MOEADPSF"}, {"MOEADMSF"}, {"MOEADM2M"}, {"MTS"}, {"BCE"}, {"HEIA"}, {"IBEA"}
const bool MLSGA_Hybrid = true;									//defines if the hybid algorithm is used or the original one
const bool MLSGA_norm_obj = false;		//Default false			//Define if the objectives for fitness evaluations should be normalised (for MSF and PSF, but true makes PSF inoperative)
const unsigned short MLSGA_n_col_b = 8;	//Default 8					//Numbers of collectives - begin ***with a choice of 4, 6, 8 (is overriden to 1 if MLSGA_Htbrid == false)
const unsigned short MLSGA_n_col_e = 8;	//Default 8					//Numbers of collectives - end - step is always 2 (is overriden to 1 if MLSGA_Htbrid == false)
const unsigned short MLSGA_n_col_elim = 1;	//Default 1				//Number of collectives elimiated ***with a choice of 1, 2  (is overriden to 0 if MLSGA_Htbrid == false)
const unsigned short MLSGA_col_elim_limit = 10; //Default 4			//Define how many genrations have to pass before col elimination occurs - if 1 elimiantion will occur every generation - MTS have override that it change to 1 - line 494 in MLSGA.cpp
const float MLSGA_elite_size = 0.1f; //Default 0.1					//How many individuals become elite <0.0-1.0>

/*Parameter setting*/
const int Pop_size_b = 600;	//Default 1000	-300				//Population size - beginning
const int Pop_size_e = 600;	//Default 1000						//Population size - end
const int Pop_size_step = 50;	//Default 200					//Step for population size change
const unsigned short n_func_b = 99;								//Index of the first function used
const unsigned short n_func_e = 99; 							//Index of the last function used
const unsigned short n_func_obj = 2;							//Number of objectives for scalable problems - DTLZ, WFG
const std::vector<std::string> func_skip{"Three_obj", "Many_obj" };	//Decides if specific functions will be skipped "Two_obj", "Three_obj", "Many_obj" - DTLZ, WFG
const unsigned short n_MLS_b = 7; //Default 7					//Index of the first MLS type used
const unsigned short n_MLS_e = 7; //Default 7					//Index of the last MLS type used
const unsigned short n_runs = 15; //Default 30 CEC'09			//Number of runs
const unsigned short n_mut = 1; //Default 1						//Number of mutation types
const unsigned short n_cross = 1; //Default 1					//Number of crossover types
const unsigned short n_select = 1; //Default 1					//Number of selection types
const unsigned short n_mode = 1; //Default 1; 8 - PAES			//Number of mode types - for PAES only
const float mut_prob_min = 0.08f; //Default 0.08				//Minimal mutation probability - for real coded it is the absolute value, for binary it is divided by number of variables
const float mut_prob_max = 0.08f; //Default 0.08				//Maximal mutation probability
const float mut_prob_step = 1.f; //Default 0.1					//Mutation probability step
const float cross_prob_min = 0.7f;	//Default 0.70				//Minimal crossover probability
const float cross_prob_max = 0.7f; //Default 0.70				//Maximal crossover probability
const float cross_prob_step = 0.10f; //Default 0.10				//Crossover probability step


/*Termination conditions - only one is used depending on T_con*/
const std::string T_con{ "nfes" };								//Which termination condition is used: "nfes" - iterations number (int max_iterations) line 55; "ngen" - number of generations (int max_generations) line 56
const int max_iterations = 300000;	//Default 300000			//max iterations
const int max_generations = 300;	//Default 300				//max genenertaions

/*Performance evaluation*/
const std::vector<std::string> perf_eval_type{ "HV"};				//Which performance evaluation methods are used "IGD" - inverted generational distance, "HV" - hypervolume


/*Dynamic functions parameter*/
const int n_steps_dyn = 10;	//Default 5							//Number of distinct steps; represents the severity of change
const int T_dyn = 10;	//Default 5								//Window where the dynamic problem remains cosntant
const std::vector<std::string> Reinit_Mode{"VP" }; //list of reinitialisation modes that have to be used "None", "BR", "CER", "VP"
const int gen_memory_size = 2;	//Default 2						//How many generations are saved and used for reinitialisation
const int JY_const_window = 0;	//Default 100					//Window where JY functions remain constatant


/*Additional GA parameters*/
const unsigned short Di_c = 20;	//Default 20					//Distribution index of crossover
const unsigned short Di_m = 20;	//Default 20					//Distribution index of mutation - suttipong mut
const double Di_m2 = 0.00001;									//Distribution index of mutation - PG mut	(1. - random mutation, 0. - no mutation)
static short fit_index_sel = 2;  //Default 2			//Which fitness is used for indiv selection
static short fit_index_col_sel = 1;  //Default 1		//Which fitness is used for col selection 
const int Binary_string_size = 1;	//Default					//Lenght of each variable in binary string
const int MPCross_size = 2;		//Default 2						//Defines how many separation points are in multi-point crossover (1 is two point crossver)

/*GA parameters - MOEAD*/
const float MOEAD_niche_multi = 0.1f;	//Default 0.1			//multiplier of the neighbour size (* pop size)
const float MOEAD_limit_multi = 0.01f;	//Default 0.01			//multiplier of the maximal number of solutions replaced(* pop size)
const float MOEAD_mating_chance = 0.9f;	//Default 0.9			//probability chance for mating
const int MOEAD_new_sol_multi = 5;		//Default 5				//How many new individuals are generated each generation = pop_size/MOEAD_new_sol_multi
const int MOEAD_evol_operation_type = 1;	//Default 1			//Define type of evol. operations; 1 - DE cross. and poly. mut.; 2 - LL cross. and mut.
const std::string MOEAD_weight_generator{ "Uniform" };	//Default "File"	//Define how the MOEA/D generate weight vectors; "File" - get from file in \Input\MOEADWeight; "Sobol" - generated by sobol - not implemented properly; "Uniform" - uniform method from MSF/PSF and mine			
const std::string MOEAD_weight_assign{ "Pop" };	//Default "Pop"			//Define how the weight wectors are assigned to the population - in case of MLS; "Pop" -Assigned to population before the collective split; "Col" - Assigned to separately to each collective; "Sep" - Each collective have separate weight vector from 0. to 1. (i.e. each covers the whole search space)

/*GA parameters - DMOEADD*/
const int DMOEADD_Dimension = 2;	//Default 2					//Number of subdimensions
const float DMOEADD_Offspring = 0.05f;	//Default 10			//Number of created offsprings
const int DMOEADD_Temp = 10000;		//Default 10000				//Temperature for DMOEADD

/*GA parameters - MTS*/
const int MTS_LocalSearch_test_amount = 1;	//Default 5				//How many times each local search is tested each generation
const int MTS_LocalSearch_amount = 1;	//Default 45				//How many local searches are made each generation
const float MTS_foreground_multp = 0.125;		//Default 0.125		//Size of the foreground
const int MTS_MPL1 = 9;	//Default 9									//Bonus 1
const int MTS_MPL2 = 2;	//Default 2									//Bonus 2
const int MTS_PF_refine_size = 1000;	//Default 100000		//Size of the PF after normal refining, when PF is too big (not the final one) for MTS only

/*GA parameters - BCE*/
//const int BCE_PC_capacity = 600;				// maximum size of the PC population - not used, now it is the size of the current population
const short BCE_mode = 2;						// setting variation in the NPC and PC evolution;  1: SBX crossover, 2: DE

/*GA parameters - HEIA*/
const float HEIA_crossover_prob = 0.5f;	//Default 0.5f		//Crossover probabilty for HEIA (between SBX and DE)

/*GA parameters - M2M*/
const short M2M_weight_method = 2;//Default 2				//Which weight method is used for weight creation (also used in other methods when nobj >= 3). 1 - reference point (not implemented); 2 - recusion point
const short M2M_class_num = 10;	//Default 10			//How many classes are defined
const std::string M2M_method{ "M2M" };					//Which evolutionary method is used from: "M2M", "NSGAII", "MOEAD", "MTS";
const float M2M_select_pro = 0.9f; //Default 0.9f
const float M2M_global_pro = 0.7f; //Default 0.7f

/*GA parameters - Transgenerational memory*/
const int TGM_size = 2;	//Default 2								//How many generations are considered (1 - maternal, 2 - grandmaternal)
const float TGM_m_min = -0.f;										//Minimal maternal effect
const float TGM_m_max = 0.f;										//Maximal maternal effect
const float TGM_m_step = 0.1f;									//Change step of maternal effect
const float TGM_g_min = -0.f;										//Minimal grandmaternal effect
const float TGM_g_max = 0.f;										//Maximal grandmaternal effect
const float TGM_g_step = 0.1f;									//Change step of grandmaternal effect

/*Precision*/
const unsigned short prec = 6;	//Default 6						//Precision of results showing
const unsigned short PF_prec = 16;	//Default 16				//Precision of PF
const int PF_real_size = 10000; //Default 10000					//Size of the real PF - for the IGD calculation
const int PF_size = 1000;  //Default 100 - CEC'09				//Output size of the calculated PF
const int PF_refine_size = 1000;	//Default 700					//Size of the PF after normal refining, when PF is too big (not the final one) - should be at least 2x PF_size - not for MTS
const unsigned short PF_res = 6; //Default 5					//Pareto front search resolution
const int c_plot_res = 500; //Default 500						//Resolution of contour plot generation (higher -> better, but slower). Define size of the net
const unsigned short c_plot_deg = 10; //Default 10				//Degree of contour plot - increase difference between best and worst points on contour plot, but high values may cover some of results
const unsigned short min_col_size_multi = 30;//Default 30		//Min size multiplier for collective size; col_min_size = pop_size / min_pop_size_multi
const unsigned short new_col_size_multi = 15;//Default 15		//Min size multiplier for collective size; col_min_size = pop_size / min_pop_size_multi
const unsigned short pow_eq_zero = 14;	//Default 14			//which decimal will decide if two numbers are equal
const float cons_check_param = -0.000001;	//Default -0.000001 - CEC09	//When constrain is violated for IGD check


/*Random*/
const int rnd_seed = 1;	//Default 1								//Seed for the pseudo random generator

/*Output*/
const std::string real_PF_out = "PF_Real";						//File name for real Pareto Front
const int temp_data_size = 100;	//Default 500					//Define amount of runs after which the temp storage will be cleaned


/*Input*/
const std::string input_name = "Input/LHCGA.csv";				//File name and path for distribution input

/*Video*/
const unsigned short frame_h = 720; //Default 720				//Frame height for video output
const unsigned short frame_w = 1280;  //Default 1280			//Frame width for video output
const unsigned short FPS = 10;	//Default 3						//Frame per second
const unsigned short frame_img_multi = 1;	//Default 1			//How many times one frame is added - increase if want to make video longer and find a suitable frame
const unsigned short frame_skip = 5;							//How many frames are skipepd - for PAES and MOEAD
#endif // !CONST_H
