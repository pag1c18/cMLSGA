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


#pragma once

//****************************************
//		 DEFINE header
//	Defining algorithm behaviour
//	paramteres are set in Const.h
//****************************************


#ifndef DEFINE_H
#define DEFINE_H
#define _CRT_SECURE_NO_WARNINGS
#define INF 1.0e14

#define VIDEO false								//Video generation
#define PERF_VAL_GEN false						//If IGD or HV value should be calculated every generation and stored to file
#define AUTO_FOLDER_COMMENT true				//If the name of folder should be created automatically and containt the run parameters
#define ENCODING "Real"							//Encoding type "Binary", "Gray", "Real"
#define PENALTY_BASED_CONSTRAINTS false			//If the penalty functions have to calculated instead of constrained ranking

#define PF_REFINE true							//If PF is refining at the end to the desired size (it will be refined anyway, if PF is too big - see const.h)
#define DEBUG false								//For debugging purposes
#define FITNESS_ALL	true						//save all fitness to the excel and show if not save only PF
#define SKIP_GRAPHS	true						//skip making graphs for runs with index > 5; for very long dynamic runs to save memory
#define EXCEL_EXCEPTION true					//do not save fitness to excel - useful for very long runs
#define ONE_OBJ_OVERRIDE false					//Calculating as a single objective (every objective have the same weight)										
#define CONTOUR_PLOT true						//Contour plot generation (FITNESS_ALL have to be defined!!)
#define REAL_RANDOM false						//if true random will be real random, if false random will be pseudo-random
#define GENERATION_GRAPHS_DATA false			//if graphs after 50,100 and 250 generations are wanted
#define FILE_INPUT false						//if values are read from file
#define SOBOL false								//if values are generated from sobol
#define MOEAD_OVERRIDE false					//if MOEAD should use crossover.cpp instead of standard DE crossover
//#define MOEAD_SOBOL false						//if namda in MOEAD are calculated from sobol
//#define MOEAD_SEPARATE_WEIGHTS false			//define if each col will have different weight vectors
#define GROUPING "Random"
#define MLS_OVERRIDE true						//if individual levle selection should calculate like NSGAII or MLS; if true NSGAII will not use MLS for selection
#define SKIP_LAST_GEN false						//if should skip the last generations (over 200k)
#define TGM false								//If transgenerational memory would oocur

#endif // !DEFINE_H
