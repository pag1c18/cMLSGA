//****************************************
//				MSLGA header
//****************************************

#pragma once
#include <Windows.h>

#include "Reinitialisation.h"
#include "SVM.h"
#include "Clustering.h"
#include "Random_Label.h"
#include "Support_Functions.h"
#include "NSGAII.h"
#include "MOEAD.h"
#include "PAES.h"
#include "DMOEADD.h"
#include "MTS.h"
#include "BCE.h"
#include "HEIA.h"
#include "M2M.h"
#include "IBEA.h"
#include "MLSGA_Add_Functions.h"
#include "Imported\Workbook.h"
#include "Video.h"
#include "IGD.h"
#include "Imported/hv_wfg.h"
//*****To DO
/*
Clean:
Class.h
Crossover.h
Fit_function.h
GA_functions.h
Mutation.h
Selection.h
*/



//****************************************
//				Notes
//****************************************
//Selection only for 1st objective fitness

//****************************************
//				Content
//****************************************
//Class_h - Contain main GA classes
//Const_H - Containts algorithm parameters
//Crossover_h - Contain crossover class
//Define_h - For defining benchmarking parameters and behaviour of algorithm
//Fit_Functions_h - Contain functions classes
//GA_Funtions_h - Storage of additional GA funtions and classes
//MLSGA_Add_Functions_h - Storage of additional MLSGA functions and classes
//Mutation_h - Contain mutation classes
//Selection_h - Contain selection classes
//Struct_h - Contain structures classes
//Support_Functions_h - Contains additional math functions
// /Imported - Imported tools (mentioned in imported content)
// /Tools - Tools


//Multi level selection types implemented
/*
MLS 1 - (f1 + f2)/2 on both levels
MLS 2 - f1 and f2 separated
MLS 3 - f1 and f2 separated among different collectives, collective var[0]
MLS 4 - col as MLS1 and indi as MLS2
MLS 5 - col as MLS2 and indi as MLS1
MLS 6 - in some col MLS2 in some MLS2R
MLS 7 - in some col MLS2 in some MLS2R in some MLS1
MLS 8 - in individual level (f1 + f2)/2 and on collective f1 or f2
MLS 9 - combination of MLS7 and MLS 8
*/


//****************************************
//			Imported content
//****************************************

//****************************************
//				Error list
//****************************************
//ERROR#01: INDIVIDUAL - CREATION
//ERROR#02: FUNCTION - CREATION
//ERROR#03: BOUNDARIES - SHOW
//ERROR#04: POPULATION - CREATION
//ERROR#05: COLLECTIVE - CREATION
//ERROR#06: PARETO FRONT- CREATION
//ERROR#07: MUTATION - CREATION
//ERROR#08: CROSSOVER - CREATION
//ERROR#09: SELECTION - CREATION
//ERROR#10: CROSSOVER - VECTOR SIZE
//ERROR#11: SELECTION - VECTOR SIZE
//ERROR#12: FUNCTION - VECTOR SIZE
//ERROR#13: POPULATION - VECTOR SIZE
//ERROR#14: SVM - LABEL
//ERROR#15: VIDEO - LABEL SIZE
//ERROR#16: ELITISM




//****************************************
//	  To Do - optimisation of code
//****************************************
//Copy the population wherever it is used (instead of calling the individuals)