/**MLSGA framework
Copyright(C) 2019  Przemyslaw A.Grudniewski and Adam J.Sobey

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see < https://www.gnu.org/licenses/>.

[1] A.J.Sobey and P.A.Grudniewski, “Re - inspiring the genetic algorithm with multi - level selection theory : Multi - level selection genetic algorithm”, Bioinspiration and Biomimetics, vol. 13, no. 5, pp. 1–13, 2018.

[2] P.A.Grudniewski and A.J.Sobey, “Behaviour of Multi - Level Selection Genetic Algorithm(MLSGA) using different individual - level selection mechanisms”, Swarm Evol.Comput., vol. 44, no.September 2018, pp. 852–862, 2018.

[3] P.A.Grudniewski and A.J.Sobey, "cMLSGA: co-evolutionary Multi-Level Selection Genetic Algorithm", 2019


If You have any questions regarding the code, please contact
Przemyslaw A. Grudniewski	at	pag1c18@soton.ac.uk		or
Adam J. Sobey				at	ajs502@soton.ac.uk


				MSLGA header




*/



#pragma once
#include <Windows.h>

#include "Reinitialisation.h"
#include "SVM.h"
#include "Clustering.h"
#include "Random_Label.h"
#include "Support_Functions.h"
#include "UNSGAIII.h"
#include "MOEAD.h"
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



//****************************************
//				Content - outdated
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
MLS 7 - in some col MLS2 in some MLS2R in some MLS1
*/


//****************************************
//			Imported content
//****************************************

//****************************************
//				Error list - outdated
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
