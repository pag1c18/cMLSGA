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
along with this program.If not, see < https://www.gnu.org/licenses/>.

   MLSGA ADD FUNCTIONS header
Storage of additional functions,  specific for MLSGA framework. Most functions were implemented to test multiple problems in series.


*/


#pragma once



#ifndef MLSGA_ADD_FUNCTIONS_H
#define MLSGA_ADD_FUNCTIONS_H

#include <fstream>
#include <string>

#include "Class.h"


/**
*Function selection. Returns the pointer to selected function class
@parameter ix - index of function
*/
function * Func_Type(short ix);

/**
*Function deactivation. Deactivate specific functions according to const.h 
*/
std::vector<bool> Func_Deactivate();

/**
Mutation selection. Returns the pointer to selected mutation class
@parameter ix - index of mutation type
@parameter Mode - Mode of MLSGA. Specifically which GA is used.
*/
mutation<short> * Mut_Type(short ix, std::string Mode);
/**
Selection type selection. Returns the pointer to selected selection class
@parameter - ix index of selection
*/
selection<individual> * Select_Type(short ix);
/**
Crossover selection. Returns the pointer to selected crossover class
@parameter ix - index of crossover
@parameter Mode - Mode of MLSGA. Specifically which GA is used.
*/
crossover<individual> * Cross_Type(short ix, std::string MODE);

/**Returning the string containing current date and time*/
std::string Date_Get();
/**
*Providining path to the file including date as a string
@param name - name of the file
@param date - date string
@param ir - index number of the file (run)
@param time - current time state. For dynamic problems
*/
std::string Name_Get(char * name, std::string & date, int ir, float time);
/**
*Providining path to the file including date as a string. For static problems
@param name - name of the file 
@param date - date string
@param ir - index number of the file (run)
*/
std::string Name_Get(char * name, std::string & date, int ir);
/**
*Providining the name of the file (without folder path). For static problems
@param name - name of the file
@param n - index number of the file
*/
std::string Name_Get(std::string name, int n);
/**
*Providining the name of the file (without folder path). For dynamic problems
@param name - name of the file
@param n - index number of the file
@param t - current time state. For dynamic problems
*/
std::string Name_Get(std::string name, int n, double t);
/**
Getting name of a folder (alone).
@param name_front - front part of the folder path
@param date - date string
@param name_back - back part of the folder path
*/
std::string F_Name_Get(char * name_front, std::string & date, char * name_back);
/**
Routine for GNUplot printing. Provides graphs of Pareto Front and search areas as an output.
@param pipe - pipe which will be used for printing. Pointer to the open GNUplot cmd.
@param ir - index of the run
@param index_r_max - index of the last run
@param date - date of the run
@param t - time state of the run. For dynamic problems
@param IGD_on - if IGD was calculated and should be printed
@param HV_on - if HV was calculated and should be printed
*/
void Print(FILE * Pipe, int ir, int index_r_max, std::string date, double t, bool IGD_on, bool HV_on);


/**
*Create the folder for storing results.
@param date - date string
*/
void Directory_Create(std::string &date);

/**
Saving the run time of optimisation. As a whole and separate modules.
@param ir - index of the run
@param IGD_val - IGD value of the current run
@param HV_val - HV value of the current run
@param min_fitness - min fitness for the current run
@param end_generation - end generation of the current run
@param iteration - current iteration
@param MODE - current GA mode
*/
void Time_Save(int ir, double IGD_val, double HV_val, double min_fitness, int end_generation, int iteration, std::string MODE);
/** Generation of the Excel output, First row.*
@param ir - index of the run
@param MODE - current GA mode
@param dyn - if function is dynamic
@param IGD_on - if IGD was calculated and should be printed
@param HV_on - if HV was calculated and should be printed
*/
void Excel_F_Row(int ir, std::vector<std::string> MODE, bool dyn, bool IGD_on, bool HV_on);
/**
Excel output generation, Second row.*
@param fcode - function class type used for individuals creation
@param gapara - GA parameters class type used for parameters
@param mcode - mutation class type used for mutation
@param scode - selection class type used for selection
@param ccode - crossover class type used for crossover
@param irun - index of the current run
@param ir - index of the current iteration
@param rimode - Reinitialisation mode used
@param mode - Mode used
@param MLS - MLS type used
*/
void Excel_S_Row(function & fcode, GA_parameters & gapara, mutation<short> & mcode, selection<individual> & scode, crossover<individual> & ccode, int irun, int ir, std::string rimode, std::vector<std::string> MODE, short MLS);
/**Excel output generation, saving GA optimsiation data*
@param success_num - how many successful runs. In case of constrained search
@param data - GA_data class type with saved GA data
@param dyn - if function is dynamic
@param IGD_on - if IGD was calculated and should be printed
@param HV_on - if HV was calculated and should be printed
*/
void Excel_GA_data(int success_num, GA_data<float> & data, bool dyn, bool IGD_on, bool HV_on);


/**Checking if everything is defined properly. Looking for errors in define.h and const.h*/
void Error_Check();

/**Show approx. iteration number and time. Works only for simple problems*/
int Show_Iter();

/**
Calculate the IGD statistics over the generations and runs, and save them to the file/
@param IGD_s - storage of the IGDs over generations and runs
@param index - index of the run
@param name - name of the current run
@param pipe - pipe which will be used for printing
*/
void IGD_Gen_Calc(std::vector<std::vector<double>> & IGD_s, int index, std::string & name, FILE * pipe);

/**
Calculate the HV statistics over the generations and runs, and save them to the file/
@param HV_s - storage of the HVs over generations and runs
@param index - index of the run
@param name - name of the current run
@param pipe - pipe which will be used for printing
*/
void HV_Gen_Calc(std::vector<std::vector<double>> & HV_s, int index, std::string & name, FILE * pipe);

/**
Create the fitness definitions for all collectives
@param col - vector of collectives
@param MLS - MLS mode selected
*/
void MLS_Col_Fit_Create(std::vector<collective> & col, short MLS);


#endif // !MLSGA_ADD_FUNCTIONS_H



