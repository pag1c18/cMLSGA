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
//    MLSGA ADD FUNCTIONS header
//Storage of additional MLSGA functions and classes
//****************************************

#ifndef MLSGA_ADD_FUNCTIONS_H
#define MLSGA_ADD_FUNCTIONS_H

#include <fstream>
#include <string>

#include "Class.h"


/*
*Function selection*
@parameter ix index of function
*/
function * Func_Type(short ix);

/*
*Function deactivation - for specific functions*
*/
std::vector<bool> Func_Deactivate();

/*
Mutation type selection*
@parameter ix index of mutation type
@parameter Mode - Mode of MLSGA
*/
mutation<short> * Mut_Type(short ix, std::string Mode);
/*
Selection type selection*
@parameter - ix index of selection
*/
selection<individual> * Select_Type(short ix);
/*
Crossover type selection*
@parameter ix - index of crossover
*/
crossover<individual> * Cross_Type(short ix, std::string MODE);
std::string Date_Get();
/*
*Getting name of a file including date as a string*
@param name - name of the file
@param date - date string
@param ir - index number of the file
@param time - current time state
*/
std::string Name_Get(char * name, std::string & date, int ir, float time);
/*
*Getting name of a file including date as a string*
@param name - name of the file
@param date - date string
@param ir - index number of the file
*/
std::string Name_Get(char * name, std::string & date, int ir);
/*
*Getting name of a file as a string*
@param name - name of the file
@param n - index number of the file
*/
std::string Name_Get(std::string name, int n);
/*
*Getting name of a file as a string*
@param name - name of the file
@param n - index number of the file
@param t - time state
*/
std::string Name_Get(std::string name, int n, double t);
/*
*Getting name of a folder*
@param name_front front part of the folder path
@param date date string
@param name_back back part of the folder path
*/
std::string F_Name_Get(char * name_front, std::string & date, char * name_back);
/*
*GNUplot printing*
@param pipe pipe which will be used for printing
@param ir index of the run
@param index_r_max index of the max run
@param date date of the run
@param t time state of the run
@param IGD_on - if IGD was calculated and should be printed
@param HV_on - if HV was calculated and should be printed
*/
void Print(FILE * Pipe, int ir, int index_r_max, std::string date, double t, bool IGD_on, bool HV_on);

void Directory_Create(std::string &date);

/*
*Time saving*
@param ir - index of the run
@param IGD_val - IGD value of the current run
@param HV_val - HV value of the current run
@param min_fitness - min fitness for the current run
@param end_generation - end generation of the current run
*/
void Time_Save(int ir, double IGD_val, double HV_val, double min_fitness, int end_generation, int iteration, std::string MODE);
/*
*Excel output generation - First row*
@param - ir index of the run
@param dyn - if function is dynamic
@param IGD_on - if IGD was calculated and should be printed
@param HV_on - if HV was calculated and should be printed
*/
void Excel_F_Row(int ir, std::vector<std::string> MODE, bool dyn, bool IGD_on, bool HV_on);
/*
*Excel output generation - Second row*
@param fcode - function class type used for individuals creation
@param gapara - GA parameters class type used for parameters
@param mcode - mutation class type used for mutation
@param scode - selection class type used for selection
@param ccode - crossover class type used for crossover
@param irun - index of the run
@param indexr - index of the iteration
@param rimode - Reinitialisation mode used
@param mode - Mode used
@param MLS - MLS type used
*/
void Excel_S_Row(function & fcode, GA_parameters & gapara, mutation<short> & mcode, selection<individual> & scode, crossover<individual> & ccode, int irun, int ir, std::string rimode, std::vector<std::string> MODE, short MLS);
/*
*Excel output generation - GA data*
@param success_num - how many successful runs
@param data GA_data - class type with data
@param dyn - if function is dynamic
@param IGD_on - if IGD was calculated and should be printed
@param HV_on - if HV was calculated and should be printed
*/
void Excel_GA_data(int success_num, GA_data<float> & data, bool dyn, bool IGD_on, bool HV_on);
/**Getting current date and time as a string**/
#endif // !MLSGA_ADD_FUNCTIONS_H

/*Checking if everything is defined properly.Looking for errors in define.h and const.h*/
void Error_Check();

/*Show approx. iteration number and time*/
int Show_Iter();

/*
*Calculate the IGD statistics over the generations and runs and save to the file*
@param IGD_s - storage of the IGDs over generations and runs
@param index - index of the run
@param name - name of the current run
@param pipe - pipe which will be used for printing
*/
void IGD_Gen_Calc(std::vector<std::vector<double>> & IGD_s, int index, std::string & name, FILE * pipe);

/*
*Calculate the HV statistics over the generations and runs and save to the file*
@param HV_s - storage of the HVs over generations and runs
@param index - index of the run
@param name - name of the current run
@param pipe - pipe which will be used for printing
*/
void HV_Gen_Calc(std::vector<std::vector<double>> & HV_s, int index, std::string & name, FILE * pipe);

/*
*Create the fitness definitions for all collectives*
@param col - vector of collectives
@param MLS - MLS mode selected
*/
void MLS_Col_Fit_Create(std::vector<collective> & col, short MLS);






