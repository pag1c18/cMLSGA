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
//    SUPPORT FUNCTIONS header
//Storage of additional math functions and classes
//****************************************


#ifndef SUPPORT_FUNCTIONS_H
#define	SUPPORT_FUNCTIONS_H
#include <vector>
#include "Class.h"

/*
*Give string with certain precision*
@param v - value to copy to the string
@param p - desired precision-> decimal spaces
*/
std::string String_Prec(float v, int p);

/**Calculation of the distance between 2 points**/
template <typename tname>
double Distance(const std::vector<tname> & vec1, const std::vector<tname> & vec2);

/**Calculation of the distance between 2 points - 2 diffent values types**/
template <typename tname,typename tname2>
double Distance2(const std::vector<tname> & vec1, const std::vector<tname2> & vec2);

/**Calculation of sign**/
template <typename tname>
short sgn(tname val);


template double Distance<double>(const std::vector<double> & vec1, const std::vector<double> & vec2);
template double Distance2<double, float>(const std::vector<double> & vec1, const std::vector<float> & vec2);

template short sgn<double>(double val);

/*
*Check Dominance*
@param ind1 - first individual
@param ind2 - 2nd individual
Return:
1 - ind1 dominates
-1 - ind2 dominates
0 - both nondominated
*/
int Dominance_Check(individual &ind1, individual &ind2, std::vector<short> &fit_indexes);
/*
*Check Dominance for MLSt*
@param ind1 - first individual
@param ind2 - 2nd individual
@param ix - inx of the collective
Return:
1 - ind1 dominates
-1 - ind2 dominates
0 - both nondominated
*/
int Dominance_Check2(individual &ind1, individual &ind2, int ix, std::vector<short> & fit_indexes);
/*
*Check Dominance for vector*
@param ind1 - first individual
@param ind2 - 2nd individual
Return:
1 - ind1 dominates
-1 - ind2 dominates
0 - both nondominated
*/
int Dominance_Check3(individual &ind1, individual &ind2);

/*
*Fitness calculation by Scalarizing Function*
@param fit - fitness vector
@param namda - weight vector
*/
double Fitness_Function(std::vector <double> &fit, std::vector <double> &namda, std::vector<short> & fit_indexes,  int iGen);
/*
*Fitness calculation by Chebycheff Scalarizing Function*
@param fit - fitness vector
@param namda - weight vector
*/
double Fitness_Function_TCH(std::vector <double> &fit, std::vector <double> &namda, std::vector<short> & fit_indexes);
/*
*Fitness calculation by PSF Scalarizing Function*
@param fit - fitness vector
@param namda - weight vector
*/
double Fitness_Function_PSF(std::vector <double> &fit, std::vector <double> &namda, std::vector<short> & fit_indexes,  int iGen);
/*
*Fitness calculation by MSF Scalarizing Function*
@param fit - fitness vector
@param namda - weight vector
*/
double Fitness_Function_MSF(std::vector <double> &fit, std::vector <double> &namda, std::vector<short> & fit_indexes,  int iGen);

/*
*Generation of the weight vector for whole population - subject to problem type* type 1-concave, 2-convex, 3-linear
@param pop_size - number of points generated
@param nobj - number of objectives (dimension)
@param type - type of weight vector 1-concave, 2-convex, 3-linear
*/
std::vector<std::vector<double>> Uniform_Weights_Generate(int pop_size, short nobj, short type);


/*
*Generation of uniformly spread points*
@param size - number of points generated
@param n_obj - number of objectives (dimension)
*/
std::vector<std::vector<double>> Points_Generate(int size, short n_obj);
/*
Routine to fast sort two asssigned vectors of equal size according to two control parameters
*/
void Min_Fast_Sort(std::vector<double> & x, std::vector<int> & idx, int n, int m);

//Normalise the objectives (according to:...)
std::vector<double> Normalize_Objective(std::vector<double> &fit);

//Clear idealpoint and nadirpoint for normalisation
void Clear_Ideal_and_Nadir(short nobj);

//Create idealpoint and nadirpoint for normalisation
void Create_Nadir(short nobj);

//update nadirpoint for given population (if flag == 1, the nadir is zeroed as the new generation occurs)
void Update_Nadirpoint(std::vector<individual>& pop, short flag);

/*
*Update the reference - idealpoint vector*
@param ind - individual used for update
*/
void Idealpoint_Update(std::vector<double> &ind);

double norm_vector(std::vector<double> &x); //from MSF/PSF

double innerproduct(std::vector<double> &vec1, std::vector<double> &vec2);

bool Termination_Check(bool param = false);

#endif // !SUPPORT_FUNCTIONS_H



