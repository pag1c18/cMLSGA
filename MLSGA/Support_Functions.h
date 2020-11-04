/**
*Copyright(C) 2019  Przemyslaw A.Grudniewski and Adam J.Sobey
*
*This file is part of the MLSGA framework
*
*The MLSGA framework is free software : you can redistribute it and/or modify
*it under the terms of the GNU General Public License as published by
*the Free Software Foundation, either version 3 of the License, or
*any later version.
*
*The MLSGA framework is distributed in the hope that it will be useful,
*but WITHOUT ANY WARRANTY; without even the implied warranty of
*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
*GNU General Public License for more details.
*
*You should have received a copy of the GNU General Public License
*along with this program.If not, see < https://www.gnu.org/licenses/>. 
*
*    SUPPORT FUNCTIONS header.
* Storage of additional math functions and classes
*/


#pragma once





#ifndef SUPPORT_FUNCTIONS_H
#define	SUPPORT_FUNCTIONS_H
#include <vector>
#include "Class.h"

/**
*Return the string of the given float with certain precision.
*@param v - value to copy to the string
*@param p - desired precision-> decimal spaces
*/
std::string String_Prec(float v, int p);

/**
*Calculation of the euclidean distance between 2 points.*/
template <typename tname>
double Distance(const std::vector<tname> & vec1, const std::vector<tname> & vec2);

/**
* Calculation of the distance between 2 points with 2 diffent values types.*/
template <typename tname,typename tname2>
double Distance2(const std::vector<tname> & vec1, const std::vector<tname2> & vec2);

/**
* Calculation of the value sign.
* Return -1 if val <0; 1 if val > 0; 0 if val == 0.
**/
template <typename tname>
short sgn(tname val);

/**
*Calculation of the euclidean distance between 2 points.*/
template double Distance<double>(const std::vector<double> & vec1, const std::vector<double> & vec2);
/**
* Calculation of the distance between 2 points with 2 diffent values types.*/
template double Distance2<double, float>(const std::vector<double> & vec1, const std::vector<float> & vec2);

template short sgn<double>(double val);

/**
*Check Dominance between two individuals.*
* @param ind1 - first individual
* @param ind2 - 2nd individual
* @param fit_indexes - vector of the fitness indexes for which dominance will be evaluated
* Return:
* 1 - ind1 dominates
* -1 - ind2 dominates
* 0 - both nondominated
*/
int Dominance_Check(individual &ind1, individual &ind2, std::vector<short> &fit_indexes);
/**
*Check Dominance  between two individuals, but without constraint check.
* @param ind1 - first individual
* @param ind2 - 2nd individual
* @param fit_indexes - vector of the fitness indexes for which dominance will be evaluated
* Return:
* 1 - ind1 dominates
* -1 - ind2 dominates
* 0 - both nondominated
*/
int Dominance_Check_NCons(individual &ind1, individual &ind2, std::vector<short> &fit_indexes);

/**
*Fitness calculation with utilisation of Scalarizing Function.
*@param fit - objective vector for which the fitness will be calculated according to given weights
*@param namda - weight vector
* @param fit_indexes - vector of the fitness indexes for which the calculation will be performed
* @param iGen - current generation (as some scalarizing function are adapting to generation number)
*/
double Fitness_Function(std::vector <double> &fit, std::vector <double> &namda, std::vector<short> & fit_indexes,  int iGen);
/**
*Fitness calculation by Chebycheff Scalarizing Function.
*@param fit - objective vector for which the fitness will be calculated according to given weights
*@param namda - weight vector
* @param fit_indexes - vector of the fitness indexes for which the calculation will be performed
*/
double Fitness_Function_TCH(std::vector <double> &fit, std::vector <double> &namda, std::vector<short> & fit_indexes);
/**
*Fitness calculation by PSF Scalarizing Function.
*@param fit - objective vector for which the fitness will be calculated according to given weights
*@param namda - weight vector
* @param fit_indexes - vector of the fitness indexes for which the calculation will be performed
* @param iGen - current generation (as some scalarizing function are adapting to generation number)
*/
double Fitness_Function_PSF(std::vector <double> &fit, std::vector <double> &namda, std::vector<short> & fit_indexes,  int iGen);
/**
*Fitness calculation by MSF Scalarizing Function.
*@param fit - objective vector for which the fitness will be calculated according to given weights
*@param namda - weight vector
* @param fit_indexes - vector of the fitness indexes for which the calculation will be performed
* @param iGen - current generation (as some scalarizing function are adapting to generation number)
*/
double Fitness_Function_MSF(std::vector <double> &fit, std::vector <double> &namda, std::vector<short> & fit_indexes,  int iGen);

/**
*Generation of the weight vector for whole population, subject to the problem type.
* @param pop_size - number of points generated (typically whole population)
* @param nobj - number of objectives (dimension of vector)
* @param type - type of problem: 1-concave, 2-convex, 3-linear
*/
std::vector<std::vector<double>> Uniform_Weights_Generate(int pop_size, short nobj, short type);

/**
*Generator of uniformly spread points.
*@param size - number of points generated
*@param n_obj - number of objectives (dimension)
*/
std::vector<std::vector<double>> Points_Generate(int size, short n_obj);


/**Normalise and return the objectives according to nadir and ideal point.
* @param fit - objective vector for which the normalisation will be performed*/
std::vector<double> Normalize_Objective(std::vector<double> &fit);

/**Clear and recreate the idealpoint and nadirpoint
* @param nobj - number of objectives (dimension)
*/
void Clear_Ideal_and_Nadir(short nobj);

/**Create nadirpoint for normalisation
* @param nobj - number of objectives (dimension)
*/
void Create_Nadir(short nobj);

/**Update nadirpoint for given population.
* @param pop - population according to which the nadirpoint will be updated.
* @param flag - indicator if the nadirpoint has to be recreated. If flag == 1, the nadir is zeroed as the new generation occurs.
*/
void Update_Nadirpoint(std::vector<individual>& pop, short flag);

/**
* Update the reference - idealpoint vector.*
@param ind - individual used for update.
*/
void Idealpoint_Update(std::vector<double> &ind);

double norm_vector(std::vector<double> &x); //from MSF/PSF

double innerproduct(std::vector<double> &vec1, std::vector<double> &vec2);

/**
* Evaluation of the termination condition.*
@param param - Define if the optimisation is dynamic. true if it is.
*/
bool Termination_Check(bool param = false);

#endif // !SUPPORT_FUNCTIONS_H



