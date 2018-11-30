#pragma once
#ifndef SOBOL_H
#define SOBOL_H
#include <vector>

/*
*Function for getting initial points according to Sobol_Sequence*
@param size size of the popualtion
@param dimensions amount of variables
*/
std::vector<std::vector<double>> Sobol_Sequence(unsigned size, unsigned  dim);



#endif // !SOBOL_H
