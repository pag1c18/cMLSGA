#pragma once
#ifndef IGD_H
#define IGD_H
#include "Class.h"
#include <vector>


/*
*IGD calculation*
@param PF Pareto Front for which IGD will be calculated
@param real_PF real Pareto Front for IGD calculation
*/
double IGD_calc(pareto_front & PF, std::vector<std::vector<double>> & real_PF);

#endif // !IGD_H
