#pragma once
#ifndef REINITIALISATION_H
#define REINITIALISATION_H
#include "Class.h"
#include <ostream>

/*
*Reinitialize the current population and Pareto Front storage according to chosen mode*
@param col - current population as a set of collectives
@param PF - current Pareto Front
@param re_mode - Reinitialisation mode
@pram mode - current hybrid mode
@param graph_v_vector - vector of graphs ofstream, for VIDEO creation
*/
void Reinit(std::vector<collective> & col, pareto_front & PF, std::string & re_mode, std::string & mode, std::vector<std::ofstream> & graph_v_vector, std::vector<std::vector<collective>> & col_memory);

#endif // !REINITIALISATION_H
