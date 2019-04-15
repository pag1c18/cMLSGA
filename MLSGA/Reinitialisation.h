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
