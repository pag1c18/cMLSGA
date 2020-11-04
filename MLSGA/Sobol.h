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

SOBOL header
Functions for population initialisation according to sobol sequence

*/

#pragma once
#ifndef SOBOL_H
#define SOBOL_H
#include <vector>

/**
*Function for getting initial points according to the Sobol_Sequence.
@param size - size of the popualtion (and thus output)
@param dimensions -  amount of variables.
*/
std::vector<std::vector<double>> Sobol_Sequence(unsigned size, unsigned  dim);



#endif // !SOBOL_H
