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

		IGD header
		Functions used for IGD calculation

*/

#pragma once
#ifndef IGD_H
#define IGD_H
#include "Class.h"
#include <vector>


/**IGD calculation
@param PF - Pareto Front for which IGD will be calculated. It is output from optimisation process.
@param real_PF - real Pareto Front for IGD calculation. If known, otherwise IGD cannot be calculated.
*/
double IGD_calc(pareto_front & PF, std::vector<std::vector<double>> & real_PF);

#endif // !IGD_H
