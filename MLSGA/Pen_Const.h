#pragma once
/*Copyright(C) 2019  Przemyslaw A.Grudniewski and Adam J.Sobey

This file is part of the MLSGA framework

The MLSGA framework is free software : you can redistribute it and /or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

The MLSGA framework is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see < https://www.gnu.org/licenses/>. */

#ifndef PEN_CONST_H
#define PEN_CONST_H
#include "Class.h"

namespace Pen_const
{
	//Add the penalty to the infeasible solutions (based on Azzouz, R., Bechikh, S., Said, L. Ben, & Trabelsi, W. (2018). Handling time-varying constraints and objectives in dynamic evolutionary multi-objective optimization. Swarm and Evolutionary Computation, 39(April 2017), 222–248)
	void Fitness_Recalc(std::vector<individual>& ind ,short ix, bool u);
	//Penalty for a single individual
	void Fitness_Recalc(individual& ind, short ix);


	//Add the penalty to the infeasible solutions (based on Azzouz, R., Bechikh, S., Said, L. Ben, & Trabelsi, W. (2018). Handling time-varying constraints and objectives in dynamic evolutionary multi-objective optimization. Swarm and Evolutionary Computation, 39(April 2017), 222–248)
	void Fitness_Recalc(collective & col, bool u);

	//initialise the global values
	void Initialise(short ncons, short n_col, short nobj);

	void Initialise_Col(int psize, short ix);

	void R_f_update(bool val, short ix);
}





#endif // !PEN_CONST_H

