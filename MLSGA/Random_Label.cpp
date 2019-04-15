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


#include "Random_Label.h"
#include "Random.h"
/*
*Generate random labels for the population*
@param pop source population for which label will be found
*/

extern short n_col;

std::vector<short> Random_Label(const population & pop)
{
	int psize = pop.Size_Show();		//size of the population
	std::vector<short> label (psize,0);			//vector for labels storage - output

	//assign labels by random
	for (int i = 0; i < psize; i++)
	{
		label[i] = Random_I(1, n_col);
	}
	
	//return the vector
	return label;
}