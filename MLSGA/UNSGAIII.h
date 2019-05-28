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
#ifndef UNSGAIII_H
#define UNSGAIII_H
#include "Class.h"

namespace UNSGAIII {
	/*
	Calculate individuals using NSGAII algorithm
	@param Col - address of a given collective
	@param iGen - current generation
	*/
	void UNSGAIII_Calc(collective & col);
	void Init(short n_obj, short n_col, population & pop);
}

#endif // !UNSGAII_H
