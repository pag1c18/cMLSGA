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
along with this program.If not, see < https://www.gnu.org/licenses/>. */

#pragma once
#ifndef RANDOM_LABEL_H
#define RANDOM_LABEL_H

#include <vector>
#include "Class.h"

/**
*Generate random labels for the population*
@param pop source population for which label will be found
*/
std::vector<short> Random_Label(const population & pop);

#endif // !RANDOM_LABEL_H
