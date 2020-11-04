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

*    STRUCT header.
* Storage of additional struct types
*/



#pragma once
#ifndef STRUCT
#define STRUCT

namespace STRUCTURES
{

	/**
	* Storage to store boundaries for variables
	*/
	struct boundaries
	{
		double upper;				///<Upper boundary
		double lower;				///<Lower boundary
	};

	/**
	* Storage to store statistical data
	*/
	template <typename tname>
	struct min_max_avg_std
	{
		tname min;			///<min value
		tname max;			///<max value
		tname avg;			///<average value
		tname std_deviation;///<standard deviation
	};
};


#endif //!STRUCT