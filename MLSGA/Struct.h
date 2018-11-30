#pragma once
#ifndef STRUCT
#define STRUCT

namespace STRUCTURES
{
	
	struct boundaries				//Structure for boundaries 
	{ 
		double upper;				//Upper boundary
		double lower;				//Lower boundary
	};

	template <typename tname>
	struct min_max_avg_std
	{
		tname min;			//min value
		tname max;			//max value
		tname avg;			//average value
		tname std_deviation;//standard deviation
	};
}


#endif //!STRUCT