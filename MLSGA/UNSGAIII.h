/*CODE BASED ON THE ORIGINAL NSGAII CODE FROM THE WEBSITE*/

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
}

#endif // !UNSGAII_H
