#pragma once
#ifndef DMOEADD_H
#define DMOEADD_H
#include "Class.h"

namespace DMOEADD {
	std::vector<individual> DMOEADD_Calc(collective & col);
}

#endif // !DMOEADD_H



/*Questions
- if they calculate crowding distance among all individuals or only in current rank
- how they calculate fitness for crowding INF
- if for sure rank is not as in NSGAII
- What selection, and what are priciples of it
- What mutation, and how to select new indi
- 
*/
