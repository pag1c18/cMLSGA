#pragma once
#ifndef RANDOM_LABEL_H
#define RANDOM_LABEL_H

#include <vector>
#include "Class.h"

/*
*Generate random labels for the population*
@param pop source population for which label will be found
*/
std::vector<short> Random_Label(const population & pop);

#endif // !RANDOM_LABEL_H
