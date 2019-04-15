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
#ifndef RANDOM_H
#define RANDOM_H

/*
*generate the random double value in the given range*
@param min min value (default 0.0)
@param max max value (default 1.0)
*/
double Random(double min = 0.0, double max = 1.0);
/*
*generate the random float value in the given range*
@param min min value (default 0.0f)
@param max max value (default 1.0f)
*/
float Random_F(float min = 0.0f, float max = 1.0f);
/*
*generate the random integer value in the given range*
@param min min value 
@param max max value 
*/
int Random_I(int min, int max);

/*
*generate the value from the normal distribution*
@param mean mean for the normal distribution
@param std standard deviation
*/
double N_Distribution(double mean, double std);

#endif // !RANDOM_H