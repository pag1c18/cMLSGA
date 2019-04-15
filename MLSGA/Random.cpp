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

#include "Random.h"
#include "Define.h"
#include "Const.h"
#include <random>

static std::random_device generator_real;				//random numbers generator device
static std::default_random_engine generator_pseudo(rnd_seed);

/*
*generate the random double value in the given range*
@param min min value (default 0.0)
@param max max value (default 1.0)
*/
double Random(double min, double max)
{
	double min_t = min;
	double max_t = max;
	if (min > max)
	{
		min_t = max;
		max_t = min;
	}
	//generate the distribution seed
	std::uniform_real_distribution<double> distribution(min_t, max_t);		//distribution seed
	double val;
	if (REAL_RANDOM == true)
	{
		//generate a random number
		val = distribution(generator_real);
	}
	else
	{
		//generate a random number
		val = distribution(generator_pseudo);
	}
	return val;
}

/*
*generate the random float value in the given range*
@param min min value (default 0.0f)
@param max max value (default 1.0f)
*/
float Random_F(float min, float max)
{
	float min_t = min;
	float max_t = max;
	if (min > max)
	{
		min_t = max;
		max_t = min;
	}
	//generate the distribution seed
	std::uniform_real_distribution<float> distribution(min_t, max_t);		//distribution seed
	float val;
	if (REAL_RANDOM == true)
	{
		//generate a random number
		val = distribution(generator_real);

		return val;
	}
	else
	{
		//generate a random number
		val = distribution(generator_pseudo);
	}
	return val;
}

/*
*generate the random integer value in the given range*
@param min min value
@param max max value
*/
int Random_I(int min, int max)
{
	int min_t = min;
	int max_t = max;
	if (min > max)
	{
		min_t = max;
		max_t = min;
	}
	//generate the distribution seed
	std::uniform_int_distribution<int> distribution(min_t, max_t);		//distribution seed
	int val;
	if (REAL_RANDOM == true)
	{
		//generate a random number
		val = distribution(generator_real);
	}
	else
	{
		//generate a random number
		val = distribution(generator_pseudo);
	}
	return val;
}
/*
*generate the value from the normal distribution*
@param mean mean for the normal distribution
@param std standard deviation
*/
double N_Distribution(double mean, double std)
{
	double val;
	//generate the distribution class
	std::normal_distribution<double> distribution(mean, std);

	//generate the number
	if (REAL_RANDOM == true)
	{
		val = distribution(generator_real);
	}
	else
	{
		val = distribution(generator_pseudo);
	}
	return val;
}
