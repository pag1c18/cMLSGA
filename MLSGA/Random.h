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