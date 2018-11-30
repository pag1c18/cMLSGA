/*!
Copyright (C) 2014, 申瑞珉 (Ruimin Shen)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#define _USE_MATH_DEFINES

#include <assert.h>
#include <math.h>
#include "Transform.h"

#define MIN( a, b ) ( ( a < b) ? a : b )
#define MAX( a, b ) ( ( a > b) ? a : b )

double Fix(double decision, double min, double max)
{
	assert(min < max);
	if (decision < min)
		return min;
	else if (decision > max)
		return max;
	else
		return decision;
}

void DegenerateTransform(double *decisionBegin, double *decisionEnd, double *degenerateBegin, double *degenerateEnd, double distance)
{
	double *degenerate;
	double *decision;
	assert(degenerateEnd - degenerateBegin > 0);
	assert(degenerateEnd - degenerateBegin <= decisionEnd - decisionBegin);
	for(degenerate = degenerateBegin, decision = decisionBegin; degenerate != degenerateEnd; ++degenerate, ++decision)
	{
		double tmp = MAX(distance, *degenerate);
		*decision = tmp * (*decision - 0.5) + 0.5;
	}
}

double FlatTransform(double decision, double region, double regionMin, double regionMax)
{
	double tmp1, tmp2;
	assert(0 <= decision && decision <= 1);
	assert(0 <= region && region <= 1);
	assert(0 <= regionMin && regionMin <= 1);
	assert(0 <= regionMax && regionMax <= 1);
	assert(regionMin < regionMax);
	assert(regionMin != 0 || region == 0);
	assert(regionMin != 0 || regionMax != 1);
	assert(regionMax != 1 || region == 1);
	assert(regionMax != 1 || regionMin != 0);
	tmp1 = MIN((double)0, floor(decision - regionMin)) * region * (regionMin - decision) / regionMin;
	tmp2 = MIN((double)0, floor(regionMax - decision)) * (1 - region) * (decision - regionMax) / (1 - regionMax);
	return Fix(region + tmp1 - tmp2, 0, 1);
}

double DependentTransform(double decision, double regionMapping, double regionFraction, double regionMin, double regionMax)
{
	double tmp;
	assert(0 <= decision && decision <= 1);
	assert(0 <= regionMapping && regionMapping <= 1);
	assert(0 < regionFraction && regionFraction < 1);
	assert(regionMin > 0);
	assert(regionMin < regionMax);
	tmp = regionFraction - (1 - 2 * regionMapping) * fabs(floor(0.5 - regionMapping) + regionFraction);
	return Fix(pow(decision, regionMin + (regionMax - regionMin) * tmp), 0, 1);
}

double PolynomialTransform(double decision, double bias)
{
	assert(0 <= decision && decision <= 1);
	assert(bias > 0 && bias != 1);
	return Fix(pow(decision, bias), 0, 1);
}

double NonSeparableTransform(double *begin, double *end, size_t degree)
{
	double numerator = 0, tmp, denominator;
	double *decision;
	assert(end - begin > 0);
	assert(0 < degree && degree <= end - begin);
	assert((end - begin) % degree == 0);
	for(decision = begin; decision != end; ++decision)
	{
		size_t i;
		const size_t nDecision = decision - begin;
		assert(0 <= *decision && *decision <= 1);
		numerator += *decision;
		for(i = 0; i <= degree - 2; ++i)
		{
			const size_t index = (nDecision + i + 1) % (end - begin);
			double *_decision = begin + index;
			numerator += fabs(*decision - *_decision);
		}
	}
	tmp = ceil((double)degree / 2);
	denominator = (end - begin) * tmp * (1 + 2 * degree - 2 * tmp) / degree;
	return Fix(numerator / denominator, 0, 1);
}

double WeightedSumTransform(double *weightBegin, double *weightEnd, double *decisionBegin)
{
	double innerProduct = 0, accumulate = 0;
	double *weight, *decision;

	assert(weightEnd - weightBegin > 0);
	for (weight = weightBegin, decision = decisionBegin; weight != weightEnd; ++weight, ++decision)
		innerProduct += *weight * *decision;
	for (weight = weightBegin; weight != weightEnd; ++weight)
		accumulate += *weight;
	return Fix(innerProduct / accumulate, 0, 1);
}

double LinearTransform(double decision, double globalMin)
{
	assert(0 <= decision && decision <= 1);
	assert(0 < globalMin && globalMin < 1);
	return Fix(fabs(decision - globalMin) / fabs(floor(globalMin - decision) + globalMin), 0, 1);
}

double DeceptiveTransform(double decision, double globalMin, double aperture, double deceptiveMin)
{
	double tmp1, tmp2;

	assert(0 <= decision && decision <= 1);
	assert(0 < globalMin && globalMin < 1);
	assert(0 < aperture && aperture < 1);
	assert(0 < deceptiveMin && deceptiveMin < 1);
	assert(globalMin - aperture > 0);
	assert(globalMin + aperture < 1);
	tmp1 = floor(decision - globalMin + aperture) * (1 - deceptiveMin + (globalMin - aperture) / aperture) / (globalMin - aperture);
	tmp2 = floor(globalMin + aperture - decision) * (1 - deceptiveMin + (1 - globalMin - aperture) / aperture) / (1 - globalMin - aperture);
	return Fix(1 + (fabs(decision - globalMin) - aperture) * (tmp1 + tmp2 + 1 / aperture), 0, 1);
}

double MultiModalTransform(double decision, const size_t nMinima, double hillSize, double globalMin)
{
	double tmp1, tmp2;

	assert(0 <= decision && decision <= 1);
	assert(nMinima >= 1);
	assert(hillSize >= 0);
	assert((4 * nMinima + 2) * M_PI >= 4 * hillSize);
	assert(0 <= globalMin && globalMin <= 1);
	tmp1 = fabs(decision - globalMin) / (2 * (floor(globalMin - decision) + globalMin));
	tmp2 = (4 * nMinima + 2) * M_PI * (0.5 - tmp1);
	return Fix((1 + cos(tmp2) + 4 * hillSize * pow(tmp1, 2)) / (hillSize + 2), 0, 1);
}
