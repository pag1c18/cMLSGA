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
#include "Shape.h"

void ConvertAngles(double *begin, double *end)
{
	double *i;
	for (i = begin; i != end; ++i)
	{
		assert(0 <= *i && *i <= 1);
		*i = *i * M_PI / 2;
	}
}

void ConvexShape(double *angleBegin, double *angleEnd, double *objectiveBegin, double *objectiveEnd, double radius)
{
	double *objective;
	assert(objectiveEnd - objectiveBegin > 0);
	{
		double temp = radius;
		double *angle;
		for (angle = angleBegin; angle != angleEnd; ++angle)
		{
			assert(0 <= *angle && *angle <= M_PI / 2);
			temp *= radius - cos(*angle);
		}
		*objectiveBegin = temp;
	}
	for (objective = objectiveBegin + 1; objective != objectiveEnd; ++objective)
	{
		size_t nObjective = objective - objectiveBegin;
		size_t index = (angleEnd - angleBegin) - nObjective;
		double *_angleEnd = angleBegin + index;
		double temp = radius;
		double *angle;
		for (angle = angleBegin; angle != _angleEnd; ++angle)
			temp *= radius - cos(*angle);
		*objective = temp * (radius - sin(*_angleEnd));
	}
}

double ConvexConcaveShape(double decision, size_t nSegments, double shape)
{
	double tmp;
	assert(0 <= decision && decision <= 1);
	assert(nSegments > 0);
	assert(shape > 0);
	tmp = 2 * nSegments * M_PI;
	return pow(1 - decision - cos(tmp * decision + M_PI / 2) / tmp, shape);
}

double DisconnectedShape(double decision, size_t nRegions, double shape, double location)
{
	double tmp;
	assert(nRegions > 0);
	assert(shape > 0);
	assert(location > 0);
	assert(0 <= decision && decision <= 1);
	tmp = nRegions * pow(decision, location) * M_PI;
	return 1 - pow(decision, shape) * pow(cos(tmp), 2);
}

void InvertConcaveShape(double *angleBegin, double *angleEnd, double *objectiveBegin, double *objectiveEnd, double radius)
{
	double *objective;
	assert(objectiveEnd - objectiveBegin > 0);
	{
		double temp = radius;
		double *angle;
		for (angle = angleBegin; angle != angleEnd; ++angle)
		{
			assert(0 <= *angle && *angle <= M_PI / 2);
			temp *= sin(*angle);
		}
		*objectiveBegin = temp;
	}
	for (objective = objectiveBegin + 1; objective != objectiveEnd; ++objective)
	{
		size_t nObjective = objective - objectiveBegin;
		size_t index = (angleEnd - angleBegin) - nObjective;
		double *_angleEnd = angleBegin + index;
		double temp = radius;
		double *angle;
		for (angle = angleBegin; angle != _angleEnd; ++angle)
			temp *= sin(*angle);
		*objective = temp * cos(*_angleEnd);
	}
}

void LinearShape(double *decisionBegin, double *decisionEnd, double *objectiveBegin, double *objectiveEnd, double distance)
{
	double *decision, *objective;
	assert(0 < decisionEnd - decisionBegin && decisionEnd - decisionBegin < objectiveEnd - objectiveBegin);
	*objectiveBegin = distance;
	for (decision = decisionBegin; decision != decisionEnd; ++decision)
		*objectiveBegin *= *decision;
	for (objective = objectiveBegin + 1; objective != objectiveEnd; ++objective)
	{
		size_t nObjective = objective - objectiveBegin;
		size_t index = decisionEnd - decisionBegin - nObjective;
		double *_decisionEnd = decisionBegin + index;
		*objective = distance;
		for (decision = decisionBegin; decision != _decisionEnd; ++decision)
			*objective *= *decision;
		*objective *= (1 - *_decisionEnd);
	}
}
