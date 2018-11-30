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

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include "WFG.h"
#include "Transform.h"
#include "Shape.h"





void WFG1(double *decisionBegin, double *decisionEnd, double *objectiveBegin, double *objectiveEnd, const size_t posGroups, const size_t distDecisions);
void WFG2(double *decisionBegin, double *decisionEnd, double *objectiveBegin, double *objectiveEnd, const size_t posGroups, const size_t distDecisions);
void WFG3(double *decisionBegin, double *decisionEnd, double *objectiveBegin, double *objectiveEnd, const size_t posGroups, const size_t distDecisions);
void WFG4(double *decisionBegin, double *decisionEnd, double *objectiveBegin, double *objectiveEnd, const size_t posGroups, const size_t distDecisions);
void WFG5(double *decisionBegin, double *decisionEnd, double *objectiveBegin, double *objectiveEnd, const size_t posGroups, const size_t distDecisions);
void WFG6(double *decisionBegin, double *decisionEnd, double *objectiveBegin, double *objectiveEnd, const size_t posGroups, const size_t distDecisions);
void WFG7(double *decisionBegin, double *decisionEnd, double *objectiveBegin, double *objectiveEnd, const size_t posGroups, const size_t distDecisions);
void WFG8(double *decisionBegin, double *decisionEnd, double *objectiveBegin, double *objectiveEnd, const size_t posGroups, const size_t distDecisions);
void WFG9(double *decisionBegin, double *decisionEnd, double *objectiveBegin, double *objectiveEnd, const size_t posGroups, const size_t distDecisions);



void WFG_Calc(short number, const std::vector<double> & decision_var, std::vector<double> & objectives)
{
	//copy the decision variables for safety
	std::vector<double> temp_variables = decision_var;

	//execute the right code
	switch (number)
	{
	case 1:
		WFG1(&temp_variables.front(), &temp_variables.back() + 1, &objectives.front(), &objectives.back() + 1, 2, 20);
		break;
	case 2:
		WFG2(&temp_variables.front(), &temp_variables.back() + 1, &objectives.front(), &objectives.back() + 1, 2, 20);
		break;
	case 3:
		WFG3(&temp_variables.front(), &temp_variables.back() + 1, &objectives.front(), &objectives.back() + 1, 2, 20);
		break;
	case 4:
		WFG4(&temp_variables.front(), &temp_variables.back() + 1, &objectives.front(), &objectives.back() + 1, 2, 20);
		break;
	case 5:
		WFG5(&temp_variables.front(), &temp_variables.back() + 1, &objectives.front(), &objectives.back() + 1, 2, 20);
		break;
	case 6:
		WFG6(&temp_variables.front(), &temp_variables.back() + 1, &objectives.front(), &objectives.back() + 1, 2, 20);
		break;
	case 7:
		WFG7(&temp_variables.front(), &temp_variables.back() + 1, &objectives.front(), &objectives.back() + 1, 2, 20);
		break;
	case 8:
		WFG8(&temp_variables.front(), &temp_variables.back() + 1, &objectives.front(), &objectives.back() + 1, 2, 20);
		break;
	case 9:
		WFG9(&temp_variables.front(), &temp_variables.back() + 1, &objectives.front(), &objectives.back() + 1, 2, 20);
		break;
	default:
		abort();
	}	

}



void Normalize(double *begin, double *end)
{
	double *i;
	for(i = begin; i != end; ++i)
	{
		double bound = 2 * (i - begin + 1);
		assert(0 <= *i && *i <= bound);
		*i = *i / bound;
	}
}

void Scale(double distance, double *begin, double *end)
{
	double *i;
	for(i = begin; i != end; ++i)
		*i = distance + 2 * (i - begin + 1) * *i;
}

void WFG1Transition1(double *begin, double *end, size_t nPosDecisions)
{
	double *i;
	assert(0 < nPosDecisions && nPosDecisions <= end - begin);
	for(i = begin + nPosDecisions; i != end; ++i)
		*i = LinearTransform(*i, 0.35);
}

void WFG1Transition2(double *begin, double *end, size_t nPosDecisions)
{
	double *i;
	assert(0 < nPosDecisions && nPosDecisions <= end - begin);
	for(i = begin + nPosDecisions; i != end; ++i)
		*i = FlatTransform(*i, 0.8, 0.75, 0.85);
}

void WFG1Transition3(double *begin, double *end)
{
	double *i;
	for(i = begin; i != end; ++i)
		*i = PolynomialTransform(*i, 0.02);
}

double WFG1Transition4(double *decisionBegin, double *decisionEnd, size_t nPosDecisions, double *posDecisionBegin, double *posDecisionEnd)
{
	size_t decisionSize, posDecisionSize, nDistDecisions, i, nUnit;
	double *weight, *posDecision, result;

	decisionSize = decisionEnd - decisionBegin;
	posDecisionSize = posDecisionEnd - posDecisionBegin;
	assert(0 < nPosDecisions && nPosDecisions <= decisionSize);
	assert(posDecisionEnd - posDecisionBegin > 0);
	assert(nPosDecisions >= posDecisionSize && nPosDecisions % posDecisionSize == 0);
	nDistDecisions = decisionSize - nPosDecisions;

	weight = (double *)malloc(decisionSize * sizeof(double));
	for(i = 0; i < decisionSize; ++i)
		weight[i] = 2 * (i + 1);
	nUnit = nPosDecisions / posDecisionSize;
	for(posDecision = posDecisionBegin; posDecision != posDecisionEnd; ++posDecision)
	{
		size_t i = posDecision - posDecisionBegin;
		size_t begin = i * nUnit;
		size_t end = (i + 1) * nUnit;
		*posDecision = WeightedSumTransform(weight + begin, weight + end, decisionBegin + begin);
	}
	if (nDistDecisions > 0)
		result = WeightedSumTransform(weight + nPosDecisions, weight + decisionSize, decisionBegin + nPosDecisions);
	else
		result = 0;
	free(weight);
	return result;
}

void WFG1Shape(double *decisionBegin, double *decisionEnd, double distance, double *objectiveBegin, double *objectiveEnd)
{
	size_t decisionSize, i;
	double *angle, *degenerate;

	decisionSize = decisionEnd - decisionBegin;
	assert(decisionSize > 0);
	angle = (double *)malloc(decisionSize * sizeof(double));
	memcpy(angle, decisionBegin, decisionSize * sizeof(double));
	degenerate = (double *)malloc(decisionSize * sizeof(double));
	for (i = 0; i < decisionSize; ++i)
		degenerate[i] = 1;
	DegenerateTransform(angle, angle + decisionSize, degenerate, degenerate + decisionSize, distance);
	ConvertAngles(angle, angle + decisionSize);
	ConvexShape(angle, angle + decisionSize, objectiveBegin, objectiveEnd - 1, 1);
	*(objectiveEnd - 1) = ConvexConcaveShape(*decisionBegin, 5, 1);
	Scale(distance, objectiveBegin, objectiveEnd);
	free(angle);
	free(degenerate);
}

void WFG1(double *decisionBegin, double *decisionEnd, double *objectiveBegin, double *objectiveEnd, const size_t posGroups, const size_t distDecisions)
{
	size_t decisionSize, nObjectives, posDecisions, posDecisionSize;
	double *_decision, *posDecision, distance;

	decisionSize = decisionEnd - decisionBegin;
	_decision = (double *)malloc(decisionSize * sizeof(double));
	memcpy(_decision, decisionBegin, decisionSize * sizeof(double));
	Normalize(_decision, _decision + decisionSize);
	nObjectives = objectiveEnd - objectiveBegin;
	posDecisions = posGroups * (nObjectives - 1);
	WFG1Transition1(_decision, _decision + decisionSize, posDecisions);
	WFG1Transition2(_decision, _decision + decisionSize, posDecisions);
	WFG1Transition3(_decision, _decision + decisionSize);
	posDecisionSize = nObjectives - 1;
	posDecision = (double *)malloc(posDecisionSize * sizeof(double));
	distance = WFG1Transition4(_decision, _decision + decisionSize, posDecisions, posDecision, posDecision + posDecisionSize);
	WFG1Shape(posDecision, posDecision + posDecisionSize, distance, objectiveBegin, objectiveEnd);
	free(_decision);
	free(posDecision);
}

void WFG2Transition1(double *begin, double *end, size_t nPosDecisions)
{
	WFG1Transition1(begin, end, nPosDecisions);
}

size_t WFG2Transition2(double *begin, double *end, size_t nPosDecisions)
{
	size_t nDistDecisions, nDecisions, i;
	assert(0 < nPosDecisions && nPosDecisions <= end - begin);
	nDistDecisions = end - begin - nPosDecisions;
	nDecisions = nPosDecisions + nDistDecisions / 2;
	for(i = nPosDecisions; i < nDecisions; ++i)
	{
		const size_t _begin = nPosDecisions + 2 * (i - nPosDecisions);
		const size_t _end = nPosDecisions + 2 * (i - nPosDecisions + 1);
		begin[i] = NonSeparableTransform(begin + _begin, begin + _end, 2);
	}
	return nDecisions;
}

double WFG2Transition3(double *decisionBegin, double *decisionEnd, size_t nPosDecisions, double *posDecisionBegin, double *posDecisionEnd)
{
	size_t nDistDecisions, i, nUnit;
	double *weight, result;
	size_t decisionSize = decisionEnd - decisionBegin;
	size_t posDecisionSize = posDecisionEnd - posDecisionBegin;
	assert(0 < nPosDecisions && nPosDecisions <= decisionSize);
	assert(posDecisionSize > 0);
	assert(nPosDecisions >= posDecisionSize && nPosDecisions % posDecisionSize == 0);
	nDistDecisions = decisionSize - nPosDecisions;

	weight = (double *)malloc(decisionSize * sizeof(double));
	for(i = 0; i < decisionSize; ++i)
		weight[i] = 1;
	nUnit = nPosDecisions / posDecisionSize;
	for(i = 0; i < posDecisionSize; ++i)
	{
		const size_t begin = i * nUnit;
		const size_t end = (i + 1) * nUnit;
		posDecisionBegin[i] = WeightedSumTransform(weight + begin, weight + end, decisionBegin + begin);
	}
	if (nDistDecisions > 0)
		result = WeightedSumTransform(weight + nPosDecisions, weight + decisionSize, decisionBegin + nPosDecisions);
	else
		result = 0;
	free(weight);
	return result;
}

void WFG2Shape(double *decisionBegin, double *decisionEnd, double distance, double *objectiveBegin, double *objectiveEnd)
{
	size_t i;
	double *angle, *degenerate;

	size_t decisionSize = decisionEnd - decisionBegin;
	assert(decisionSize > 0);
	angle = (double *)malloc(decisionSize * sizeof(double));
	memcpy(angle, decisionBegin, decisionSize * sizeof(double));
	degenerate = (double *)malloc(decisionSize * sizeof(double));
	for (i = 0; i < decisionSize; ++i)
		degenerate[i] = 1;
	DegenerateTransform(angle, angle + decisionSize, degenerate, degenerate + decisionSize, distance);
	ConvertAngles(angle, angle + decisionSize);
	ConvexShape(angle, angle + decisionSize, objectiveBegin, objectiveEnd - 1, 1);
	*(objectiveEnd - 1) = DisconnectedShape(*decisionBegin, 5, 1, 1);
	Scale(distance, objectiveBegin, objectiveEnd);
	free(angle);
	free(degenerate);
}

void WFG2(double *decisionBegin, double *decisionEnd, double *objectiveBegin, double *objectiveEnd, const size_t posGroups, const size_t distDecisions)
{
	size_t nObjectives, posDecisions, posDecisionSize;
	double *posDecision, distance;

	size_t decisionSize = decisionEnd - decisionBegin;
	double *_decision = (double *)malloc(decisionSize * sizeof(double));
	memcpy(_decision, decisionBegin, decisionSize * sizeof(double));
	Normalize(_decision, _decision + decisionSize);
	nObjectives = objectiveEnd - objectiveBegin;
	posDecisions = posGroups * (nObjectives - 1);
	WFG2Transition1(_decision, _decision + decisionSize, posDecisions);
	decisionSize = WFG2Transition2(_decision, _decision + decisionSize, posDecisions);
	posDecisionSize = nObjectives - 1;
	posDecision = (double *)malloc(posDecisionSize * sizeof(double));
	distance = WFG2Transition3(_decision, _decision + decisionSize, posDecisions, posDecision, posDecision + posDecisionSize);
	WFG2Shape(posDecision, posDecision + posDecisionSize, distance, objectiveBegin, objectiveEnd);
	free(_decision);
	free(posDecision);
}

void WFG3Transition1(double *begin, double *end, size_t nPosDecisions)
{
	WFG2Transition1(begin, end, nPosDecisions);
}

size_t WFG3Transition2(double *begin, double *end, size_t nPosDecisions)
{
	return WFG2Transition2(begin, end, nPosDecisions);
}

double WFG3Transition3(double *decisionBegin, double *decisionEnd, size_t nPosDecisions, double *posDecisionBegin, double *posDecisionEnd)
{
	return WFG2Transition3(decisionBegin, decisionEnd, nPosDecisions, posDecisionBegin, posDecisionEnd);
}

void WFG3Shape(double *decisionBegin, double *decisionEnd, double distance, double *objectiveBegin, double *objectiveEnd)
{
	size_t i;
	double *_decision, *degenerate;

	size_t decisionSize = decisionEnd - decisionBegin;
	assert(decisionSize > 0);
	_decision = (double *)malloc(decisionSize * sizeof(double));
	memcpy(_decision, decisionBegin, decisionSize * sizeof(double));
	degenerate = (double *)malloc(decisionSize * sizeof(double));
	degenerate[0] = 1;
	for (i = 1; i < decisionSize; ++i)
		degenerate[i] = 0;
	DegenerateTransform(_decision, _decision + decisionSize, degenerate, degenerate + decisionSize, distance);
	LinearShape(_decision, _decision + decisionSize, objectiveBegin, objectiveEnd, 1);
	Scale(distance, objectiveBegin, objectiveEnd);
	free(_decision);
	free(degenerate);
}

void WFG3(double *decisionBegin, double *decisionEnd, double *objectiveBegin, double *objectiveEnd, const size_t posGroups, const size_t distDecisions)
{
	size_t nObjectives, posDecisions, posDecisionSize;
	double *posDecision, distance;

	size_t decisionSize = decisionEnd - decisionBegin;
	double *_decision = (double *)malloc(decisionSize * sizeof(double));
	memcpy(_decision, decisionBegin, decisionSize * sizeof(double));
	Normalize(_decision, _decision + decisionSize);
	nObjectives = objectiveEnd - objectiveBegin;
	posDecisions = posGroups * (nObjectives - 1);
	WFG3Transition1(_decision, _decision + decisionSize, posDecisions);
	decisionSize = WFG3Transition2(_decision, _decision + decisionSize, posDecisions);
	posDecisionSize = nObjectives - 1;
	posDecision = (double *)malloc(posDecisionSize * sizeof(double));
	distance = WFG3Transition3(_decision, _decision + decisionSize, posDecisions, posDecision, posDecision + posDecisionSize);
	WFG3Shape(posDecision, posDecision + posDecisionSize, distance, objectiveBegin, objectiveEnd);
	free(_decision);
	free(posDecision);
}

void WFG4Transition1(double *begin, double *end)
{
	double *i;
	for(i = begin; i != end; ++i)
		*i = MultiModalTransform(*i, 30, 10, 0.35);
}

double WFG4Transition2(double *decisionBegin, double *decisionEnd, size_t nPosDecisions, double *posDecisionBegin, double *posDecisionEnd)
{
	return WFG2Transition3(decisionBegin, decisionEnd, nPosDecisions, posDecisionBegin, posDecisionEnd);
}

void WFG4Shape(double *decisionBegin, double *decisionEnd, double distance, double *objectiveBegin, double *objectiveEnd)
{
	size_t i;
	double *angle, *degenerate;

	size_t decisionSize = decisionEnd - decisionBegin;
	assert(decisionSize > 0);
	angle = (double *)malloc(decisionSize * sizeof(double));
	memcpy(angle, decisionBegin, decisionSize * sizeof(double));
	degenerate = (double *)malloc(decisionSize * sizeof(double));
	for (i = 0; i < decisionSize; ++i)
		degenerate[i] = 1;
	DegenerateTransform(angle, angle + decisionSize, degenerate, degenerate + decisionSize, distance);
	ConvertAngles(angle, angle + decisionSize);
	InvertConcaveShape(angle, angle + decisionSize, objectiveBegin, objectiveEnd, 1);
	Scale(distance, objectiveBegin, objectiveEnd);
	free(angle);
	free(degenerate);
}

void WFG4(double *decisionBegin, double *decisionEnd, double *objectiveBegin, double *objectiveEnd, const size_t posGroups, const size_t distDecisions)
{
	size_t nObjectives, posDecisions, posDecisionSize;
	double *posDecision, distance;

	size_t decisionSize = decisionEnd - decisionBegin;
	double *_decision = (double *)malloc(decisionSize * sizeof(double));
	memcpy(_decision, decisionBegin, decisionSize * sizeof(double));
	Normalize(_decision, _decision + decisionSize);
	nObjectives = objectiveEnd - objectiveBegin;
	posDecisions = posGroups * (nObjectives - 1);
	WFG4Transition1(_decision, _decision + decisionSize);
	posDecisionSize = nObjectives - 1;
	posDecision = (double *)malloc(posDecisionSize * sizeof(double));
	distance = WFG4Transition2(_decision, _decision + decisionSize, posDecisions, posDecision, posDecision + posDecisionSize);
	WFG4Shape(posDecision, posDecision + posDecisionSize, distance, objectiveBegin, objectiveEnd);
	free(_decision);
	free(posDecision);
}

void WFG5Transition1(double *begin, double *end)
{
	double *i;
	for(i = begin; i != end; ++i)
		*i = DeceptiveTransform(*i, 0.35, 0.001, 0.05);
}

double WFG5Transition2(double *decisionBegin, double *decisionEnd, size_t nPosDecisions, double *posDecisionBegin, double *posDecisionEnd)
{
	return WFG2Transition3(decisionBegin, decisionEnd, nPosDecisions, posDecisionBegin, posDecisionEnd);
}

void WFG5Shape(double *decisionBegin, double *decisionEnd, double distance, double *objectiveBegin, double *objectiveEnd)
{
	WFG4Shape(decisionBegin, decisionEnd, distance, objectiveBegin, objectiveEnd);
}

void WFG5(double *decisionBegin, double *decisionEnd, double *objectiveBegin, double *objectiveEnd, const size_t posGroups, const size_t distDecisions)
{
	size_t nObjectives, posDecisions, posDecisionSize;
	double *posDecision, distance;

	size_t decisionSize = decisionEnd - decisionBegin;
	double *_decision = (double *)malloc(decisionSize * sizeof(double));
	memcpy(_decision, decisionBegin, decisionSize * sizeof(double));
	Normalize(_decision, _decision + decisionSize);
	nObjectives = objectiveEnd - objectiveBegin;
	posDecisions = posGroups * (nObjectives - 1);
	WFG5Transition1(_decision, _decision + decisionSize);
	posDecisionSize = nObjectives - 1;
	posDecision = (double *)malloc(posDecisionSize * sizeof(double));
	distance = WFG5Transition2(_decision, _decision + decisionSize, posDecisions, posDecision, posDecision + posDecisionSize);
	WFG5Shape(posDecision, posDecision + posDecisionSize, distance, objectiveBegin, objectiveEnd);
	free(_decision);
	free(posDecision);
}

void WFG6Transition1(double *begin, double *end, size_t nPosDecisions)
{
	WFG1Transition1(begin, end, nPosDecisions);
}

double WFG6Transition2(double *decisionBegin, double *decisionEnd, size_t nPosDecisions, double *posDecisionBegin, double *posDecisionEnd)
{
	size_t nDistDecisions, nUnit, i;
	double result;

	size_t decisionSize = decisionEnd - decisionBegin;
	size_t posDecisionSize = posDecisionEnd - posDecisionBegin;
	assert(0 < nPosDecisions && nPosDecisions <= decisionSize);
	assert(posDecisionSize > 0);
	assert(nPosDecisions >= posDecisionSize && nPosDecisions % posDecisionSize == 0);
	nDistDecisions = decisionSize - nPosDecisions;

	nUnit = nPosDecisions / posDecisionSize;
	for(i = 0; i < posDecisionSize; ++i)
	{
		const size_t begin = i * nUnit;
		const size_t end = (i + 1) * nUnit;
		posDecisionBegin[i] = NonSeparableTransform(decisionBegin + begin, decisionBegin + end, nUnit);
	}
	if (nDistDecisions > 0)
		result = NonSeparableTransform(decisionBegin + nPosDecisions, decisionBegin + decisionSize, decisionSize - nPosDecisions);
	else
		result = 0;
	return result;
}

void WFG6Shape(double *decisionBegin, double *decisionEnd, double distance, double *objectiveBegin, double *objectiveEnd)
{
	WFG4Shape(decisionBegin, decisionEnd, distance, objectiveBegin, objectiveEnd);
}

void WFG6(double *decisionBegin, double *decisionEnd, double *objectiveBegin, double *objectiveEnd, const size_t posGroups, const size_t distDecisions)
{
	size_t nObjectives, posDecisions, posDecisionSize;
	double *posDecision, distance;

	size_t decisionSize = decisionEnd - decisionBegin;
	double *_decision = (double *)malloc(decisionSize * sizeof(double));
	memcpy(_decision, decisionBegin, decisionSize * sizeof(double));
	Normalize(_decision, _decision + decisionSize);
	nObjectives = objectiveEnd - objectiveBegin;
	posDecisions = posGroups * (nObjectives - 1);
	WFG6Transition1(_decision, _decision + decisionSize, posDecisions);
	posDecisionSize = nObjectives - 1;
	posDecision = (double *)malloc(posDecisionSize * sizeof(double));
	distance = WFG6Transition2(_decision, _decision + decisionSize, posDecisions, posDecision, posDecision + posDecisionSize);
	WFG6Shape(posDecision, posDecision + posDecisionSize, distance, objectiveBegin, objectiveEnd);
	free(_decision);
	free(posDecision);
}

void WFG7Transition1(double *begin, double *end, size_t nPosDecisions)
{
	size_t i;
	double *weight;

	size_t decisionSize = end - begin;
	assert(0 < nPosDecisions && nPosDecisions <= decisionSize);
	weight = (double *)malloc(decisionSize * sizeof(double));
	for(i = 0; i < decisionSize; ++i)
		weight[i] = 1;
	for(i = 0; i < nPosDecisions; ++i)
	{
		size_t _begin = i + 1;
		if (_begin < decisionSize)
		{
			double regionMapping = WeightedSumTransform(weight + _begin, weight + decisionSize, begin + _begin);
			begin[i] = DependentTransform(begin[i], regionMapping, 0.98 / 49.98, 0.02, 50);
		}
	}
	free(weight);
}

void WFG7Transition2(double *begin, double *end, size_t nPosDecisions)
{
	WFG1Transition1(begin, end, nPosDecisions);
}

double WFG7Transition3(double *decisionBegin, double *decisionEnd, size_t nPosDecisions, double *posDecisionBegin, double *posDecisionEnd)
{
	return WFG2Transition3(decisionBegin, decisionEnd, nPosDecisions, posDecisionBegin, posDecisionEnd);
}

void WFG7Shape(double *decisionBegin, double *decisionEnd, double distance, double *objectiveBegin, double *objectiveEnd)
{
	WFG4Shape(decisionBegin, decisionEnd, distance, objectiveBegin, objectiveEnd);
}

void WFG7(double *decisionBegin, double *decisionEnd, double *objectiveBegin, double *objectiveEnd, const size_t posGroups, const size_t distDecisions)
{
	size_t nObjectives, posDecisions, posDecisionSize;
	double *posDecision, distance;

	size_t decisionSize = decisionEnd - decisionBegin;
	double *_decision = (double *)malloc(decisionSize * sizeof(double));
	memcpy(_decision, decisionBegin, decisionSize * sizeof(double));
	Normalize(_decision, _decision + decisionSize);
	nObjectives = objectiveEnd - objectiveBegin;
	posDecisions = posGroups * (nObjectives - 1);
	WFG7Transition1(_decision, _decision + decisionSize, posDecisions);
	WFG7Transition2(_decision, _decision + decisionSize, posDecisions);
	posDecisionSize = nObjectives - 1;
	posDecision = (double *)malloc(posDecisionSize * sizeof(double));
	distance = WFG7Transition3(_decision, _decision + decisionSize, posDecisions, posDecision, posDecision + posDecisionSize);
	WFG7Shape(posDecision, posDecision + posDecisionSize, distance, objectiveBegin, objectiveEnd);
	free(_decision);
	free(posDecision);
}

void WFG8Transition1(double *begin, double *end, size_t nPosDecisions)
{
	size_t i;
	double *weight, *_decision;

	size_t decisionSize = end - begin;
	assert(0 < nPosDecisions && nPosDecisions <= decisionSize);
	weight = (double *)malloc(decisionSize * sizeof(double));
	for(i = 0; i < decisionSize; ++i)
		weight[i] = 1;
	_decision = (double *)malloc(decisionSize * sizeof(double));
	memcpy(_decision, begin, decisionSize * sizeof(double));
	for(i = nPosDecisions; i < decisionSize; ++i)
	{
		double regionMapping = WeightedSumTransform(weight, weight + i, _decision);
		begin[i] = DependentTransform(_decision[i], regionMapping, 0.98 / 49.98, 0.02, 50);
	}
	free(weight);
	free(_decision);
}

void WFG8Transition2(double *begin, double *end, size_t nPosDecisions)
{
	WFG1Transition1(begin, end, nPosDecisions);
}

double WFG8Transition3(double *decisionBegin, double *decisionEnd, size_t nPosDecisions, double *posDecisionBegin, double *posDecisionEnd)
{
	return WFG2Transition3(decisionBegin, decisionEnd, nPosDecisions, posDecisionBegin, posDecisionEnd);
}

void WFG8Shape(double *decisionBegin, double *decisionEnd, double distance, double *objectiveBegin, double *objectiveEnd)
{
	WFG4Shape(decisionBegin, decisionEnd, distance, objectiveBegin, objectiveEnd);
}

void WFG8(double *decisionBegin, double *decisionEnd, double *objectiveBegin, double *objectiveEnd, const size_t posGroups, const size_t distDecisions)
{
	size_t nObjectives, posDecisions, posDecisionSize;
	double *posDecision, distance;

	size_t decisionSize = decisionEnd - decisionBegin;
	double *_decision = (double *)malloc(decisionSize * sizeof(double));
	memcpy(_decision, decisionBegin, decisionSize * sizeof(double));
	Normalize(_decision, _decision + decisionSize);
	nObjectives = objectiveEnd - objectiveBegin;
	posDecisions = posGroups * (nObjectives - 1);
	WFG8Transition1(_decision, _decision + decisionSize, posDecisions);
	WFG8Transition2(_decision, _decision + decisionSize, posDecisions);
	posDecisionSize = nObjectives - 1;
	posDecision = (double *)malloc(posDecisionSize * sizeof(double));
	distance = WFG8Transition3(_decision, _decision + decisionSize, posDecisions, posDecision, posDecision + posDecisionSize);
	WFG8Shape(posDecision, posDecision + posDecisionSize, distance, objectiveBegin, objectiveEnd);
	free(_decision);
	free(posDecision);
}

void WFG9Transition1(double *begin, double *end)
{
	size_t i;
	size_t decisionSize = end - begin;
	double *weight = (double *)malloc(decisionSize * sizeof(double));
	for(i = 0; i < decisionSize; ++i)
		weight[i] = 1;
	for(i = 0; i < decisionSize - 1; ++i)
	{
		size_t _begin = i + 1;
		double regionMapping = WeightedSumTransform(weight + _begin, weight + decisionSize, begin + _begin);
		begin[i] = DependentTransform(begin[i], regionMapping, 0.98 / 49.98, 0.02, 50);
	}
	free(weight);
}

void WFG9Transition2(double *begin, double *end, size_t nPosDecisions)
{
	size_t i;
	size_t decisionSize = end - begin;
	assert(0 < nPosDecisions && nPosDecisions <= decisionSize);
	for(i = 0; i < nPosDecisions; ++i)
		begin[i] = DeceptiveTransform(begin[i], 0.35, 0.001, 0.05);
	for(i = nPosDecisions; i < decisionSize; ++i)
		begin[i] = MultiModalTransform(begin[i], 30, 95, 0.35);
}

double WFG9Transition3(double *decisionBegin, double *decisionEnd, size_t nPosDecisions, double *posDecisionBegin, double *posDecisionEnd)
{
	return WFG6Transition2(decisionBegin, decisionEnd, nPosDecisions, posDecisionBegin, posDecisionEnd);
}

void WFG9Shape(double *decisionBegin, double *decisionEnd, double distance, double *objectiveBegin, double *objectiveEnd)
{
	WFG4Shape(decisionBegin, decisionEnd, distance, objectiveBegin, objectiveEnd);
}

void WFG9(double *decisionBegin, double *decisionEnd, double *objectiveBegin, double *objectiveEnd, const size_t posGroups, const size_t distDecisions)
{
	size_t nObjectives, posDecisions, posDecisionSize;
	double *posDecision, distance;

	size_t decisionSize = decisionEnd - decisionBegin;
	double *_decision = (double *)malloc(decisionSize * sizeof(double));
	memcpy(_decision, decisionBegin, decisionSize * sizeof(double));
	Normalize(_decision, _decision + decisionSize);
	nObjectives = objectiveEnd - objectiveBegin;
	posDecisions = posGroups * (nObjectives - 1);
	WFG9Transition1(_decision, _decision + decisionSize);
	WFG9Transition2(_decision, _decision + decisionSize, posDecisions);
	posDecisionSize = nObjectives - 1;
	posDecision = (double *)malloc(posDecisionSize * sizeof(double));
	distance = WFG9Transition3(_decision, _decision + decisionSize, posDecisions, posDecision, posDecision + posDecisionSize);
	WFG9Shape(posDecision, posDecision + posDecisionSize, distance, objectiveBegin, objectiveEnd);
	free(_decision);
	free(posDecision);
}
