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

#pragma once

#include <stddef.h>

void DegenerateTransform(double *decisionBegin, double *decisionEnd, double *degenerateBegin, double *degenerateEnd, double distance);
//Bias
double FlatTransform(double decision, double region, double regionMin, double regionMax);
double DependentTransform(double decision, double regionMapping, double regionFraction, double regionMin, double regionMax);
double PolynomialTransform(double decision, double bias);
//Reduction
double NonSeparableTransform(double *begin, double *end, size_t degree);
double WeightedSumTransform(double *weightBegin, double *weightEnd, double *decisionBegin);
//Shift
double LinearTransform(double decision, double globalMin);
double DeceptiveTransform(double decision, double globalMin, double aperture, double deceptiveMin);
double MultiModalTransform(double decision, const size_t nMinima, double hillSize, double globalMin);
