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

void ConvertAngles(double *begin, double *end);
void ConvexShape(double *angleBegin, double *angleEnd, double *objectiveBegin, double *objectiveEnd, double radius);
double ConvexConcaveShape(double decision, size_t nSegments, double shape);
double DisconnectedShape(double decision, size_t nRegions, double shape, double location);
void InvertConcaveShape(double *angleBegin, double *angleEnd, double *objectiveBegin, double *objectiveEnd, double radius);
void LinearShape(double *decisionBegin, double *decisionEnd, double *objectiveBegin, double *objectiveEnd, double distance);
