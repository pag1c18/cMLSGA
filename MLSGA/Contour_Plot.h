/**
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
along with this program.If not, see < https://www.gnu.org/licenses/>. 

		CONTOUR_PLOT header
		Function for contour plot making


*/

#pragma once

#ifndef CONTOUR_PLOT_H
#define CONTOUR_PLOT_H


/**
*Contour plot generation
@param indexr - index of the run for the current contour plot
@param t - current time step. For dynamic functions.
*/
void Contour_Plot(int indexr, double t);
/**Calculate the max z value for the contour plot*/
double Get_Max_CP();
#endif // !CONTOUR_PLOT_H
