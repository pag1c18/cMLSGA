#pragma once
//****************************************
//			CONTOUR_PLOT header
//		Function for conrour plot making
//****************************************
#ifndef CONTOUR_PLOT_H
#define CONTOUR_PLOT_H


/*
*Contour plot generation*
@param indexr index of the run for the current contour plot
*/
void Contour_Plot(int indexr);
/**returning max z value for the contour plot**/
double Get_Max_CP();
#endif // !CONTOUR_PLOT_H
