#define _CRT_SECURE_NO_WARNINGS

#include "TGM.h"

extern std::vector<float> TGM_vect;												//current grandmaternal effect

void TGM_Calc(std::vector<std::vector<double>> &TGM_fit, const std::vector<double> & fit)
{
	int fit_size = fit.size();
	std::vector<double> new_fit = fit;
	for (int i_tgm = 0; i_tgm < TGM_size; i_tgm++)
	{
		for (int i_fit = 0; i_fit < fit_size; i_fit++)
			new_fit[i_fit] += TGM_fit[i_tgm + 1][i_fit] * TGM_vect[i_tgm];
	}
	TGM_fit[0] = new_fit;
}