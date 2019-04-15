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