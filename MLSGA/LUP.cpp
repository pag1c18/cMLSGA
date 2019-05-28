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


#include <math.h>
#include "LUP.h"





std::vector<double>  Decompose_Solve(double **A, double *b, int N);


std::vector<double> LU_Solve(std::vector<std::vector<double>> matrix_A, std::vector<double> matrix_B)
{
	std::vector<double> x_out;
	int mat_size = matrix_A.size();

	if (mat_size != matrix_B.size() || mat_size != matrix_A[0].size())
		abort();

	double **A = (double**)malloc(sizeof(double*) * mat_size);

	double *b = (double*)malloc(sizeof(double) * mat_size);


	for (int i = 0; i < mat_size; i++)
	{
		A[i] = (double*)malloc(sizeof(double) * mat_size);

		b[i] = matrix_B[i];


		for (int j = 0; j < mat_size; j++)
			A[i][j] = matrix_A[i][j];
	}
	


	x_out = Decompose_Solve(A, b, mat_size);

	for (int i = 0; i < mat_size; i++)
	{
		free(A[i]);
	}
	free(A);
	free(b);

	return x_out;

}





std::vector<double>  Decompose_Solve(double **A, double *b, int N)
{
	int *P = (int*)malloc(sizeof(int) * N);

	double *x = (double*)malloc(sizeof(double) * N);

	double *ptr;

	for (int i_obj = 0; i_obj <= N; i_obj++)
		P[i_obj] = i_obj;

	for (int i_obj = 0; i_obj < N; i_obj++) 
	{


		double maxA = 0.;
		int imax = i_obj;

		for (int k = i_obj; k < N; k++)
		{
			if (abs(A[k][i_obj]) > maxA)
			{
				maxA = abs(A[k][i_obj]);
				imax = k;
			}
		}


		if (imax != i_obj) 
		{
			int j = P[i_obj];
			P[i_obj] = P[imax];
			P[imax] = j;

			ptr = A[i_obj];
			A[i_obj] = A[imax];
			A[imax] = ptr;

			P[N]++;
		}

		for (int j = i_obj + 1; j < N; j++) 
		{
			A[j][i_obj] /= A[i_obj][i_obj];

			for (int k = i_obj + 1; k < N; k++)
				A[j][k] -= A[j][i_obj] * A[i_obj][k];
		}
	}

	std::vector<double> x_out;
	for (int i_obj = 0; i_obj < N; i_obj++) 
	{
		x[i_obj] = b[P[i_obj]];

		for (int i_obj_2 = 0; i_obj_2 < i_obj; i_obj_2++)
			x[i_obj] -= A[i_obj][i_obj_2] * x[i_obj_2];
	}

	for (int i_obj = N - 1; i_obj >= 0; i_obj--) 
	{
		for (int i_obj_2 = i_obj + 1; i_obj_2 < N; i_obj_2++)
			x[i_obj] -= A[i_obj][i_obj_2] * x[i_obj_2];

		x[i_obj] = x[i_obj] / A[i_obj][i_obj];
	}

	for (int i_obj = 0; i_obj < N; i_obj++)
	{
		x_out.push_back(x[i_obj]);
	}
	free(x);
	//free(P);

	return x_out;

}
