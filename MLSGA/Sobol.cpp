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


#include "Sobol.h"
#include <fstream>
#include "Random.h"
#include <sstream>


std::vector<std::vector<double>> Sobol_Sequence(unsigned size, unsigned dim)
{
	std::fstream file;
	file.open("Input/dane.txt");
	std::ofstream file2;
	file2.open("Input/out.csv");
	std::vector<unsigned> C;
	unsigned s = 15;			//number of bits
	std::vector<std::vector<double>> output;	//output vector
	std::vector<std::vector<unsigned>> m;

	//Check if size not exceed max nubmber of bits
	if ((unsigned)ceil(log((double)size) / log(2.0)) > s)
		abort();

	//Calculate C values
	C.push_back(1);
	for (unsigned i = 1; i < size; i++)
	{
		C.push_back(1);
		unsigned val = i;
		while (val & 1)
		{
			val >>= 1;
			C[i]++;
		}
	}
	for (unsigned i = 0; i < 1800; i++)
	{
		std::vector<unsigned> m_temp;
		for (unsigned j = 0; j < s; j++)
		{
			unsigned temp;
			file >> temp;
			m_temp.push_back(temp);
		}
		m.push_back(m_temp);
	}

	//Caclculate the sobol sequence
	for (unsigned i = 0; i < dim; i++)
	{
		std::vector<double> temp_vect;		//temporary vector for storage of output values
		std::vector<unsigned> v, x;
		int r_i = Random_I(0, 1799);
		for (unsigned j = 0; j < s; j++)
		{
			v.push_back(m[r_i][j] << (32 - (j + 1)));
		}
		x.push_back(0);
		temp_vect.push_back(0.0);
		file2 << (double)x[0] / pow(2.0, 32) << ", ";
		for (unsigned j = 1; j < size; j++)
		{
			x.push_back(x[j - 1] ^ v[C[j - 1]-1]);
			temp_vect.push_back((double)x[j] / pow(2.0, 32));
			file2 << (double)x[j] / pow(2.0, 32) << ", ";
		}
		file2 << std::endl;
		output.push_back(temp_vect);
	}
	file.close();
	file2.close();
	return output;
}