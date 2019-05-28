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

#include "MLSGA_Add_Functions.h"
#include "Support_Functions.h"
#include <Windows.h>
#include <ctime>
#include "Contour_Plot.h"
#include "Imported\Workbook.h"

/*#include "imported.h"
#include "SVM.h"
#include "Contour_Plot.h"
#include "Video.h"
#include "Clustering.h"
#include "IGD.h"
#include "GA_Functions.h"
#include "Random_Label.h"*/



int cons_viol_count;		//How many times constrains have been violated in the current run
extern short n_col;
extern int pop_size;

extern std::vector<float> TGM_vect;												//current grandmaternal effect


function * Func_Type(short ix)
{
	//1-6 ZDT 1-6
	if (ix == 1)
	{
		ZDT1 * zdt1 = new ZDT1;
		return zdt1;
	}
	else if (ix == 2)
	{
		ZDT2 * zdt2 = new ZDT2;
		return zdt2;
	}
	else if (ix == 3)
	{
		ZDT3 * zdt3 = new ZDT3;
		return zdt3;
	}
	else if (ix == 4)
	{
		ZDT4 * zdt4 = new ZDT4;
		return zdt4;
	}
	else if (ix == 5)
	{
		//not implemented!!
		//ZDT5 * zdt5 = new ZDT5;
		//return zdt5;
		return NULL;
	}
	else if (ix == 6)
	{
		ZDT6 * zdt6 = new ZDT6;
		return zdt6;
	}
	//7-13 DTLZ 1-7
	else if (ix == 7)
	{
		DTLZ1 *dtlz1 = new DTLZ1;
		return dtlz1;
	}
	else if (ix == 8)
	{
		DTLZ2 *dtlz2 = new DTLZ2;
		return dtlz2;
	}
	else if (ix == 9)
	{
		DTLZ3 *dtlz3 = new DTLZ3;
		return dtlz3;
	}
	else if (ix == 10)
	{
		DTLZ4 *dtlz4 = new DTLZ4;
		return dtlz4;
	}
	else if (ix == 11)
	{
		DTLZ5 *dtlz5 = new DTLZ5;
		return dtlz5;
	}
	else if (ix == 12)
	{
		DTLZ6 *dtlz6 = new DTLZ6;
		return dtlz6;
	}
	else if (ix == 13)
	{
		DTLZ7 *dtlz7 = new DTLZ7;
		return dtlz7;
	}
	//14-20 MOP 1-7
	else if (ix == 14)
	{
		MOP1 *mop1 = new MOP1;
		return mop1;
	}
	else if (ix == 15)
	{
		MOP2 *mop2 = new MOP2;
		return mop2;
	}
	else if (ix == 16)
	{
		MOP3 *mop3 = new MOP3;
		return mop3;
	}
	else if (ix == 17)
	{
		MOP4 *mop4 = new MOP4;
		return mop4;
	}
	else if (ix == 18)
	{
		MOP5 *mop5 = new MOP5;
		return mop5;
	}
	else if (ix == 19)
	{
		MOP6 *mop6 = new MOP6;
		return mop6;
	}
	else if (ix == 20)
	{
		MOP7 *mop7 = new MOP7;
		return mop7;
	}
	//21-30 UF1-10
	else if (ix == 21)
	{
		UF1 *U1 = new UF1;
		return U1;
	}
	else if (ix == 22)
	{
		UF2 *U2 = new UF2;
		return U2;
	}
	else if (ix == 23)
	{
		UF3 *U3 = new UF3;
		return U3;
	}
	else if (ix == 24)
	{
		UF4 *U4 = new UF4;
		return U4;
	}
	else if (ix == 25)
	{
		UF5 *U5 = new UF5;
		return U5;
	}
	else if (ix == 26)
	{
		UF6 *U6 = new UF6;
		return U6;
	}
	else if (ix == 27)
	{
		UF7 *U7 = new UF7;
		return U7;
	}
	else if (ix == 28)
	{
		UF8 *U8 = new UF8;
		return U8;
	}
	else if (ix == 29)
	{
		UF9 *U9 = new UF9;
		return U9;
	}
	else if (ix == 30)
	{
		UF10 *U10 = new UF10;
		return U10;
	}
	//31-39 WFG 1-9
	else if (ix == 31)
	{
		WFG1 *wfg1 = new WFG1;
		return wfg1;

	}
	else if (ix == 32)
	{
		WFG2 *wfg2 = new WFG2;
		return wfg2;
	}
	else if (ix == 33)
	{
		WFG3 *wfg3 = new WFG3;
		return wfg3;
	}
	else if (ix == 34)
	{
		WFG4 *wfg4 = new WFG4;
		return wfg4;
	}
	else if (ix == 35)
	{
		WFG5 *wfg5 = new WFG5;
		return wfg5;
	}
	else if (ix == 36)
	{
		WFG6 *wfg6 = new WFG6;
		return wfg6;
	}
	else if (ix == 37)
	{
		WFG7 *wfg7 = new WFG7;
		return wfg7;
	}
	else if (ix == 38)
	{
		WFG8 *wfg8 = new WFG8;
		return wfg8;
	}
	else if (ix == 39)
	{
		WFG9 *wfg9 = new WFG9;
		return wfg9;
	}
	//40-49 IMB 1-10
	else if (ix == 40)
	{
		IMB1 *imb1 = new IMB1;
		return imb1;
	}
	else if (ix == 41)
	{
		IMB2 *imb2 = new IMB2;
		return imb2;
	}
	else if (ix == 42)
	{
		IMB3 *imb3 = new IMB3;
		return imb3;
	}
	else if (ix == 43)
	{
		IMB4 *imb4 = new IMB4;
		return imb4;
	}
	else if (ix == 44)
	{
		IMB5 *imb5 = new IMB5;
		return imb5;
	}
	else if (ix == 45)
	{
		IMB6 *imb6 = new IMB6;
		return imb6;
	}
	else if (ix == 46)
	{
		IMB7 *imb7 = new IMB7;
		return imb7;
	}
	else if (ix == 47)
	{
		IMB8 *imb8 = new IMB8;
		return imb8;
	}
	else if (ix == 48)
	{
		IMB9 *imb9 = new IMB9;
		return imb9;
	}
	else if (ix == 49)
	{
		IMB10 *imb10 = new IMB10;
		return imb10;
	}
	//50-59 CF 1-10
	else if (ix == 50)
	{
		CF1 *cf1 = new CF1;
		return cf1;
	}
	else if (ix == 51)
	{
		CF2 *cf2 = new CF2;
		return cf2;
	}
	else if (ix == 52)
	{
		CF3 *cf3 = new CF3;
		return cf3;
	}
	else if (ix == 53)
	{
		CF4 *cf4 = new CF4;
		return cf4;
	}
	else if (ix == 54)
	{
		CF5 *cf5 = new CF5;
		return cf5;
	}
	else if (ix == 55)
	{
		CF6 *cf6 = new CF6;
		return cf6;
	}
	else if (ix == 56)
	{
		CF7 *cf7 = new CF7;
		return cf7;
	}
	else if (ix == 57)
	{
		CF8 *cf8 = new CF8;
		return cf8;
	}
	else if (ix == 58)
	{
		CF9 *cf9 = new CF9;
		return cf9;
	}
	else if (ix == 59)
	{
		CF10 *cf10 = new CF10;
		return cf10;
	}
	//60-61 DTLZ 8-9
	else if (ix == 60)
	{
		DTLZ8 *dtlz8 = new DTLZ8;
		return dtlz8;
	}
	else if (ix == 61)
	{
		DTLZ9 *dtlz9 = new DTLZ9;
		return dtlz9;
	}
	//62-65 IMB 11-14
	else if (ix == 62)
	{
		IMB11 *imb11 = new IMB11;
		return imb11;
	}
	else if (ix == 63)
	{
		IMB12 *imb12 = new IMB12;
		return imb12;
	}
	else if (ix == 64)
	{
		IMB13 *imb13 = new IMB13;
		return imb13;
	}
	else if (ix == 65)
	{
		IMB14 *imb14 = new IMB14;
		return imb14;
	}
	//66-73 UDF 1-9
	else if (ix == 66)
	{
		UDF1 *UD1 = new UDF1;
		return UD1;
	}
	else if (ix == 67)
	{
		UDF2 *UD2 = new UDF2;
		return UD2;
	}
	else if (ix == 68)
	{
		UDF3 *UD3 = new UDF3;
		return UD3;
	}
	else if (ix == 69)
	{
		UDF4 *UD4 = new UDF4;
		return UD4;
	}
	else if (ix == 70)
	{
		UDF5 *UD5 = new UDF5;
		return UD5;
	}
	else if (ix == 71)
	{
		UDF6 *UD6 = new UDF6;
		return UD6;
	}
	else if (ix == 72)
	{
		UDF8 *UD8 = new UDF8;
		return UD8;
	}
	else if (ix == 73)
	{
		UDF9 *UD9 = new UDF9;
		return UD9;
	}
	//74-83 JY 1-10
	else if (ix == 74)
	{
		JY1 *J1 = new JY1;
		return J1;
	}
	else if (ix == 75)
	{
		JY2 *J2 = new JY2;
		return J2;
	}
	else if (ix == 76)
	{
		JY3 *J3 = new JY3;
		return J3;
	}
	else if (ix == 77)
	{
		//not implemented
		//JY4 *J4 = new JY4;
		//return J4;
		return NULL;
	}
	else if (ix == 78)
	{
		JY5 *J5 = new JY5;
		return J5;
	}
	else if (ix == 79)
	{
		JY6 *J6 = new JY6;
		return J6;
	}
	else if (ix == 80)
	{
		JY7 *J7 = new JY7;
		return J7;
	}
	else if (ix == 81)
	{
		JY8 *J8 = new JY8;
		return J8;
	}
	else if (ix == 82)
	{
		//not implemented
		//JY9 *J9 = new JY9;
		//return J9;
		return NULL;
	}
	else if (ix == 83)
	{
		//not implemented
		//JY10 *J10 = new JY10;
		//return J10;
		return NULL;
	}
	//84-88 FDA 1-5
	else if (ix == 84)
	{
		FDA1 *F1 = new FDA1;
		return F1;
	}
	else if (ix == 85)
	{
		FDA2 *F2 = new FDA2;
		return F2;
	}
	else if (ix == 86)
	{
		FDA3 *F3 = new FDA3;
		return F3;
	}
	else if (ix == 87)
	{
		//not implemented
		//FDA4 *F4 = new FDA4;
		//return F4;
		return NULL;
	}
	else if (ix == 88)
	{
		//not implemented
		//FDA5 *F5 = new FDA5;
		//return F5;
		return NULL;
	}
	//89-98 CDF1-10
	else if (ix == 89)
	{
		CDF1 *CD1 = new CDF1;
		return CD1;
	}
	else if (ix == 90)
	{
		CDF2 *CD2 = new CDF2;
		return CD2;
	}
	else if (ix == 91)
	{
		CDF3 *CD3 = new CDF3;
		return CD3;
	}
	else if (ix == 92)
	{
		CDF4 *CD4 = new CDF4;
		return CD4;
	}
	else if (ix == 93)
	{
		CDF5 *CD5 = new CDF5;
		return CD5;
	}
	else if (ix == 94)
	{
		CDF6 *CD6 = new CDF6;
		return CD6;
	}
	else if (ix == 95)
	{
		CDF7 *CD7 = new CDF7;
		return CD7;
	}
	else if (ix == 96)
	{
		CDF8 *CD8 = new CDF8;
		return CD8;
	}
	else if (ix == 97)
	{
		CDF9 *CD9 = new CDF9;
		return CD9;
	}
	else if (ix == 98)
	{
		CDF10 *CD10 = new CDF10;
		return CD10;
	}
	else if (ix == 201)
	{
		GEN *gen = new GEN;
		return gen;
	}
	else if (ix == 202)
	{
		GEN_pat *gen2 = new GEN_pat;
		return gen2;
	}
	else
	{
		std::cout << "*****************************\n"
			<< "Function #" << ix << " not initialised\n"
			<< "Program will terminate"
			<< "\n*****************************\n";
		system("pause");
		abort();
	}

}

/*
*Function deactivation - for specific functions*
*/
std::vector<bool> Func_Deactivate()
{
	std::vector<bool> output(n_func_e - n_func_b + 1, true);	//Output vector of active functions, all are active

	//check which functions have to be skipped
	bool skip_two = false, skip_three = false, skip_mobj = false;
	for (int i = 0; i < func_skip.size(); i++)
	{
		if (func_skip[i] == "Two_obj")
			skip_two = true;
		else if (func_skip[i] == "Three_obj")
			skip_three = true;
		else if (func_skip[i] == "Many_obj")
			skip_mobj = true;
	}
	
	for (int i_func = n_func_b; i_func <= n_func_e; i_func++)
	{
		//Deactivate  functions that are not implemented
		if (Func_Type(i_func) == NULL)
		{
			output[i_func - n_func_b] = false;
			continue;
		}
		//Deactivate the functions based on func_skip - defined in const_h
		if (skip_two)
		{
			if (Func_Type(i_func)[0].Objs() == 2)
			{
				output[i_func - n_func_b] = false;
				continue;
			}
		}
		if (skip_three)
		{
			if (Func_Type(i_func)[0].Objs() == 3)
			{
				output[i_func - n_func_b] = false;
				continue;
			}
		}
		if (skip_mobj)
		{
			if (Func_Type(i_func)[0].Objs() > 3)
			{
				output[i_func - n_func_b] = false;
				continue;
			}
		}
	}

	return output;
}


mutation<short> * Mut_Type(short ix, std::string Mode)
{
	if (Mode == "MTS")
	{
		mutation<short> *MTS_M = new mutation<short>(0,"MTS");
		return MTS_M;
	}
	if (ENCODING == "Real")
	{
		if (ix == 1)
		{

			mutation_1<short> *M1 = new mutation_1<short>;
			return M1;
		}
		if (ix == 2)
		{
			mutation_1b<short> *M1 = new mutation_1b<short>;
			return M1;
		}
		if (ix == 3)
		{
			mutation_2<short> *M2 = new mutation_2<short>;
			return M2;
		}
		if (ix == 4)
		{
			mutation_3<short> * M3 = new mutation_3<short>;
			return M3;
		}
		else
		{
			std::cout << "*****************************\n"
				<< "Mutation #" << ix << " not initialised\n"
				<< "Program will terminate"
				<< "\n*****************************\n";
			system("pause");
			abort();
		}
	}
	else if (ENCODING == "Binary" || ENCODING == "Gray")
	{
		if (ix == 1)
		{
			mutation_B1<short> *M1 = new mutation_B1<short>;
			return M1;
		}
		else
		{
			std::cout << "*****************************\n"
				<< "Mutation #" << ix << " not initialised\n"
				<< "Program will terminate"
				<< "\n*****************************\n";
			system("pause");
			abort();
		}
	}
	else
	{
		std::cout << "*****************************\n"
			<< "Mutation #" << ix << " not initialised\n"
			<< "Program will terminate"
			<< "\n*****************************\n";
		system("pause");
		abort();
	}
}

selection<individual> * Select_Type(short ix)
{
	if (ix == 1)
	{
		roulette_wheel<individual> *S1 = new roulette_wheel<individual>;
		return S1;


	}
	else
	{
		std::cout << "*****************************\n"
			<< "Selection #" << ix << " not initialised"
			<< "Program will terminate"
			<< "\n*****************************\n";
		system("pause");
		abort();
	}
}
crossover<individual> * Cross_Type(short ix, std::string MODE)
{
	if (ENCODING == "Real")
	{
		if (ix == 1)
		{
			if (MODE == "Normal")
			{
				crossover_1<individual> *C1 = new crossover_1<individual>;
				return C1;
			}

			else if (MODE == "UNSGAIII" || MODE == "MOEAD" || MODE == "MOEADMSF" || MODE == "MOEADPSF" || MODE == "MOEADM2M" || MODE == "PAES" || MODE == "DMOEADD" || MODE == "BCE" || MODE == "HEIA" || MODE == "IBEA")
			{
				crossover_1b<individual> *C1 = new crossover_1b<individual>;
				return C1;
			}
			else if (MODE == "MTS")
			{
				crossover<individual> *MTS_C = new crossover<individual>(0, "MTS");
				return MTS_C;
			}
			else
			{
				std::cout << "*****************************\n"
					<< "Crossover #" << ix << " for MODE " << MODE << " not initialised\n"
					<< "Program will terminate"
					<< "\n*****************************\n";
				system("pause");
				abort();
			}
		}
		else
		{
			std::cout << "*****************************\n"
				<< "Crossover #" << ix << " not initialised\n"
				<< "Program will terminate"
				<< "\n*****************************\n";
			system("pause");
			abort();
		}
	}
	else if (ENCODING == "Binary" || ENCODING == "Gray")
	{
		if (MODE == "MTS")
		{
			crossover<individual> *MTS_C = new crossover<individual>(0, "MTS");
			return MTS_C;
		}
		else if (ix == 1)
		{
			crossover_B1<individual> *C1 = new crossover_B1<individual>;
			return C1;

		}
		else if (ix == 2)
		{
			crossover_B2<individual> *C2 = new crossover_B2<individual>;
			return C2;

		}
		else
		{
			std::cout << "*****************************\n"
				<< "Crossover #" << ix << " not initialised\n"
				<< "Program will terminate"
				<< "\n*****************************\n";
			system("pause");
			abort();
		}
	}
	else
	{
		std::cout << "*****************************\n"
			<< "Crossover #" << ix << " not initialised\n"
			<< "Program will terminate"
			<< "\n*****************************\n";
		system("pause");
		abort();
	}
}



std::string Date_Get()
{
	time_t t = time(0);
	struct tm * c_date = localtime(&t);
	int hour_t;
	std::string out, month, day, hour, min;
	if (c_date->tm_mon + 1 < 10)
		month = "0" + std::to_string(c_date->tm_mon + 1);
	else
		month = std::to_string(c_date->tm_mon + 1);

	if (c_date->tm_mday < 10)
		day = "0" + std::to_string(c_date->tm_mday);
	else
		day = std::to_string(c_date->tm_mday);

	hour_t = c_date->tm_hour;

	if (hour_t < 10)
		hour = "0" + std::to_string(hour_t);
	else
		hour = std::to_string(hour_t);

	if (c_date->tm_min < 10)
		min = "0" + std::to_string(c_date->tm_min);
	else
		min = std::to_string(c_date->tm_min);

	out = std::to_string(c_date->tm_year + 1900);
	out += "_" + month;
	out += "_" + day;
	out += "_" + hour;
	out += "_" + min + " ";
	return out;
}
std::string Name_Get(char * name, std::string & date, int indexr, float t)
{
	std::string name_temp;
	if (indexr < 0)
		name_temp = date + " " + name;
	else
		name_temp = date + " #" + std::to_string(indexr) + "_" + String_Prec(t,1) + "_" + name;
	return name_temp;
	/*const char * out = name_temp.c_str();
	return out;*/
}
std::string Name_Get(char * name, std::string & date, int indexr)
{
	std::string name_temp;
	if (indexr < 0)
		name_temp = date + " " + name;
	else
		name_temp = date + " #" + std::to_string(indexr) + "_" + name;
	return name_temp;
	/*const char * out = name_temp.c_str();
	return out;*/
}
std::string Name_Get(std::string name, int n)
{
	std::string name_temp = name;
	name_temp += "_" + std::to_string(n) + ".x1";
	return name_temp;
	/*const char * out1 = temp.c_str();
	return out1;*/
}

std::string Name_Get(std::string name, int n, double t)
{
	std::string name_temp = name;
	name_temp += "_" + std::to_string(n) + "_" + String_Prec(t, 1) + ".x1";
	return name_temp;
}


std::string F_Name_Get(char * name_front, std::string & date, char * name_back)
{
	std::string name_temp;
	name_temp = name_front + date;
	name_temp += name_back;
	return name_temp;
}

void Print(FILE * Pipe, int ir, int index_r_max, std::string date, double t, bool IGD_on, bool HV_on)
{
	if (SKIP_GRAPHS == true)
	{
		if ((ir - 1) % n_runs > 4)
			return;
	}


	std::string print_temp;
	print_temp = "set style line 1 lt rgb \"red\" \n ";	//Style line for points
	print_temp += "set style line 2 lt rgb \"blue\" \n ";	//Style line for PF
	print_temp += "set style line 3 lc rgb \"green\" pt 6 ps 0.25 \n ";	//Style line for real PF
	print_temp += "set style line 4 lc rgb \"green\" pt 7 ps 1 \n ";	//Style line for real PF legend
	fprintf(Pipe, print_temp.c_str());
	print_temp.clear();

	if (CONTOUR_PLOT == true && FITNESS_ALL == true)
		Contour_Plot(ir,t);
	// CONTOUR_PLOT
	//FILE *plotPipe = _popen("C:\\gnuplot\\bin\\gnuplot.exe -persistent", "w");
	fprintf(Pipe, "set xlabel \"Objective 1\"\nset ylabel \"Objective 2\"\n");
	print_temp = "set term wxt 0 size ";
	print_temp += std::to_string(frame_w);
	print_temp += ",";
	print_temp += std::to_string(frame_h);
	print_temp += "\n";
	fprintf(Pipe, print_temp.c_str());
	print_temp.clear();

	/* Plotting all results - Graph_all.jpg */
	if (FITNESS_ALL == true)
	{
		print_temp = "set title \"Individuals over generations\"\n";
		fprintf(Pipe, print_temp.c_str());
		print_temp.clear();
		/*if (index_r_max == 1)
		{
			print_temp = "set term wxt 0 size ";
			print_temp += std::to_string(frame_w);
			print_temp += ",";
			print_temp += std::to_string(frame_h);
			print_temp += "\nplot \"Temp\\\\";
			print_temp += Name_Get("graph", ir);
			print_temp += "\" using 1:2 with points ls 1";
			print_temp += " title \"All individuals\"";
			//Print PF if exisiting
			if (ONE_OBJ_OVERRIDE != true)
			{
				print_temp += ", ";
				print_temp += F_Name_Get("\"Results\\\\", date, "\\\\");
				print_temp += Name_Get("PF.csv", date, ir);
				print_temp += "\" using 1:2 with points ls 2";
				print_temp += " title \"Achieved PF\"";
			}
			print_temp += ", \"Temp\\\\";
			print_temp += Name_Get(real_PF_out, ir, t);
			print_temp += "\" using 1:2 with points ls 3 notitle";
			print_temp += ", 1/0 with points ls 4 title \"Real PF\"";
			print_temp += "\n";
			fprintf(Pipe, print_temp.c_str());
			print_temp.clear();
		}*/
		print_temp = "set term png size ";
		print_temp += std::to_string(frame_w);
		print_temp += ",";
		print_temp += std::to_string(frame_h);
		print_temp += "\nset output ";
		print_temp += F_Name_Get("\"Results\\\\", date, "\\\\");
		print_temp += Name_Get("Graph_all.jpg", std::string(), ir, t);
		print_temp += "\nplot \"Temp\\\\";
		print_temp += Name_Get("graph", ir, t);
		print_temp += "\" using 1:2 with points ls 1";
		print_temp += " title \"All individuals\"";
		//Print PF if exisiting
		if (ONE_OBJ_OVERRIDE != true)
		{
			print_temp += ", ";
			print_temp += F_Name_Get("\"Results\\\\", date, "\\\\");
			print_temp += Name_Get("PF.csv", std::string(), ir, t);
			print_temp += "\" using 1:2 with points ls 2";
			print_temp += " title \"Achieved PF\"";
		}
		print_temp += ", \"Temp\\\\";
		print_temp += Name_Get(real_PF_out, ir, t);
		print_temp += "\" using 1:2 with points ls 3 notitle";
		print_temp += ", 1/0 with points ls 4 title \"Real PF\"";
		print_temp += "\n";
		fprintf(Pipe, print_temp.c_str());
		print_temp.clear();
	} // FITNESS_ALL


	  /* Plotting Contour Plot - Contour_Plot.jpg */
	if (CONTOUR_PLOT == true && FITNESS_ALL == true)
	{
		//fprintf(Pipe, "set term wxt 2 \nset palette color positive\nset pm3d map \nsplot \"Graph2.x1\" \n");
		print_temp = "set zr[0:";
		print_temp += std::to_string(Get_Max_CP());
		print_temp += "] \n";
		fprintf(Pipe, print_temp.c_str());
		print_temp.clear();
		fprintf(Pipe, "set term wxt 2 \nset palette define ( 0 \"black\", 1 \"blue\", 6 \"red\", 11 \"orange\", 14 \"yellow\", 16 \"green\") \nset pm3d map \n");
		/*if (index_r_max == 1)
		{
			print_temp = "splot \"Temp\\\\";
			print_temp += Name_Get("graph2", ir);
			print_temp += "\" \n";
			fprintf(Pipe, print_temp.c_str());
			print_temp.clear();
		}*/
		print_temp = "set term png font arial 30 size 3840,2160 \nset output ";
		print_temp += F_Name_Get("\"Results\\\\", date, "\\\\");
		print_temp += Name_Get("Contour_Plot.jpg", std::string(), ir,t);
		print_temp += "\nsplot \"Temp\\\\";
		print_temp += Name_Get("CP", ir, t);
		print_temp += "\" \n";
		print_temp += "set term wxt 2 \nunset pm3d \n";
		fprintf(Pipe, print_temp.c_str());
		print_temp.clear();
	} // CONTOUR_PLOT

	  /* Plotting Pareto Front only - PF.jpg */
	if (ONE_OBJ_OVERRIDE != true)
	{
		print_temp = "set title \"Achieved Pareto Front\"\n";
		fprintf(Pipe, print_temp.c_str());
		print_temp.clear();
		/*if (index_r_max == 1)
		{
			print_temp = "set term wxt 1 \n set xr[*:*] \nset yr[*:*] \nplot ";
			print_temp += F_Name_Get("\"Results\\\\", date, "\\\\");
			print_temp += Name_Get("PF.csv", date, ir);
			print_temp += "\" using 1:2 with points ls 2";
			print_temp += " title \"Achieved PF\"";
			print_temp += ", \"Temp\\\\";
			print_temp += Name_Get(real_PF_out, ir, t);
			print_temp += "\" using 1:2 with points ls 3 notitle";
			print_temp += ", 1/0 with points ls 4 title \"Real PF\"";
			print_temp += "\n";
			fprintf(Pipe, print_temp.c_str());
			print_temp.clear();
		}*/
		print_temp = "set term png size ";
		print_temp += std::to_string(frame_w);
		print_temp += ",";
		print_temp += std::to_string(frame_h);
		print_temp += "\nset output ";
		print_temp += F_Name_Get("\"Results\\\\", date, "\\\\");
		print_temp += Name_Get("PF.jpg", std::string(), ir, t);
		print_temp += "\nset xr[*:*] \nset yr[*:*] \nplot ";
		print_temp += F_Name_Get("\"Results\\\\", date, "\\\\");
		print_temp += Name_Get("PF.csv", std::string(), ir, t);
		print_temp += "\" using 1:2 with points ls 2";
		print_temp += " title \"Achieved PF\"";
		print_temp += ", \"Temp\\\\";
		print_temp += Name_Get(real_PF_out, ir, t);
		print_temp += "\" using 1:2 with points ls 3 notitle";
		print_temp += ", 1/0 with points ls 4 title \"Real PF\"";
		print_temp += "\n";
		fprintf(Pipe, print_temp.c_str());
		print_temp.clear();
	}
	else if (ONE_OBJ_OVERRIDE == true)
	{
		fprintf(Pipe, "set xlabel \"Generations\"\nset ylabel \"Fitness\"\n");
		/*if (index_r_max == 1)
		{
			print_temp = "set term wxt 1 \n set xr[*:*] \nset yr[*:*] \nplot ";
			print_temp += F_Name_Get("\"Results\\\\", date, "\\\\");
			print_temp += Name_Get("PF.csv", date, ir);
			print_temp += "\" using 1:2 with points ls 2\n";
			fprintf(Pipe, print_temp.c_str());
			print_temp.clear();
		}*/
		print_temp = "set term png size ";
		print_temp += std::to_string(frame_w);
		print_temp += ",";
		print_temp += std::to_string(frame_h);
		print_temp += "\nset output ";
		print_temp += F_Name_Get("\"Results\\\\", date, "\\\\");
		print_temp += Name_Get("PF.jpg", std::string(), ir);
		print_temp += "\nset xr[*:*] \nset yr[*:*] \nplot ";
		print_temp += F_Name_Get("\"Results\\\\", date, "\\\\");
		print_temp += Name_Get("PF.csv", std::string(), ir);
		print_temp += "\" using 1:2 with points ls 2\n";
		fprintf(Pipe, print_temp.c_str());
		print_temp.clear();
	}
	/*Plotting the IGD graph*/
	if (PERF_VAL_GEN == true)
	{
		if (IGD_on == true)
		{
			fprintf(Pipe, "set xlabel \"Generation\"\nset ylabel \"IGD\"\n");
			print_temp = "set title \"IGD over generations\"\n";
			fprintf(Pipe, print_temp.c_str());
			print_temp.clear();

			print_temp = "set term png size ";
			print_temp += std::to_string(frame_w);
			print_temp += ",";
			print_temp += std::to_string(frame_h);
			print_temp += "\nset output ";
			print_temp += F_Name_Get("\"Results\\\\", date, "\\\\");
			print_temp += Name_Get("IGD.jpg", std::string(), ir);
			print_temp += "\nset xr[*:*] \nset yr[*:*] \nplot ";
			print_temp += F_Name_Get("\"Results\\\\", date, "\\\\");
			print_temp += Name_Get("IGD.csv", std::string(), ir);
			print_temp += "\" using 1:2 with points ls 4";
			print_temp += " title \"IGD\"";
			print_temp += "\n";
			fprintf(Pipe, print_temp.c_str());
			print_temp.clear();
		}
		if (HV_on == true)
		{
			fprintf(Pipe, "set xlabel \"Generation\"\nset ylabel \"HV\"\n");
			print_temp = "set title \"HV over generations\"\n";
			fprintf(Pipe, print_temp.c_str());
			print_temp.clear();

			print_temp = "set term png size ";
			print_temp += std::to_string(frame_w);
			print_temp += ",";
			print_temp += std::to_string(frame_h);
			print_temp += "\nset output ";
			print_temp += F_Name_Get("\"Results\\\\", date, "\\\\");
			print_temp += Name_Get("HV.jpg", std::string(), ir);
			print_temp += "\nset xr[*:*] \nset yr[*:*] \nplot ";
			print_temp += F_Name_Get("\"Results\\\\", date, "\\\\");
			print_temp += Name_Get("HV.csv", std::string(), ir);
			print_temp += "\" using 1:2 with points ls 4";
			print_temp += " title \"HV\"";
			print_temp += "\n";
			fprintf(Pipe, print_temp.c_str());
			print_temp.clear();
		}
	}

	fflush(Pipe);
}

void Directory_Create(std::string &date)
{
	if (DEBUG != true)
	{
		//Create the automatic folder name - if desired
		if (AUTO_FOLDER_COMMENT == true)
		{
			//Add the MLSt part
			//Check if MLSGA occur
			bool MLSGA_occ = false;					//indicator to MLSGA occurance
			if (MLSGA_Hybrid == true)
			{
				MLSGA_occ = true;
			}			
			//Add the MLSt part only when MLSGA occur
			if (MLSGA_occ == true)
			{

				date += " MLS";
				//If the only one value - add only one value
				if (MLSGA_n_MLS_b == MLSGA_n_MLS_e)
					date += std::to_string(MLSGA_n_MLS_b);
				//Else add the range
				else
					date += std::to_string(MLSGA_n_MLS_b) + "-" + std::to_string(MLSGA_n_MLS_e);
			}
			//Add the selected 
			std::string temp_date;
			for (int i = 0; i < GA_mode.size(); i++)
			{
				temp_date += " ";
				for (int iCoevol = 0; iCoevol < GA_mode[i].size(); iCoevol++)
				{
					temp_date += GA_mode[i][iCoevol];
					if (iCoevol != GA_mode[i].size() - 1)
						temp_date += "v";
				}
			}
			if (temp_date.size() >= 20)
				temp_date = " Many";

			date += temp_date;


			//Add population size
			//If the only one value - add only one value
			if (Pop_size_b == Pop_size_e)
				date += " " + std::to_string(Pop_size_b);
			//Else add the range
			else
				date += " " + std::to_string(Pop_size_b) + "-" + std::to_string(Pop_size_e) + "s" + std::to_string(Pop_size_step);
			date += "pop";

			if (T_con == "nfes")
			{
				//Add number of iterations
				//Add k for each 10^3 multiplication
				//Copy the value
				int temp_iteration_n = max_iterations;		//Value of iterations to save - integers
				int k_n = 0;								//Multiplication of 1000
				int temp_iteration_n_m = 0;					//Modulo
				//Derivate every 1000
				while (temp_iteration_n >= 1000)
				{
					//calculate modulo
					temp_iteration_n_m = temp_iteration_n % 1000;
					temp_iteration_n /= 1000;
					k_n++;
				}
				//Add to the string
				//Save the interger part
				date += " " + std::to_string(temp_iteration_n);
				//Save the rest
				if (k_n > 0);
				{
					//Add the rest - if exist
					if (temp_iteration_n_m != 0)
						date += "." + std::to_string(temp_iteration_n_m);
					//Add the k
					for (int i = 0; i < k_n; i++)
						date += "k";
					date += "_";
				}
				//Add the title
				date += "iter";
			}
			else if (T_con == "ngen")
			{
				date += " " + std::to_string(max_generations) + "gen";
			}
			else if (T_con == "ntime")
			{
				date += " " + std::to_string(max_time) + "sec";
			}

			//Save the number of collectives
			//If the only one value - add only one value
			if (MLSGA_Hybrid == true)
			{
				if (MLSGA_n_col_b == MLSGA_n_col_e)
					date += " " + std::to_string(MLSGA_n_col_b);
				//Else add the range
				else
					date += " " + std::to_string(MLSGA_n_col_b) + "-" + std::to_string(MLSGA_n_col_e);
				date += "col";

				//Save the number of collective skip
				date += " " + std::to_string(MLSGA_col_elim_limit) + "elimit";
			}			

			//Save the number of runs - if greater than 1
			if (n_runs > 1)
				date += " " + std::to_string(n_runs) + "run";
		}
		//Create the directory for results - if not existing
		CreateDirectoryA("Results/", NULL);
		//Create the directory for temp - if not existing
		CreateDirectoryA("Temp/", NULL);
		//Create the directory for current results
		std::string temp1 = F_Name_Get("Results/", date, "/");
		const char * temp = temp1.c_str();
		if (!CreateDirectoryA(temp, NULL))
		{
			std::cout << "**********************************\n"
				<< "Folder exists or cannot be created\nProgram will terminate"
				<< "\n**********************************\n";
			system("pause");
			abort();
		}
	}
}

extern time_t mut_t;										//Time of Mutation
extern time_t pop_t;										//Time of Population initialization,
extern time_t SVM_t;										//Time of SVM
extern time_t run_t;										//Time of run
extern time_t col_t;										//Time of collective generation
extern time_t elit_t;										//Time of Elitism
extern time_t cross_t;										//Time of Crossover
extern time_t selec_t;										//Time of Selection
extern time_t PF_t;											//Time of PF creation
extern time_t save_t;										//Time of results saving
extern float run_time;
extern float GA_time;
extern std::ofstream file;
SimpleXlsx::CWorksheet sheet1, sheet2, sheet3, sheet4;						//For saving data to excel
SimpleXlsx::CWorkbook book1;												//Book - for all data	

void Time_Save(int ir, double IGD_val, double HV_val, double min_fitness, int end_generation,int iteration, std::string MODE)
{
	/*************Time calculation and saving************/

	
	float mut_time = (mut_t * 1000 / CLOCKS_PER_SEC);
	mut_time /= 1000;
	float cross_time = (cross_t * 1000 / CLOCKS_PER_SEC);
	cross_time /= 1000;
	float selec_time = (selec_t * 1000 / CLOCKS_PER_SEC);
	selec_time /= 1000;
	float pop_time = (pop_t * 1000 / CLOCKS_PER_SEC);
	pop_time /= 1000;
	float pf_time = (PF_t * 1000 / CLOCKS_PER_SEC);
	pf_time /= 1000;
	float svm_time = (SVM_t * 1000 / CLOCKS_PER_SEC);
	svm_time /= 1000;
	float col_time = (col_t * 1000 / CLOCKS_PER_SEC);
	col_time /= 1000;
	float elit_time = (elit_t * 1000 / CLOCKS_PER_SEC);
	elit_time /= 1000;
	float save_time = (save_t * 1000 / CLOCKS_PER_SEC);
	save_time /= 1000;
	float sum_time = run_time - mut_time - pop_time - pf_time - svm_time - col_time - elit_time;
	if (MODE == "normal")
		GA_time = (pop_time + mut_time + col_time + elit_time - save_time);
	else
		GA_time = (pop_time + mut_time + col_time + elit_time + cross_time + selec_time);
	file << "Time of " << ir << "#: " << run_time << std::endl;
	file << "Time of: mutation " << mut_time << std::endl
		<< "crossover: " << cross_time << std::endl
		<< "selection: " << selec_time << std::endl
		<< "Collective operations: " << col_time << std::endl
		<< "population creation: " << (pop_time - cross_time - selec_time) << std::endl
		<< "elitism: " << elit_time << std::endl
		<< "Sum of GA time: " << GA_time << std::endl
		<< "PF: " << pf_time << std::endl
		<< "SVM/Clustering: " << svm_time << std::endl
		<< "Other: " << sum_time << std::endl
		<< "IGD: " << IGD_val << std::endl
		<< "HV: " << HV_val << std::endl
		<< "Saving: " << save_time << std::endl << std::endl;
	std::cout << "\n\nTime of " << ir << "#: " << run_time << std::endl;
	std::cout << "Time of: mutation " << mut_time << std::endl
		<< "crossover: " << cross_time << std::endl
		<< "selection: " << selec_time << std::endl
		<< "Collective operations: " << col_time << std::endl
		<< "population creation: " << (pop_time - cross_time - selec_time) << std::endl
		<< "elitism: " << elit_time << std::endl
		<< "Sum of GA time: " << GA_time << std::endl
		<< "PF: " << pf_time << std::endl
		<< "SVM/Clustering: " << svm_time << std::endl
		<< "Other: " << sum_time << std::endl
		<< "IGD: " << IGD_val << std::endl
		<< "HV: " << HV_val << std::endl
		<< "Saving: " << save_time << std::endl << std::endl;;

	sheet3.BeginRow();
	sheet3.AddCell("#" + std::to_string((int)((ir - 1) / n_runs) + 1));
	sheet3.AddCell(std::to_string((int)((ir - 1) % n_runs + 1)));
	sheet3.AddCell(std::to_string(run_time));
	sheet3.AddCell(std::to_string(cross_time));
	sheet3.AddCell(std::to_string(selec_time));
	sheet3.AddCell(std::to_string(mut_time));
	sheet3.AddCell(std::to_string(col_time));
	sheet3.AddCell(std::to_string(elit_time));
	sheet3.AddCell(std::to_string((pop_time - cross_time - selec_time)));
	sheet3.AddCell(std::to_string(GA_time));
	sheet3.AddCell(std::to_string(pf_time));
	sheet3.AddCell(std::to_string(svm_time));
	sheet3.AddCell(std::to_string(save_time));
	sheet3.AddCell(std::to_string(sum_time));
	if (ONE_OBJ_OVERRIDE != true)
	{
		sheet3.AddCell(std::to_string(IGD_val));
		sheet3.AddCell(std::to_string(HV_val));
	}
	sheet3.AddCell(std::to_string(end_generation));
	sheet3.AddCell(std::to_string(iteration));
	if (ONE_OBJ_OVERRIDE == true)
	{
		sheet3.AddCell(std::to_string(min_fitness));

	}
	sheet3.EndRow();


}


/*
*Excel output generation - First row*
@param ir - index of the run
@param MODE - mode of run - for PAES only
@param dyn - if function is dynamic
*/
void Excel_F_Row(int indexr, std::vector<std::string> MODE, bool dyn, bool IGD_on, bool HV_on)
{

	std::vector<SimpleXlsx::CellDataStr> excel_begin;				//data storage for the first row for the excel sheet - data tab for each run

																	//push back titles to the data storage vector
	excel_begin.push_back("Index");
	excel_begin.push_back("Run");
	excel_begin.push_back("Modes");
	excel_begin.push_back("RI Mode");
	excel_begin.push_back("Termination");
	excel_begin.push_back("T value");
	excel_begin.push_back("Selection");
	excel_begin.push_back("Crossover");
	excel_begin.push_back("Crossover rate");
	excel_begin.push_back("Mutation");
	excel_begin.push_back("Mutation rate");
	excel_begin.push_back("Function");
	if (TGM == true)
	{
		excel_begin.push_back("Maternal eff");
		excel_begin.push_back("Gmaternal eff");
	}
	excel_begin.push_back("Population");
	excel_begin.push_back("Generations");
	if (MLSGA_Hybrid == true)
	{
		excel_begin.push_back("Collectives");
		excel_begin.push_back("Collectives E");
		excel_begin.push_back("Fit for indi");
		excel_begin.push_back("Fit for col");
		excel_begin.push_back("Elitism rate (%)");
		excel_begin.push_back("MLS");
		excel_begin.push_back("Col elim limit");
	}
	excel_begin.push_back("Number of constrains");
	excel_begin.push_back("niche multi");
	excel_begin.push_back("limit multi");
	excel_begin.push_back("mating chance");
	excel_begin.push_back("One objective override");
	excel_begin.push_back("PF_Refine");
	excel_begin.push_back("Real random");
	excel_begin.push_back("File input");
	excel_begin.push_back("Sobol");
	if (indexr == 1)
	{
		std::vector<SimpleXlsx::CellDataStr> excel_begin2;		//data storage for the first row for the excel sheet - Index tab
		std::vector<SimpleXlsx::CellDataStr> excel_begin3;		//data storage for the first row for the excel sheet - Time and IGD values tab
		std::vector<SimpleXlsx::CellDataStr> excel_begin4;		//data storage for the first row for the excel sheet - GA data tab

																//Create the first sheet of the excel file - GA parameters
		sheet2 = book1.AddSheet("INDEX");

		//copy the data storage vector
		excel_begin2 = excel_begin;

		//Add row of data to the sheet
		sheet2.AddRow(excel_begin2);

		//push back titles to the data storage vector
		excel_begin3.push_back("Time of (Index):");
		excel_begin3.push_back("Run (Index)");
		excel_begin3.push_back("Run");
		excel_begin3.push_back("Crossover");
		excel_begin3.push_back("Selection");
		excel_begin3.push_back("Mutation");
		excel_begin3.push_back("Collective operations");
		excel_begin3.push_back("Elitism");
		excel_begin3.push_back("Population creation");
		excel_begin3.push_back("Sum of GA time(Only)");
		excel_begin3.push_back("PF");
		excel_begin3.push_back("SVM/Clustering");
		excel_begin3.push_back("Saving");
		excel_begin3.push_back("Other (Tools)");
		if (ONE_OBJ_OVERRIDE != true)					//IGD only for multi objective
		{
			excel_begin3.push_back("IGD");
			excel_begin3.push_back("HV");
		}
		excel_begin3.push_back("End generation");
		excel_begin3.push_back("End iteration");
		if (ONE_OBJ_OVERRIDE == true)					//Fitness only for one objective
		{
			excel_begin3.push_back("Fitness");

		}
		//Create the second sheet of the excel file - time values for each run
		sheet3 = book1.AddSheet("Time_IGD");

		//Add row of data to the sheet
		sheet3.AddRow(excel_begin3);

		//push back titles to the data storage vector
		excel_begin4.push_back("Index");
		excel_begin4.push_back("Min time (run)");
		excel_begin4.push_back("Max time (run)");
		excel_begin4.push_back("Average time (run)");
		excel_begin4.push_back("Std deviation (time run)");
		excel_begin4.push_back("Min time (GA only)");
		excel_begin4.push_back("Max time (GA only)");
		excel_begin4.push_back("Average time (GA only)");
		excel_begin4.push_back("Std deviation (GA time only)");
		if (ONE_OBJ_OVERRIDE != true)				//IGD only for multi objective
		{
			if (IGD_on == true)
			{
				excel_begin4.push_back("Min IGD");
				excel_begin4.push_back("Max IGD");
				excel_begin4.push_back("Average IGD");
				excel_begin4.push_back("Std deviation (IGD)");
				if (dyn == true)
				{
					excel_begin4.push_back("Min IGD2");
					excel_begin4.push_back("Max IGD2");
					excel_begin4.push_back("Average IGD2");
					excel_begin4.push_back("Std deviation (IGD2)");
				}

			}
			if (HV_on == true)
			{
				excel_begin4.push_back("Min HV");
				excel_begin4.push_back("Max HV");
				excel_begin4.push_back("Average HV");
				excel_begin4.push_back("Std deviation (HV)");
				if (dyn == true)
				{
					excel_begin4.push_back("Min HV2");
					excel_begin4.push_back("Max HV2");
					excel_begin4.push_back("Average HV2");
					excel_begin4.push_back("Std deviation (HV2)");
				}
			}
			
		}
		excel_begin4.push_back("Min generation");
		excel_begin4.push_back("Max generation");
		excel_begin4.push_back("Average generation");
		excel_begin4.push_back("Std deviation (generations)");
		excel_begin4.push_back("Min iteration");
		excel_begin4.push_back("Max iteration");
		excel_begin4.push_back("Average iteration");
		excel_begin4.push_back("Std deviation (iterations)");
		if (ONE_OBJ_OVERRIDE == true)				//Fitness only for one objective
		{
			excel_begin4.push_back("Min fitness");
			excel_begin4.push_back("Max fitness");
			excel_begin4.push_back("Average fitness");
			excel_begin4.push_back("Std deviation (fitness)");
			excel_begin4.push_back("Succes rate [%]");
		}
		excel_begin4.push_back("Constrain violation count");
		

		//Create the third sheet of the excel file - GA data (min,max,avg,std deviation)
		sheet4 = book1.AddSheet("GA_data");

		//Add row of data to the sheet
		sheet4.AddRow(excel_begin4);
	}

	if (EXCEL_EXCEPTION == false)

	{
		//Create the next sheet of the excel file - fitness data for each run
		sheet1 = book1.AddSheet(std::to_string(indexr));

		//Add row of data to the sheet
		sheet1.AddRow(excel_begin);
	}



}

void Excel_S_Row(function & fcode, GA_parameters & gapara, mutation<short> & mcode, selection<individual> & scode, crossover<individual> & ccode, int irun, int ir, std::string rimode, std::vector<std::string> MODE, short MLS)
{
	std::string Modes;
	bool Normal_on = false, MOEAD_on = false, NSGA_on = false;
	for (int iMode = 0; iMode < MODE.size(); iMode++)
	{
		Modes += MODE[iMode] + " ";
		if (MODE[iMode] == "Normal")
			Normal_on = true;
		else if (MODE[iMode] == "UNSGAIII")
			NSGA_on = true;
		else if (MODE[iMode] == "MOEAD" || MODE[iMode] == "MOEADMSF" || MODE[iMode] == "MOEADPSF" || MODE[iMode] == "MOEADM2M" || MODE[iMode] == "PAES")
			MOEAD_on = true;
	}

	/**************Excel second row*******************/
	sheet1.BeginRow();
	sheet1.AddCell("#" + std::to_string((int)((ir - 1) / n_runs) + 1));
	sheet1.AddCell(std::to_string(irun + 1));
	sheet1.AddCell(Modes);
	sheet1.AddCell(rimode);
	if (T_con == "nfes")
	{
		sheet1.AddCell("Iterations");
		sheet1.AddCell(std::to_string(max_iterations));
	}
	else if (T_con == "ngen")
	{
		sheet1.AddCell("Generations");
		sheet1.AddCell(std::to_string(max_generations));
	}
	else if (T_con == "ntime")
	{
		sheet1.AddCell("Time");
		sheet1.AddCell(std::to_string(max_time));
	}
	if (Normal_on == true)
		sheet1.AddCell(scode.Name_Show());
	else
		sheet1.AddCell(Modes);

	sheet1.AddCell(ccode.Name_Show());		
	
	sheet1.AddCell(std::to_string(gapara.Cross_Prob()));
	sheet1.AddCell(mcode.Name_Show());
	sheet1.AddCell(std::to_string(gapara.Mut_Prob()));
	sheet1.AddCell(fcode.Name_Show());
	if (TGM == true)
	{
		sheet1.AddCell(std::to_string(TGM_vect[0]));
		sheet1.AddCell(std::to_string(TGM_vect[1]));
	}
	sheet1.AddCell(std::to_string(gapara.Pop_Size()));
	sheet1.AddCell(std::to_string(gapara.Max_Gen()));
	if (MLSGA_Hybrid == true)
	{
		sheet1.AddCell(std::to_string(n_col));
		sheet1.AddCell(std::to_string(MLSGA_n_col_elim));
		sheet1.AddCell("f" + std::to_string(fit_index_sel));
		sheet1.AddCell("f" + std::to_string(fit_index_col_sel));
		sheet1.AddCell(std::to_string(MLSGA_elite_size));
		sheet1.AddCell(std::to_string(MLS));
		sheet1.AddCell(std::to_string(MLSGA_col_elim_limit));
	}	
	sheet1.AddCell(std::to_string(fcode.Cons()));

	if (MOEAD_on == true)
	{
		sheet1.AddCell(std::to_string(MOEAD_niche_multi));
	}
	else
	{
		sheet1.AddCell("w/o");
	}
	if (MOEAD_on == true)
	{
		sheet1.AddCell(std::to_string(MOEAD_limit_multi));
		sheet1.AddCell(std::to_string(MOEAD_mating_chance));
	}
	else
	{
		sheet1.AddCell("w/o");
		sheet1.AddCell("w/o");
	}
	if (ONE_OBJ_OVERRIDE == false)
		sheet1.AddCell("false");
	else
		sheet1.AddCell("true");
	if (PF_REFINE == false)
		sheet1.AddCell("false");
	else
		sheet1.AddCell("true");
	if (REAL_RANDOM == false)
		sheet1.AddCell("false");
	else
		sheet1.AddCell("true");
	if (FILE_INPUT == false)
		sheet1.AddCell("false");
	else
		sheet1.AddCell("true");
	if (SOBOL == false)
		sheet1.AddCell("false");
	else
		sheet1.AddCell("true");
	sheet1.EndRow();

	sheet2.BeginRow();
	sheet2.AddCell("#" + std::to_string((int)((ir - 1) / n_runs) + 1));
	sheet2.AddCell(std::to_string(irun + 1));
	sheet2.AddCell(Modes);
	sheet2.AddCell(rimode);
	if (T_con == "nfes")
	{
		sheet2.AddCell("Iterations");
		sheet2.AddCell(std::to_string(max_iterations));
	}
	else if (T_con == "ngen")
	{
		sheet2.AddCell("Generations");
		sheet2.AddCell(std::to_string(max_generations));
	}
	else if (T_con == "ntime")
	{
		sheet2.AddCell("Time");
		sheet2.AddCell(std::to_string(max_time));
	}
	if (Normal_on == true)
		sheet2.AddCell(scode.Name_Show());
	else
		sheet2.AddCell(Modes);
	
	sheet2.AddCell(ccode.Name_Show());
	sheet2.AddCell(std::to_string(gapara.Cross_Prob()));

	sheet2.AddCell(mcode.Name_Show());
	sheet2.AddCell(std::to_string(gapara.Mut_Prob()));
	sheet2.AddCell(fcode.Name_Show());
	if (TGM == true)
	{
		sheet2.AddCell(std::to_string(TGM_vect[0]));
		sheet2.AddCell(std::to_string(TGM_vect[1]));
	}
	sheet2.AddCell(std::to_string(gapara.Pop_Size()));
	sheet2.AddCell(std::to_string(gapara.Max_Gen()));
	if (MLSGA_Hybrid == true)
	{
		sheet2.AddCell(std::to_string(n_col));
		sheet2.AddCell(std::to_string(MLSGA_n_col_elim));
		sheet2.AddCell("f" + std::to_string(fit_index_sel));
		sheet2.AddCell("f" + std::to_string(fit_index_col_sel));
		sheet2.AddCell(std::to_string(MLSGA_elite_size));
		sheet2.AddCell(std::to_string(MLS));
		sheet2.AddCell(std::to_string(MLSGA_col_elim_limit));
	}
	sheet2.AddCell(std::to_string(fcode.Cons()));
	if (MOEAD_on == true)
	{
		sheet2.AddCell(std::to_string(MOEAD_niche_multi));
	}
	else
	{
		sheet2.AddCell("w/o");
	}
	if (MOEAD_on == true)
	{
		sheet2.AddCell(std::to_string(MOEAD_limit_multi));
		sheet2.AddCell(std::to_string(MOEAD_mating_chance));
	}
	else
	{
		sheet2.AddCell("w/o");
		sheet2.AddCell("w/o");
	}
	if (ONE_OBJ_OVERRIDE == false)
		sheet2.AddCell("false");
	else
		sheet2.AddCell("true");
	if (PF_REFINE == false)
		sheet2.AddCell("false");
	else
		sheet2.AddCell("true");
	if (REAL_RANDOM == false)
		sheet2.AddCell("false");
	else
		sheet2.AddCell("true");
	if (FILE_INPUT == false)
		sheet2.AddCell("false");
	else
		sheet2.AddCell("true");
	if (SOBOL == false)
		sheet2.AddCell("false");
	else
		sheet2.AddCell("true");
	sheet2.EndRow();

}

void Excel_GA_data(int success_num, GA_data<float> & data, bool dyn, bool IGD_on, bool HV_on)
{
	static int index;							//index of the data
	index++;
	//begin new row in the sheet
	sheet4.BeginRow();
	//add data to the sheet
	sheet4.AddCell("#" + std::to_string(index));
	sheet4.AddCell(std::to_string(data.Show_Time_Struct().min));
	sheet4.AddCell(std::to_string(data.Show_Time_Struct().max));
	sheet4.AddCell(std::to_string(data.Show_Time_Struct().avg));
	sheet4.AddCell(std::to_string(data.Show_Time_Struct().std_deviation));
	sheet4.AddCell(std::to_string(data.Show_GA_Time_Struct().min));
	sheet4.AddCell(std::to_string(data.Show_GA_Time_Struct().max));
	sheet4.AddCell(std::to_string(data.Show_GA_Time_Struct().avg));
	sheet4.AddCell(std::to_string(data.Show_GA_Time_Struct().std_deviation));
	if (ONE_OBJ_OVERRIDE != true)
	{
		if (IGD_on == true)
		{
			sheet4.AddCell(std::to_string(data.Show_IGD_Struct().min));
			sheet4.AddCell(std::to_string(data.Show_IGD_Struct().max));
			sheet4.AddCell(std::to_string(data.Show_IGD_Struct().avg));
			sheet4.AddCell(std::to_string(data.Show_IGD_Struct().std_deviation));
			if (dyn == true)
			{
				sheet4.AddCell(std::to_string(data.Show_IGD2_Struct().min));
				sheet4.AddCell(std::to_string(data.Show_IGD2_Struct().max));
				sheet4.AddCell(std::to_string(data.Show_IGD2_Struct().avg));
				sheet4.AddCell(std::to_string(data.Show_IGD2_Struct().std_deviation));
			}
		}
		if (HV_on == true)
		{
			sheet4.AddCell(std::to_string(data.Show_HV_Struct().min));
			sheet4.AddCell(std::to_string(data.Show_HV_Struct().max));
			sheet4.AddCell(std::to_string(data.Show_HV_Struct().avg));
			sheet4.AddCell(std::to_string(data.Show_HV_Struct().std_deviation));
			if (dyn == true)
			{
				sheet4.AddCell(std::to_string(data.Show_HV2_Struct().min));
				sheet4.AddCell(std::to_string(data.Show_HV2_Struct().max));
				sheet4.AddCell(std::to_string(data.Show_HV2_Struct().avg));
				sheet4.AddCell(std::to_string(data.Show_HV2_Struct().std_deviation));
			}
		}
	}
	sheet4.AddCell(std::to_string(data.Show_Generation_Struct().min));
	sheet4.AddCell(std::to_string(data.Show_Generation_Struct().max));
	sheet4.AddCell(std::to_string(data.Show_Generation_Struct().avg));
	sheet4.AddCell(std::to_string(data.Show_Generation_Struct().std_deviation));
	sheet4.AddCell(std::to_string(data.Show_Iteration_Struct().min));
	sheet4.AddCell(std::to_string(data.Show_Iteration_Struct().max));
	sheet4.AddCell(std::to_string(data.Show_Iteration_Struct().avg));
	sheet4.AddCell(std::to_string(data.Show_Iteration_Struct().std_deviation));
	//if one objective optimisatiom add fitness
	if (ONE_OBJ_OVERRIDE == true)
	{
		sheet4.AddCell(std::to_string(data.Show_Fitness_Struct().min));
		sheet4.AddCell(std::to_string(data.Show_Fitness_Struct().max));
		sheet4.AddCell(std::to_string(data.Show_Fitness_Struct().avg));
		sheet4.AddCell(std::to_string(data.Show_Fitness_Struct().std_deviation));
		sheet4.AddCell(std::to_string((float)((float)success_num / (float)n_runs)*100.f));
	}
	sheet4.AddCell(std::to_string(cons_viol_count));

	//end sheet
	sheet4.EndRow();
}


/*Checking if everything is defined properly.Looking for errors in define.h and const.h*/
void Error_Check()
{
	if (max_iterations != 300000)
	{
		std::cout << "***************Warning*******************\n"
			<< "Wrong number of iterations chosen. Do not met CEC09 standards. Check Pop size and max gen. Const.h\nDo You want to continue anyway?"
			<< "\n**********************************\n";
		system("pause");
	}

	for (int i = 0; i < GA_mode.size(); i++)
	{
		if (GA_mode[i].size() > MLSGA_n_col_b || (GA_mode[i].size() > 1 && MLSGA_Hybrid != true))
		{
			{
				std::cout << "**********************************\n"
					<< "Wrong mode chosen. Number of different coevolutionary algorithms cannot be bigger than number of collectives. Check Const.h\nProgram will terminate"
					<< "\n**********************************\n";
				system("pause");
				abort();
			}
		}


		for (int iCoevol = 0; iCoevol < GA_mode[i].size(); iCoevol++)
		{
			std::string MODE = GA_mode[i][iCoevol];
			//Check if correct mode is selected
			if (MODE == "Normal")
			{
				/*if (MLSGA_Hybrid != true)
				{
					std::cout << "**********************************\n"
						<< "Wrong mode chosen. Check Const.h\nProgram will terminate"
						<< "\n**********************************\n";
					system("pause");
					abort();
				}*/


				//Check number of collectives
				if ((MLSGA_n_col_b != 4 && MLSGA_n_col_b != 6 && MLSGA_n_col_b != 8) && (MLSGA_n_col_e != 4 && MLSGA_n_col_e != 6 && MLSGA_n_col_e != 8) && GROUPING == "SVM")
				{
					std::cout << "**********************************\n"
						<< "Wrong number of collective chosen for SVM. Check Const.h\nProgram will terminate"
						<< "\n**********************************\n";
					system("pause");
					abort();
				}
				//Check iterations number
				/*if (max_iterations != 300000)
				{
					std::cout << "***************Warning*******************\n"
						<< "Wrong number of iterations chosen. Do not met CEC09 standards. Check Pop size and max gen. Const.h\nDo You want to continue anyway?"
						<< "\n**********************************\n";
					system("pause");
				}*/
			}
			else if (MODE == "MOEAD" || MODE == "MOEADMSF" || MODE == "MOEADPSF" || MODE == "MOEADM2M")
			{
				//Check number of objectives
				if (ONE_OBJ_OVERRIDE != false)
				{
					std::cout << "**********************************\n"
						<< "Wrong parameters chosen. Cannot use MOEAD for one objective Check Define.h\nProgram will terminate"
						<< "\n**********************************\n";
					system("pause");
					abort();
				}
				//Check number of collectives
				if ((MLSGA_n_col_b != 1 && MLSGA_n_col_e != 1) && ((MLSGA_n_col_b != 4 && MLSGA_n_col_b != 6 && MLSGA_n_col_b != 8) && (MLSGA_n_col_e != 4 && MLSGA_n_col_e != 6 && MLSGA_n_col_e != 8) && GROUPING == "SVM"))
				{
					std::cout << "**********************************\n"
						<< "Wrong number of collective chosen. Check Const.h\nProgram will terminate"
						<< "\n**********************************\n";
					system("pause");
					abort();
				}
				if (((Pop_size_b != 300 && Pop_size_b != 400 && Pop_size_b != 500 && Pop_size_b != 600 && Pop_size_b != 800 && Pop_size_b != 1000) || (Pop_size_e != 300 && Pop_size_e != 400 && Pop_size_e != 500 && Pop_size_e != 600 && Pop_size_e != 800 && Pop_size_e != 1000)) && MOEAD_weight_generator == "File")
				{
					std::cout << "**********************************\n"
						<< "Wrong Population size. MOEAD works only for 300, 400, 500, 600, 800 & 1000. Check Const.h\nProgram will terminate"
						<< "\n**********************************\n";
					system("pause");
					abort();
				}
				//Check iterations number
				/*if (max_iterations != 300000)
				{
					std::cout << "***************Warning*******************\n"
						<< "Wrong number of iterations chosen. Do not met CEC09 standards. Check Pop size and max gen in Const.h\nDo You want to continue anyway?"
						<< "\n**********************************\n";
					system("pause");
				}*/
			}
			else if (MODE == "UNSGAIII")
			{
				//Check number of objectives
				if (ONE_OBJ_OVERRIDE != false)
				{
					std::cout << "**********************************\n"
						<< "Wrong parameters chosen. Cannot use MOEAD for one objective Check Define.h\nProgram will terminate"
						<< "\n**********************************\n";
					system("pause");
					abort();
				}
				//Check number of collectives
				if ((MLSGA_n_col_b != 1 && MLSGA_n_col_e != 1) && ((MLSGA_n_col_b != 4 && MLSGA_n_col_b != 6 && MLSGA_n_col_b != 8) && (MLSGA_n_col_e != 4 && MLSGA_n_col_e != 6 && MLSGA_n_col_e != 8) && GROUPING == "SVM"))
				{
					std::cout << "**********************************\n"
						<< "Wrong number of collective chosen. Check Const.h\nProgram will terminate"
						<< "\n**********************************\n";
					system("pause");
					abort();
				}
				//Check iterations number
				/*if (max_iterations != 300000)
				{
					std::cout << "***************Warning*******************\n"
						<< "Wrong number of iterations chosen. Do not met CEC09 standards. Check Pop size and max gen. Const.h\nDo You want to continue anyway?"
						<< "\n**********************************\n";
					system("pause");
				}*/
			}
			else if (MODE == "MTS")
			{
				//Check number of objectives
				if (ONE_OBJ_OVERRIDE != false)
				{
					std::cout << "**********************************\n"
						<< "Wrong parameters chosen. Cannot use MOEAD for one objective Check Define.h\nProgram will terminate"
						<< "\n**********************************\n";
					system("pause");
					abort();
				}
				//Check iterations number
				/*if (max_iterations != 300000)
				{
					std::cout << "***************Warning*******************\n"
						<< "Wrong number of iterations chosen. Do not met CEC09 standards. Check Pop size and max gen. Const.h\nDo You want to continue anyway?"
						<< "\n**********************************\n";
					system("pause");
				}*/
			}
			else if (MODE == "DMOEADD")
			{
				//Check number of objectives
				if (ONE_OBJ_OVERRIDE != false)
				{
					std::cout << "**********************************\n"
						<< "Wrong parameters chosen. Cannot use MOEAD for one objective Check Define.h\nProgram will terminate"
						<< "\n**********************************\n";
					system("pause");
					abort();
				}
				//Check number of collectives
				if ((MLSGA_n_col_b != 1 && MLSGA_n_col_e != 1) && ((MLSGA_n_col_b != 4 && MLSGA_n_col_b != 6 && MLSGA_n_col_b != 8) && (MLSGA_n_col_e != 4 && MLSGA_n_col_e != 6 && MLSGA_n_col_e != 8) && GROUPING == "SVM"))
				{
					std::cout << "**********************************\n"
						<< "Wrong number of collective chosen. Check Const.h\nProgram will terminate"
						<< "\n**********************************\n";
					system("pause");
					abort();
				}
				//Check iterations number
				/*if (max_iterations != 300000)
				{
					std::cout << "***************Warning*******************\n"
						<< "Wrong number of iterations chosen. Do not met CEC09 standards. Check Pop size and max gen. Const.h\nDo You want to continue anyway?"
						<< "\n**********************************\n";
					system("pause");
				}*/
			}
			else if (MODE == "PAES")
			{
				
				//Check number of objectives
				if (ONE_OBJ_OVERRIDE != false)
				{
					std::cout << "**********************************\n"
						<< "Wrong parameters chosen. Cannot use MOEAD for one objective Check Define.h\nProgram will terminate"
						<< "\n**********************************\n";
					system("pause");
					abort();
				}
				//Check number of collectives
				if (n_col != 1 && n_col != 4 && n_col != 6 && n_col != 8)
				{
					std::cout << "**********************************\n"
						<< "Wrong number of collective chosen. Check Const.h\nProgram will terminate"
						<< "\n**********************************\n";
					system("pause");
					abort();
				}
				if (n_cross != 1)
				{
					std::cout << "**********************************\n"
						<< "Wrong number of crossover types. Cannot be greater than 1. MOEAD have internal crossover. Check Const.h\nProgram will terminate"
						<< "\n**********************************\n";
					system("pause");
					abort();
				}
				if ((pop_size != 300 && pop_size != 400 && pop_size != 500 && pop_size != 600 && pop_size != 800 && pop_size != 1000) && MOEAD_weight_generator == "File")
				{
					std::cout << "**********************************\n"
						<< "Wrong Population size. MOEAD works only for 300, 400, 500, 600, 800 & 1000. Check Const.h\nProgram will terminate"
						<< "\n**********************************\n";
					system("pause");
					abort();
				}
				//Check iterations number
				/*if (max_iterations != 300000)
				{
					std::cout << "***************Warning*******************\n"
						<< "Wrong number of iterations chosen. Do not met CEC09 standards. Check Pop size and max gen in Const.h\nDo You want to continue anyway?"
						<< "\n**********************************\n";
					system("pause");
				}*/
			}
			else if (MODE == "BCE")
			{
				//Check number of objectives
				if (ONE_OBJ_OVERRIDE != false)
				{
					std::cout << "**********************************\n"
						<< "Wrong parameters chosen. Cannot use BCE for one objective Check Define.h\nProgram will terminate"
						<< "\n**********************************\n";
					system("pause");
					abort();
				}
				//Check number of collectives
				/*if ((n_col_b != 1 && n_col_e != 1))
				{
					std::cout << "**********************************\n"
						<< "Wrong number of collective chosen. Cannot use BCE with MLSGA. Check Const.h\nProgram will terminate"
						<< "\n**********************************\n";
					system("pause");
					abort();
				}*/
				//Check iterations number
				/*if (max_iterations != 300000)
				{
					std::cout << "***************Warning*******************\n"
						<< "Wrong number of iterations chosen. Do not met CEC09 standards. Check Pop size and max gen. Const.h\nDo You want to continue anyway?"
						<< "\n**********************************\n";
					system("pause");
				}*/
			}
			else if (MODE == "HEIA")
			{
				//Check number of objectives
				if (ONE_OBJ_OVERRIDE != false)
				{
					std::cout << "**********************************\n"
						<< "Wrong parameters chosen. Cannot use BCE for one objective Check Define.h\nProgram will terminate"
						<< "\n**********************************\n";
					system("pause");
					abort();
				}
				//Check number of collectives
				/*if ((n_col_b != 1 && n_col_e != 1))
				{
				std::cout << "**********************************\n"
				<< "Wrong number of collective chosen. Cannot use BCE with MLSGA. Check Const.h\nProgram will terminate"
				<< "\n**********************************\n";
				system("pause");
				abort();
				}*/
				//Check iterations number
				/*if (max_iterations != 300000)
				{
					std::cout << "***************Warning*******************\n"
						<< "Wrong number of iterations chosen. Do not met CEC09 standards. Check Pop size and max gen. Const.h\nDo You want to continue anyway?"
						<< "\n**********************************\n";
					system("pause");
				}*/
			}
			else if (MODE == "IBEA")
			{
				//Check number of objectives
				if (ONE_OBJ_OVERRIDE != false)
				{
					std::cout << "**********************************\n"
						<< "Wrong parameters chosen. Cannot use BCE for one objective Check Define.h\nProgram will terminate"
						<< "\n**********************************\n";
					system("pause");
					abort();
				}
			}
			else
			{
				std::cout << "**********************************\n"
					<< "Wrong MODE chosen. Check Define.h\nProgram will terminate"
					<< "\n**********************************\n";
				system("pause");
				abort();
			}
		}
	}

	//Check reinitialisation mode
	if (n_func_e >= 66)
	{
		for (int i = 0; i < Reinit_Mode.size(); i++)
		{
			if (Reinit_Mode[i] == "None")
				;
			else if (Reinit_Mode[i] == "Random")
				;
			else if (Reinit_Mode[i] == "BR")
				;
			else if (Reinit_Mode[i] == "CER")
				;
			else if (Reinit_Mode[i] == "VP")
				;
			else
			{
				std::cout << "**********************************\n"
					<< "Wrong Reinitialisation Mode set. There is no \""
					<< Reinit_Mode[i]
					<< "\" Check Const.h\nProgram will terminate"
					<< "\n**********************************\n";
				system("pause");
				abort();
			}
		}
	}


	//Check fitness saving related functions
	if (FITNESS_ALL == false)
	{
		//Check if video have data
		if (VIDEO == true)
		{
			std::cout << "**********************************\n"
				<< "Wrong VIDEO and FITNESS_ALL parameters chosen. Cannot make video without saving fitness. Check Define.h\nProgram will terminate"
				<< "\n**********************************\n";
			system("pause");
			abort();
		}
	}
	//Check initialisation method
	if (FILE_INPUT == true && SOBOL == true)
	{
		std::cout << "**********************************\n"
			<< "Wrong initialisation parameters chosen. Cannot initialise in 2 different ways. Check Define.h\nProgram will terminate"
			<< "\n**********************************\n";
		system("pause");
		abort();
	}
	//Check if GA parameters are not out of bounds
	for (int i = n_func_b; i <= n_func_e; i++)
	{
		function * f_temp = Func_Type(i);
		delete f_temp;
	}
	for (int i = 1; i <= n_mut; i++)
		Mut_Type(i, "Normal");
	for (int i = 0; i < GA_mode.size(); i++)
	{
		for (int iCoevol = 0; iCoevol < GA_mode[i].size(); iCoevol++)
		{
			std::string MODE = GA_mode[i][iCoevol];
			for (int iC = 1; iC <= n_cross; iC++)
				Cross_Type(iC, MODE);
		}
	}
	for (int i = 1; i <= n_select; i++)
	{		
		Select_Type(i);
	}

	//Check mutation rates
	if (mut_prob_max < mut_prob_min || mut_prob_min < 0 || (mut_prob_max > 1 && ENCODING == "Real"))
	{
		std::cout << "**********************************\n"
			<< "Wrong mutations parameters set. Check Const.h\nProgram will terminate"
			<< "\n**********************************\n";
		system("pause");
		abort();
	}

	//Check crossover rates
	if (cross_prob_max < cross_prob_min || cross_prob_min < 0 || cross_prob_max > 1)
	{
		std::cout << "**********************************\n"
			<< "Wrong mutations parameters set. Check Const.h\nProgram will terminate"
			<< "\n**********************************\n";
		system("pause");
		abort();
	}

	//Check elite size
	if (MLSGA_elite_size >= 1 || MLSGA_elite_size < 0)
	{
		std::cout << "**********************************\n"
			<< "Wrong elite parameters set. Check Const.h\nProgram will terminate"
			<< "\n**********************************\n";
		system("pause");
		abort();
	}

	//Check collective elimiantion
	if (MLSGA_n_col_elim >= MLSGA_n_col_b || MLSGA_n_col_elim < 0)
	{
		std::cout << "**********************************\n"
			<< "Wrong number of eliminated collectives set. Check Const.h\nProgram will terminate"
			<< "\n**********************************\n";
		system("pause");
		abort();
	}

	//Check population size
	if (Pop_size_b <= MLSGA_n_col_e * 10)
	{
		std::cout << "**********************************\n"
			<< "Too small population set. Check Const.h\nProgram will terminate"
			<< "\n**********************************\n";
		system("pause");
		abort();
	}
	if (SOBOL == true)
	{
		if (MLSGA_n_col_e == MLSGA_n_col_b == 1 && Pop_size_e > 32768)
		{
			std::cout << "**********************************\n"
				<< "Too big population set for SOBOL. MAX is 16834. Check Const.h\nProgram will terminate"
				<< "\n**********************************\n";
			system("pause");
			abort();
		}
		else if (Pop_size_e * 2 / (MLSGA_n_col_b) > 32768)
		{
			std::cout << "**********************************\n"
				<< "Too big population set for SOBOL. Check Const.h\nProgram will terminate"
				<< "\n**********************************\n";
			system("pause");
			abort();
		}
	}

}

int Show_Iter()
{
	int n_c_prob = ((cross_prob_max - cross_prob_min) / cross_prob_step) + 1;		//number of iterations for crossover probability rate
	int n_m_prob = ((mut_prob_max - mut_prob_min) / mut_prob_step) + 1;			//number of iterations for mutation probability rate
	int n_pop_dyn = ((Pop_size_e - Pop_size_b) / Pop_size_step) + 1;				//number of iterations for population size
	int n_TGM_mat = 1;						//number of iterations for maternal effect
	int n_TGM_gmat = 1;						//number of iterations for grandmaternal effect

	if (TGM == true)
	{
		n_TGM_mat = ((TGM_m_max - TGM_m_min) / TGM_m_step) + 1;
		n_TGM_gmat = ((TGM_g_max - TGM_g_min) / TGM_g_step) + 1;
	}

	//calculate the number of runs 
	short n_col_e = 1, n_col_b = 1, n_MLS_b = 1, n_MLS_e = 1;
	if (MLSGA_Hybrid == true)
	{
		n_col_b = MLSGA_n_col_b;
		n_col_e = MLSGA_n_col_e;
		n_MLS_b = MLSGA_n_MLS_b;
		n_MLS_e = MLSGA_n_MLS_e;
	}
	long long index_r_max = n_runs * GA_mode.size() * Reinit_Mode.size() * n_select * n_cross * n_mut * (n_func_e - n_func_b + 1)*(n_MLS_e - n_MLS_b + 1)* ((n_col_e - n_col_b) / 2 + 1) * n_pop_dyn*n_c_prob * n_m_prob  *n_TGM_mat * n_TGM_gmat;	//Number of iterations - MAX
	
	//show number of runs
	std::cout << "\nApprox. #" << index_r_max << " iterations;\n";

	//calculate and show approx. time
	long long approx_time = index_r_max * (long long)max_iterations / 5000LL;				//approx time in sec (overall)
	int approx_time_d = approx_time / 86400;			//approx. time - day part	
	int approx_time_h = (approx_time % 86400) / 3600;	//approx. time - hour part
	int approx_time_m = (approx_time % 3600) / 60;		//approx. time - minute part
	int approx_time_s = approx_time % 60;				//approx. time - second  part
	std::cout << "Approx time: ";
	if (approx_time_d > 0)
		std::cout << approx_time_d << "d ";
	if (approx_time_h > 0)
		std::cout << approx_time_h << "h ";
	if (approx_time_m > 0)
		std::cout << approx_time_m << "m ";
	if (approx_time_s > 0)
		std::cout << approx_time_s << "s ";
	std::cout << std::endl;
	system("pause");

	return index_r_max;
}

void IGD_Gen_Calc(std::vector<std::vector<double>> & IGD_s, int index, std::string & name, FILE * Pipe)
{
	std::ofstream IGD_file;			//file for saving the IGD values

	//Open the output file
	IGD_file.open(F_Name_Get("Results/", name, "/") + Name_Get("IGD_o.csv", std::string(), index));

	//Get the max generation size
	int gen_max_size = 0;
	int run_size = IGD_s.size();
	for (int i = 0; i < IGD_s.size(); i++)
		if (gen_max_size < IGD_s[i].size())
			gen_max_size = IGD_s[i].size();
	//Create new inverted vector of normalised size
	std::vector<std::vector<double>> IGD_s_inv(gen_max_size, std::vector<double>(run_size, -2.));

	//Assign the values
	for (int i = 0; i < run_size; i++)
	{
		for (int j = 0; j < IGD_s[i].size(); j++)
			IGD_s_inv[j][i] = IGD_s[i][j];
	}

	//Calculate the average, min, max, std for each generation
	for (int i = 0; i < gen_max_size; i++)
	{
		double IGD_min, IGD_max, IGD_avg, IGD_std;		//storage values
		int IGD_size = 0;
		IGD_max = IGD_avg = IGD_std = 0;
		IGD_min = 0;

		//Calculate the average, min, max for each run
		for (int j = 0; j < run_size; j++)
		{
			//if empty do not count
			if (IGD_s_inv[i][j] <= -1)
				continue;
			//Check if is max
			if (j == 0)
				IGD_max = IGD_min = IGD_s_inv[i][j];
			else if (IGD_max < IGD_s_inv[i][j])
				IGD_max = IGD_s_inv[i][j];
			//Check if is min
			else if (IGD_min > IGD_s_inv[i][j])
				IGD_min = IGD_s_inv[i][j];

			//Add to avg
			IGD_avg += IGD_s_inv[i][j];
			IGD_size++;
		}

		//Calcualte the actual average
		IGD_avg /= IGD_size;

		//Calculate the std for each run
		for (int j = 0; j < run_size; j++)
		{
			//if empty do not count
			if (IGD_s_inv[i][j] <= -1)
				continue;
			IGD_std += pow(IGD_s_inv[i][j] - IGD_avg, 2);
		}
		//Calculate the actual std
		IGD_std = sqrt(IGD_std / IGD_size);

		//save values to the file
		IGD_file << i << " " << IGD_min << " " << IGD_max << " " << IGD_avg << " " << IGD_std << std::endl;
	}

	//Close the IGD file
	IGD_file.close();

	//print the graph
	std::string print_temp;
	//Set styles
	print_temp += "set style line 5 lc rgb \"blue\" pt 7 ps 1 \n ";	//Style line for min
	print_temp += "set style line 6 lc rgb \"green\" pt 7 ps 1 \n ";	//Style line for avg
	print_temp += "set style line 7 lc rgb \"red\" pt 7 ps 1 \n ";	//Style line for max
	fprintf(Pipe, print_temp.c_str());
	print_temp.clear();

	//Set labels and titles
	fprintf(Pipe, "set xlabel \"Generation\"\nset ylabel \"IGD\"\n");
	fprintf(Pipe, "set title \"IGD over generations\"\n");

	//Set output and inputs
	print_temp = "set term png size ";
	print_temp += std::to_string(frame_w);
	print_temp += ",";
	print_temp += std::to_string(frame_h);
	print_temp += "\nset output ";
	print_temp += F_Name_Get("\"Results\\\\", name, "\\\\");
	print_temp += Name_Get("IGD_o.jpg", std::string(), index);
	print_temp += "\nset xr[*:*] \nset yr[*:*] \nplot ";
	print_temp += F_Name_Get("\"Results\\\\", name, "\\\\");
	print_temp += Name_Get("IGD_o.csv", std::string(), index);
	print_temp += "\" using 1:2 with points ls 5";
	print_temp += " title \"IGD_min\",";
	print_temp += F_Name_Get("\"Results\\\\", name, "\\\\");
	print_temp += Name_Get("IGD_o.csv", std::string(), index);
	print_temp += "\" using 1:4 with points ls 6";
	print_temp += " title \"IGD_avg\",";
	print_temp += F_Name_Get("\"Results\\\\", name, "\\\\");
	print_temp += Name_Get("IGD_o.csv", std::string(), index);
	print_temp += "\" using 1:3 with points ls 7";
	print_temp += " title \"IGD_max\"";
	print_temp += "\n";
	fprintf(Pipe, print_temp.c_str());
	print_temp.clear();
	
	fflush(Pipe);
}

void HV_Gen_Calc(std::vector<std::vector<double>> & HV_s, int index, std::string & name, FILE * Pipe)
{
	std::ofstream HV_file;			//file for saving the HV values

									//Open the output file
	HV_file.open(F_Name_Get("Results/", name, "/") + Name_Get("HV_o.csv", std::string(), index));

	//Get the max generation size
	int gen_max_size = 0;
	int run_size = HV_s.size();
	for (int i = 0; i < HV_s.size(); i++)
		if (gen_max_size < HV_s[i].size())
			gen_max_size = HV_s[i].size();
	//Create new inverted vector of normalised size
	std::vector<std::vector<double>> HV_s_inv(gen_max_size, std::vector<double>(run_size, -2.));

	//Assign the values
	for (int i = 0; i < run_size; i++)
	{
		for (int j = 0; j < HV_s[i].size(); j++)
			HV_s_inv[j][i] = HV_s[i][j];
	}

	//Calculate the average, min, max, std for each generation
	for (int i = 0; i < gen_max_size; i++)
	{
		double HV_min, HV_max, HV_avg, HV_std;		//storage values
		int HV_size = 0;
		HV_max = HV_avg = HV_std = 0;
		HV_min = 0;

		//Calculate the average, min, max for each run
		for (int j = 0; j < run_size; j++)
		{
			//if empty do not count
			if (HV_s_inv[i][j] <= -1)
				continue;
			//Check if is max
			if (j == 0)
				HV_max = HV_min = HV_s_inv[i][j];
			else if (HV_max < HV_s_inv[i][j])
				HV_max = HV_s_inv[i][j];
			//Check if is min
			else if (HV_min > HV_s_inv[i][j])
				HV_min = HV_s_inv[i][j];

			//Add to avg
			HV_avg += HV_s_inv[i][j];
			HV_size++;
		}

		//Calcualte the actual average
		HV_avg /= HV_size;

		//Calculate the std for each run
		for (int j = 0; j < run_size; j++)
		{
			//if empty do not count
			if (HV_s_inv[i][j] <= -1)
				continue;
			HV_std += pow(HV_s_inv[i][j] - HV_avg, 2);
		}
		//Calculate the actual std
		HV_std = sqrt(HV_std / HV_size);

		//save values to the file
		HV_file << i << " " << HV_min << " " << HV_max << " " << HV_avg << " " << HV_std << std::endl;
	}

	//Close the HV file
	HV_file.close();

	//print the graph
	std::string print_temp;
	//Set styles
	print_temp += "set style line 5 lc rgb \"blue\" pt 7 ps 1 \n ";	//Style line for min
	print_temp += "set style line 6 lc rgb \"green\" pt 7 ps 1 \n ";	//Style line for avg
	print_temp += "set style line 7 lc rgb \"red\" pt 7 ps 1 \n ";	//Style line for max
	fprintf(Pipe, print_temp.c_str());
	print_temp.clear();

	//Set labels and titles
	fprintf(Pipe, "set xlabel \"Generation\"\nset ylabel \"HV\"\n");
	fprintf(Pipe, "set title \"HV over generations\"\n");

	//Set output and inputs
	print_temp = "set term png size ";
	print_temp += std::to_string(frame_w);
	print_temp += ",";
	print_temp += std::to_string(frame_h);
	print_temp += "\nset output ";
	print_temp += F_Name_Get("\"Results\\\\", name, "\\\\");
	print_temp += Name_Get("HV_o.jpg", std::string(), index);
	print_temp += "\nset xr[*:*] \nset yr[*:*] \nplot ";
	print_temp += F_Name_Get("\"Results\\\\", name, "\\\\");
	print_temp += Name_Get("HV_o.csv", std::string(), index);
	print_temp += "\" using 1:2 with points ls 5";
	print_temp += " title \"HV_min\",";
	print_temp += F_Name_Get("\"Results\\\\", name, "\\\\");
	print_temp += Name_Get("HV_o.csv", std::string(), index);
	print_temp += "\" using 1:4 with points ls 6";
	print_temp += " title \"HV_avg\",";
	print_temp += F_Name_Get("\"Results\\\\", name, "\\\\");
	print_temp += Name_Get("HV_o.csv", std::string(), index);
	print_temp += "\" using 1:3 with points ls 7";
	print_temp += " title \"HV_max\"";
	print_temp += "\n";
	fprintf(Pipe, print_temp.c_str());
	print_temp.clear();

	fflush(Pipe);
}

void MLS_Col_Fit_Create(std::vector<collective> & col, short MLS)
{
	short nobj = col[0].FCode_Show()[0].Objs();
	short ncol = col.size();

	if (MLS == 0)
	{
		std::vector<short> temp_fit_index_level;
		for (short i_obj = 1; i_obj <= nobj; i_obj++)
			temp_fit_index_level.push_back(i_obj);

		std::vector<std::vector<short>> temp_fit_index{ temp_fit_index_level,temp_fit_index_level };

		col[0].fit_index = temp_fit_index;

		return;
	}

	if (nobj == 2)
	{
		if (MLS == 1)
		{
			std::vector<std::vector<short>> temp_fit_index{ {1,2},{1,2} };
			for (short i_col = 0; i_col < ncol; i_col++)
			{
				col[i_col].fit_index = temp_fit_index;
			}
		}
		else if (MLS == 2)
		{
			std::vector<std::vector<short>> temp_fit_index_n{ { fit_index_sel },{ fit_index_col_sel } };
			std::vector<std::vector<short>> temp_fit_index{ { fit_index_sel, fit_index_col_sel },{ fit_index_col_sel } };
			for (short i_col = 0; i_col < ncol; i_col++)
			{
				if (col[i_col].Mode_Show() == "Normal")
					col[i_col].fit_index = temp_fit_index_n;
				else
					col[i_col].fit_index = temp_fit_index;
			}
		}
		else if (MLS == 3)
		{
			std::vector<std::vector<short>> temp_fit_index_n{ { fit_index_col_sel },{ fit_index_sel } };
			std::vector<std::vector<short>> temp_fit_index{ { fit_index_sel, fit_index_col_sel },{ fit_index_sel } };
			for (short i_col = 0; i_col < ncol; i_col++)
			{
				if (col[i_col].Mode_Show() == "Normal")
					col[i_col].fit_index = temp_fit_index_n;
				else
					col[i_col].fit_index = temp_fit_index;
			}
		}
		else if (MLS == 7)
		{
			std::vector<std::vector<short>> temp_fit_index1{ { fit_index_sel, fit_index_col_sel },{ fit_index_sel, fit_index_col_sel } };
			std::vector<std::vector<short>> temp_fit_index2{ { fit_index_sel, fit_index_col_sel },{fit_index_col_sel } };
			std::vector<std::vector<short>> temp_fit_index2_n{ { fit_index_sel},{ fit_index_col_sel } };
			std::vector<std::vector<short>> temp_fit_index3{ { fit_index_sel, fit_index_col_sel },{ fit_index_sel } };
			std::vector<std::vector<short>> temp_fit_index3_n{ { fit_index_col_sel },{ fit_index_sel } };

			for (short i_col = 0; i_col < ncol; i_col++)
			{
				short col_index = col[i_col].Index_Show();
				if (col_index % 3 == 1)
				{
					col[i_col].fit_index = temp_fit_index1;
				}
				else if (col_index % 3 == 0)
				{
					if (col[i_col].Mode_Show() == "Normal")
						col[i_col].fit_index = temp_fit_index2_n;
					else
						col[i_col].fit_index = temp_fit_index2;
				}
				else
				{
					if (col[i_col].Mode_Show() == "Normal")
						col[i_col].fit_index = temp_fit_index3_n;
					else
						col[i_col].fit_index = temp_fit_index3;
				}
			}
		}
		else if (MLS == 8)
		{
			std::vector<std::vector<short>> temp_fit_index1{ { fit_index_sel, fit_index_col_sel },{ fit_index_sel, fit_index_col_sel } };
			std::vector<std::vector<short>> temp_fit_index2{ { fit_index_sel, fit_index_col_sel },{fit_index_col_sel } };
			std::vector<std::vector<short>> temp_fit_index2_n{ { fit_index_sel},{ fit_index_sel, fit_index_col_sel } };
			std::vector<std::vector<short>> temp_fit_index3{ { fit_index_sel, fit_index_col_sel },{ fit_index_sel } };
			std::vector<std::vector<short>> temp_fit_index3_n{ { fit_index_col_sel },{ fit_index_sel, fit_index_col_sel } };

			for (short i_col = 0; i_col < ncol; i_col++)
			{
				short col_index = col[i_col].Index_Show();
				if (col_index % 3 == 1)
				{
					col[i_col].fit_index = temp_fit_index1;
				}
				else if (col_index % 3 == 0)
				{
					if (col[i_col].Mode_Show() == "Normal")
						col[i_col].fit_index = temp_fit_index2_n;
					else
						col[i_col].fit_index = temp_fit_index2;
				}
				else
				{
					if (col[i_col].Mode_Show() == "Normal")
						col[i_col].fit_index = temp_fit_index3_n;
					else
						col[i_col].fit_index = temp_fit_index3;
				}
			}
		}
		else if (MLS == 9)
		{
			std::vector<std::vector<short>> temp_fit_index1{ { fit_index_sel, fit_index_col_sel },{ fit_index_sel, fit_index_col_sel } };
			std::vector<std::vector<short>> temp_fit_index2{ { fit_index_sel, fit_index_col_sel },{fit_index_col_sel } };
			std::vector<std::vector<short>> temp_fit_index2_n{ { fit_index_sel, fit_index_col_sel } ,{ fit_index_col_sel } };
			std::vector<std::vector<short>> temp_fit_index3{ { fit_index_sel, fit_index_col_sel },{ fit_index_sel } };
			std::vector<std::vector<short>> temp_fit_index3_n{ { fit_index_sel, fit_index_col_sel } ,{ fit_index_sel } };

			for (short i_col = 0; i_col < ncol; i_col++)
			{
				short col_index = col[i_col].Index_Show();
				if (col_index % 3 == 1)
				{
					col[i_col].fit_index = temp_fit_index1;
				}
				else if (col_index % 3 == 0)
				{
					if (col[i_col].Mode_Show() == "Normal")
						col[i_col].fit_index = temp_fit_index2_n;
					else
						col[i_col].fit_index = temp_fit_index2;
				}
				else
				{
					if (col[i_col].Mode_Show() == "Normal")
						col[i_col].fit_index = temp_fit_index3_n;
					else
						col[i_col].fit_index = temp_fit_index3;
				}
			}
		}
		else
			abort();
	}
	else
	{
		if (MLS == 7)
		{
			std::vector<short> temp_fit;
			for (short i = 0; i < nobj; i++)
				temp_fit.push_back(i + 1);
			std::vector<std::vector<short>> temp_fit_index1{ temp_fit,temp_fit };
			for (short i_col = 0; i_col < ncol; i_col++)
			{
				short col_index = col[i_col].Index_Show() - 1;

				short sw = col_index % (nobj + 1);
				if (sw == 0)
				{
					col[i_col].fit_index = temp_fit_index1;
				}
				else
				{
					col[i_col].fit_index = { temp_fit,{sw} };
				}
			}
		}
		else
			abort();
	}
	/*else if (nobj == 3)
	{
		if (MLS == 1)
		{
			std::vector<std::vector<short>> temp_fit_index{ { 1,2,3 },{ 1,2,3 } };
			for (short i_col = 0; i_col < ncol; i_col++)
			{
				col[i_col].fit_index = temp_fit_index;
			}
		}
		else if (MLS == 2)
		{
			std::vector<std::vector<short>> temp_fit_index1{ { 1, 2, 3 },{ 1 } };
			std::vector<std::vector<short>> temp_fit_index2{ { 1, 2, 3 },{ 2 } };
			std::vector<std::vector<short>> temp_fit_index3{ { 1, 2, 3 },{ 3 } };
			for (short i_col = 0; i_col < ncol; i_col++)
			{
				short col_index = col[i_col].Index_Show();
				if (col_index % 3 == 1)
				{
					col[i_col].fit_index = temp_fit_index1;
				}
				else if (col_index % 3 == 0)
				{
					col[i_col].fit_index = temp_fit_index2;
				}
				else
				{
					col[i_col].fit_index = temp_fit_index3;
				}
			}
		}
		else if (MLS == 3)
		{
			std::vector<std::vector<short>> temp_fit_index1{ { 1, 2 },{ 2,3 } };
			std::vector<std::vector<short>> temp_fit_index2{ { 2, 3 },{ 1,3 } };
			std::vector<std::vector<short>> temp_fit_index3{ { 1, 3 },{ 1,2 } };
			for (short i_col = 0; i_col < ncol; i_col++)
			{
				short col_index = col[i_col].Index_Show();
				if (col_index % 3 == 1)
				{
					col[i_col].fit_index = temp_fit_index1;
				}
				else if (col_index % 3 == 0)
				{
					col[i_col].fit_index = temp_fit_index2;
				}
				else
				{
					col[i_col].fit_index = temp_fit_index3;
				}
			}
		}
		else if (MLS == 7)
		{
			std::vector<std::vector<short>> temp_fit_index1{ { 1,2,3 },{ 1,2,3 } };
			std::vector<std::vector<short>> temp_fit_index2{ { 1, 2 },{ 2,3 } };
			std::vector<std::vector<short>> temp_fit_index3{ { 2, 3 },{ 1,3 } };
			std::vector<std::vector<short>> temp_fit_index4{ { 1, 3 },{ 1,2 } };
			for (short i_col = 0; i_col < ncol; i_col++)
			{
				short col_index = col[i_col].Index_Show();
				if (col_index % 4 == 1)
				{
					col[i_col].fit_index = temp_fit_index1;
				}
				else if (col_index % 4 == 2)
				{
					col[i_col].fit_index = temp_fit_index2;
				}
				else if (col_index % 4 == 3)
				{
					col[i_col].fit_index = temp_fit_index3;
				}
				else
				{
					col[i_col].fit_index = temp_fit_index4;
				}
			}
		}
		else if (MLS == 8)
		{
			std::vector<std::vector<short>> temp_fit_index1{ { 1,2,3 },{ 1,2,3 } };
			std::vector<std::vector<short>> temp_fit_index2{ { 1, 2 },{ 3 } };
			std::vector<std::vector<short>> temp_fit_index3{ { 2, 3 },{ 1 } };
			std::vector<std::vector<short>> temp_fit_index4{ { 1, 3 },{ 2 } };
			for (short i_col = 0; i_col < ncol; i_col++)
			{
				short col_index = col[i_col].Index_Show();
				if (col_index % 4 == 1)
				{
					col[i_col].fit_index = temp_fit_index1;
				}
				else if (col_index % 4 == 2)
				{
					col[i_col].fit_index = temp_fit_index2;
				}
				else if (col_index % 4 == 3)
				{
					col[i_col].fit_index = temp_fit_index3;
				}
				else
				{
					col[i_col].fit_index = temp_fit_index4;
				}
			}
		}
		else if (MLS == 9)
		{
		std::vector<std::vector<short>> temp_fit_index1{ { 1,2,3 },{ 1,2,3 } };
		std::vector<std::vector<short>> temp_fit_index2{ { 1,2,3 },{ 3 } };
		std::vector<std::vector<short>> temp_fit_index3{ { 1,2,3 },{ 1 } };
		std::vector<std::vector<short>> temp_fit_index4{ { 1,2,3 },{ 2 } };
		for (short i_col = 0; i_col < ncol; i_col++)
		{
			short col_index = col[i_col].Index_Show();
			if (col_index % 4 == 1)
			{
				col[i_col].fit_index = temp_fit_index1;
			}
			else if (col_index % 4 == 2)
			{
				col[i_col].fit_index = temp_fit_index2;
			}
			else if (col_index % 4 == 3)
			{
				col[i_col].fit_index = temp_fit_index3;
			}
			else
			{
				col[i_col].fit_index = temp_fit_index4;
			}
		}
		}
		else
			abort();
	}
	else if (nobj == 5)
	{
		if (MLS == 7)
		{
			std::vector<std::vector<short>> temp_fit_index1{ { 1,2,3,4,5 },{ 1,2,3,4,5 } };
			std::vector<std::vector<short>> temp_fit_index2{ { 1, 2, 3 },{ 3,4,5 } };
			std::vector<std::vector<short>> temp_fit_index3{ { 2, 3, 4 },{ 1,4,5 } };
			std::vector<std::vector<short>> temp_fit_index4{ { 3, 4, 5 },{ 1,2,5 } };
			std::vector<std::vector<short>> temp_fit_index5{ { 1, 4, 5 },{ 1,2,3 } };
			std::vector<std::vector<short>> temp_fit_index6{ { 1, 2, 5 },{ 2,3,4 } };
			std::vector<std::vector<short>> temp_fit_index7{ { 1, 2, 3 },{ 1,4,5 } };
			std::vector<std::vector<short>> temp_fit_index8{ { 3, 4, 5 },{ 1,2,3 } };
			for (short i_col = 0; i_col < ncol; i_col++)
			{
				short col_index = col[i_col].Index_Show();
				if (col_index % 8 == 1)
				{
					col[i_col].fit_index = temp_fit_index1;
				}
				else if (col_index % 8 == 2)
				{
					col[i_col].fit_index = temp_fit_index2;
				}
				else if (col_index % 8 == 3)
				{
					col[i_col].fit_index = temp_fit_index3;
				}
				else if (col_index % 8 == 4)
				{
					col[i_col].fit_index = temp_fit_index4;
				}
				else if (col_index % 8 == 5)
				{
					col[i_col].fit_index = temp_fit_index5;
				}
				else if (col_index % 8 == 6)
				{
					col[i_col].fit_index = temp_fit_index6;
				}
				else if (col_index % 8 == 7)
				{
					col[i_col].fit_index = temp_fit_index7;
				}
				else
				{
					col[i_col].fit_index = temp_fit_index8;
				}
			}
		}
		else if (MLS == 8)
		{
			std::vector<std::vector<short>> temp_fit_index1{ { 1,2,3,4,5 },{ 1,2,3,4,5 } };
			std::vector<std::vector<short>> temp_fit_index2{ { 1,2,3,4 },{ 5 } };
			std::vector<std::vector<short>> temp_fit_index3{ { 1,2,3,5 },{ 4 } };
			std::vector<std::vector<short>> temp_fit_index4{ { 1,2,4,5 },{ 3 } };
			std::vector<std::vector<short>> temp_fit_index5{ { 1,3,4,5 },{ 2 } };
			std::vector<std::vector<short>> temp_fit_index6{ { 2,3,4,5 },{ 1 } };
			std::vector<std::vector<short>> temp_fit_index7{ { 1,3,4 },{ 2,5 } };
			std::vector<std::vector<short>> temp_fit_index8{ { 2,3,5 },{ 1,4 } };
			for (short i_col = 0; i_col < ncol; i_col++)
			{
				short col_index = col[i_col].Index_Show();
				if (col_index % 8 == 1)
				{
					col[i_col].fit_index = temp_fit_index1;
				}
				else if (col_index % 8 == 2)
				{
					col[i_col].fit_index = temp_fit_index2;
				}
				else if (col_index % 8 == 3)
				{
					col[i_col].fit_index = temp_fit_index3;
				}
				else if (col_index % 8 == 4)
				{
					col[i_col].fit_index = temp_fit_index4;
				}
				else if (col_index % 8 == 5)
				{
					col[i_col].fit_index = temp_fit_index5;
				}
				else if (col_index % 8 == 6)
				{
					col[i_col].fit_index = temp_fit_index6;
				}
				else if (col_index % 8 == 7)
				{
					col[i_col].fit_index = temp_fit_index7;
				}
				else
				{
					col[i_col].fit_index = temp_fit_index8;
				}
			}
		}
		else if (MLS == 9)
		{
			std::vector<std::vector<short>> temp_fit_index1{ { 1,2,3,4,5 },{ 1,2,3,4,5 } };
			std::vector<std::vector<short>> temp_fit_index2{ { 1,2,3,4,5 },{ 5 } };
			std::vector<std::vector<short>> temp_fit_index3{ { 1,2,3,4,5 },{ 4 } };
			std::vector<std::vector<short>> temp_fit_index4{ { 1,2,3,4,5 },{ 3 } };
			std::vector<std::vector<short>> temp_fit_index5{ { 1,2,3,4,5 },{ 2 } };
			std::vector<std::vector<short>> temp_fit_index6{ { 1,2,3,4,5 },{ 1 } };
			std::vector<std::vector<short>> temp_fit_index7{ { 1,3,4 },{ 2,5 } };
			std::vector<std::vector<short>> temp_fit_index8{ { 2,3,5 },{ 1,4 } };
			for (short i_col = 0; i_col < ncol; i_col++)
			{
				short col_index = col[i_col].Index_Show();
				if (col_index % 8 == 1)
				{
					col[i_col].fit_index = temp_fit_index1;
				}
				else if (col_index % 8 == 2)
				{
					col[i_col].fit_index = temp_fit_index2;
				}
				else if (col_index % 8 == 3)
				{
					col[i_col].fit_index = temp_fit_index3;
				}
				else if (col_index % 8 == 4)
				{
					col[i_col].fit_index = temp_fit_index4;
				}
				else if (col_index % 8 == 5)
				{
					col[i_col].fit_index = temp_fit_index5;
				}
				else if (col_index % 8 == 6)
				{
					col[i_col].fit_index = temp_fit_index6;
				}
				else if (col_index % 8 == 7)
				{
					col[i_col].fit_index = temp_fit_index7;
				}
				else
				{
					col[i_col].fit_index = temp_fit_index8;
				}
			}
		}
		else
			abort();
	}
	else
		abort();*/
}




