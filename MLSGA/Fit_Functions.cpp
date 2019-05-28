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

#include "Fit_Functions.h"
#include "Support_Functions.h"
#include "Define.h"
#include "Imported/WFG.h"
#include <sstream>

/*
*Min fitness calculation for one objective optimisation*
@param pareto_front real pareto front for which fitness will be calculated
*/
double function::Min_Fitness_Get(const std::vector<std::vector<double>> & pareto_front)
{
	//Check if it is one objective
	if (ONE_OBJ_OVERRIDE != true)
		abort();
	double min_fit;						//min fitness for current function
	//loop for whole real pareto front
	for (int i = 0; i < pareto_front.size(); i++)
	{
		double fit_temp = 0;					//temporary fitness/average of fitnesses
		//loop for number of objectives
		for (short j = 0; j < num_objs; j++)
		{
			//calculate average fitness
			fit_temp += pareto_front[i][j] / (double)num_objs;
		}
		//check if fit_temp is smaller than current min, if is assign to min
		if (i == 0)
			min_fit = fit_temp;
		else if (fit_temp < min_fit)
			min_fit = fit_temp;
	}
	return min_fit;
}

std::vector<std::vector<double>> function::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving
	std::ifstream PF_real_input;
								//check if the index is correct
	PF_real_input.open("Input/PF/" + name_func + ".dat");
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

	if (!PF_real_input.good())
	{
		std::cout << "ERROR#12: FUNCTION - PF_READ; FILE NOT OPEN";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}
	double temp_val;
	std::vector<double> PF_real_val_temp;
	while (PF_real_input >> temp_val)
	{
		//push value to temporary vector
		PF_real_val_temp.push_back(temp_val);
		//Push value to file
		PF_real << temp_val;
		//Send values to output vector
		if (PF_real_val_temp.size() == Objs())			//if size of temprary vector is equal to nubmer of objectives
		{
			//push to output vector
			PF_real_val.push_back(PF_real_val_temp);
			//clear the temporary vector
			PF_real_val_temp.clear();
			//end line in output file
			PF_real << std::endl;
		}
		else
			//push space to the file
			PF_real << " ";
	}

	//close the PF files
	PF_real.close();
	PF_real_input.close();

	return PF_real_val;
}

void function::T_Vector_Update(bool s) 
{
	//get the size of the current vector
	short v_size = t_vector.size();

	//if t_vector have to be zeroed 
	if (!s)
	{
		//create the new vector of given size
		t_vector = std::vector<float>(v_size, 0.f);

		//terminate the function
		return;
	}
	//Choose which t value have to be updated by random
	int t_vector_val = Random_I(0, v_size - 1);

	//update the t_vector
	t_vector[t_vector_val] += 1.f / n_steps_dyn;
}

/********************************************
			UNCONSTRAINED FUNCTIONS
*********************************************/

/**************ZDT1****************************/
std::vector<double> ZDT1::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;			//fitness vector - output
	double sum_f2 = 0;						//temp value for the fitness calculation
	double g_x = 0;							//temp values for the fitness calculation

											//calculate temp values
	for (int i = 1; i < Vars(); i++)
		sum_f2 += code[i];
	g_x = 1 + 9 * sum_f2 / ((double)Vars() - 1);

	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0]);								//fitness 1
	fitness.push_back(g_x*(1 - sqrt((code[0] / g_x))));				//fitness 2

																	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}
void ZDT1::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set boundaries for variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 1,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 7,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> ZDT1::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1 - sqrt(i));					//PF value
		std::vector<double> temp_vect;				//temporary vector of fitness

													//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}
/**************ZDT2****************************/
std::vector<double> ZDT2::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;				//fitness vector - output
	double sum_f2 = 0;							//temp value for the fitness calculation
	double g_x = 0;								//temp value for the fitness calculation

												//calculate temp values
	for (int i = 1; i < Vars(); i++)
		sum_f2 += code[i];
	g_x = 1 + 9 * sum_f2 / ((double)Vars() - 1);

	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0]);										//fitness 1
	fitness.push_back(g_x*(1 - pow((code[0] / g_x), 2)));				//fitness 2

																		//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}
void ZDT2::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set boundaries for variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}
	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 1,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 7,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> ZDT2::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;			//ofstream file for the PF saving

									//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1 - pow(i, 2));		//PF value
		std::vector<double> temp_vect;		//temporary vector of fitness

											//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}
/**************ZDT3****************************/
std::vector<double> ZDT3::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;				//fitness vector - output
	double sum_f2 = 0;							//temp value for the fitness calculation
	double g_x = 0;								//temp value for the fitness calculation
	double f2, f1;								//temp values for the fitness

												//calculate fitness 1
	f1 = code[0];

	//calculate temp values
	for (int i = 1; i < Vars(); i++)
		sum_f2 += code[i];
	g_x = 1 + 9 * sum_f2 / ((double)Vars() - 1);

	//calculate the fitness and push it to the fitness vector
	f2 = g_x*(1 - sqrt(f1 / g_x) - f1 / g_x*sin(10 * pi*f1));
	fitness.push_back(f1);									//fitness 1
	fitness.push_back(f2);				//fitness 2

										//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}
void ZDT3::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set boundaries for variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}
	//check if the amount of boudaries is correct

	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 1,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 7,-1 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> ZDT3::Plot_PF(int indexr, int size)				//non continous
{
	std::ofstream PF_real;			//ofstream file for the PF saving

									//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//calculate the step of the i value for the loop
	double i_step = ((0.8518328654 - 0.8233317983) + (0.6525117038 - 0.6183967944) + (0.4538821041 - 0.4093136748) + (0.2577623634 - 0.1822287280) + (0.0830015349)) / (double)(size - 1);

	//calculate the real PF points
	for (double i = 0; i <= 0.0830015349; i += i_step)
	{
		//calculate the PF value
		double temp = (1 - sqrt(i) - i*sin(10 * pi*i));		//PF value

		std::vector<double> temp_vect;		//temporary vector of fitness

											//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//calculate the real PF points
	for (double i = 0.1822287280; i <= 0.2577623634; i += i_step)
	{
		//calculate the PF value
		double temp = (1 - sqrt(i) - i*sin(10 * pi*i));		//PF value

		std::vector<double> temp_vect;		//temporary vector of fitness

											//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		//push the fitness vector to the PF vector

		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//calculate the real PF points
	for (double i = 0.4093136748; i <= 0.4538821041; i += i_step)
	{
		//calculate the PF value
		double temp = (1. - sqrt(i) - i*sin(10 * pi*i));		//PF value

		std::vector<double> temp_vect;		//temporary vector of fitness

											//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//calculate the real PF points
	for (double i = 0.6183967944; i <= 0.6525117038; i += i_step)
	{
		//calculate the PF value
		double temp = (1 - sqrt(i) - i*sin(10 * pi*i));		//PF value

		std::vector<double> temp_vect;		//temporary vector of fitness

											//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//calculate the real PF points
	for (double i = 0.8233317983; i <= 0.8518328654; i += i_step)
	{
		//calculate the PF value
		double temp = (1 - sqrt(i) - i*sin(10 * pi*i));		//PF value

		std::vector<double> temp_vect;		//temporary vector of fitness

											//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();
	return PF_real_val;
}
/**************ZDT4****************************/
std::vector<double> ZDT4::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;			//fitness vector - output
	double sum_f2 = 0;						//temp values for the fitness calculation
	double g_x = 0;							//temp values for the fitness calculation
	double f2, f1;							//temp values of fitness

											//calculate temp values
	for (int i = 1; i < Vars(); i++)
		sum_f2 += (pow(code[i], 2) - 10 * cos(4 * pi*code[i]));
	g_x = 1 + 10 * ((double)Vars() - 1) + sum_f2;

	//calculate the fitness and push it to the fitness vector
	f1 = code[0];
	f2 = g_x*(1 - sqrt(f1 / g_x));
	fitness.push_back(f1);									//fitness 1
	fitness.push_back(f2);									//fitness 2

															//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}
void ZDT4::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 5;
		temp.lower = -5;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 1,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 150,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> ZDT4::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1 - sqrt(i));			//PF value
		std::vector<double> temp_vect;		//temporary vector of fitness

											//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}
/**************ZDT5****************************/
std::vector<double> ZDT5::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;					//fitness vector - output
	double sum_f2 = 0;								//temp value for the fitness calculation
	double g_x = 0;									//temp value for the fitness calculation

													//calculate the fitness and push it to the fitness vector
	for (int i = 1; i < Vars(); i++)
		sum_f2 += code[i];
	g_x = 1 + 9 * sum_f2 / ((double)Vars() - 1);

	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0]);										//fitness 1
	fitness.push_back(g_x*(1 - sqrt((code[0] / g_x))));				//fitness 2

																	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}
void ZDT5::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set boundaries for variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 1,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> ZDT5::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;				//ofstream file for the PF saving

										//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1 - sqrt(i));			//PF value
		std::vector<double> temp_vect;		//temporary vector of fitness

											//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}
/**************ZDT6****************************/
std::vector<double> ZDT6::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;				//fitness vector - output
	double sum_f2 = 0;							//temp value for the fitness calculation
	double f1, f2;								//temp values of fitness 
	double g_x = 0;								//temp value for the fitness calculation

												//calculate temp value
	for (int i = 1; i < Vars(); i++)
		sum_f2 += code[i];
	g_x = 1 + 9 * pow((sum_f2 / ((double)Vars() - 1)), 0.25);

	//calculate the fitness and push it to the fitness vector
	f1 = 1 - exp(-4.0 * code[0])*pow(sin(6 * pi*code[0]), 6);
	f2 = g_x * (1 - pow((f1 / g_x), 2));			//From paper
													//f2 = 1 - pow((f1 / g_x), 2);					//From website - ETH - wrong!!
	fitness.push_back(f1);										//fitness 1
	fitness.push_back(f2);										//fitness 2

																//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}
void ZDT6::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set boundaries for variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 1,0};								// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> ZDT6::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;			//ofstream file for the PF saving

									//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//calculate the step for the loop
	double i_min = 0.2807753191;							//min i value
	double i_max = 1.0;										//max i value
	double i_step = (i_max - i_min) / (double)(size - 1);	//step of the i value change

															//calculate the real PF points
	for (double i = i_min; i <= i_max; i += i_step)
	{
		//calculate the PF value
		double temp = (1.f - (i*i));			//PF value
		std::vector<double> temp_vect;		//temporary vector of fitness

											//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}
/**************DTLZ1****************************/
std::vector<double> DTLZ1::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;			//fitness vector - output
	double sum_f = 0;						//temp value for the fitness calculation
	double g_x = 0;							//temp values for the fitness calculation

	short nobj = this->Objs();



	for (int i = nobj - 1; i < Vars(); i++)
	{
		sum_f += pow((code[i] - 0.5), 2.0) - cos(20 * pi*(code[i] - 0.5));
	}
	g_x = 100 * (sum_f + Vars() - nobj + 1) + 1.0;
	sum_f = g_x;
	for (int j = 0; j<nobj - 1; j++)
	{
		sum_f = sum_f * code[j];
	}
	fitness.push_back(0.5 * sum_f);

	for (int i = 1; i < nobj; i++)
	{
		sum_f = g_x;
		for (int j = 0; j<nobj - 1 - i; j++)
		{
			sum_f = sum_f * code[j];
		}
		sum_f = sum_f * (1.0 - code[nobj - 1 - i]);
		fitness.push_back(0.5 * sum_f);
	}


																	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}
void DTLZ1::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set boundaries for variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 100,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> DTLZ1::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

	//Get the uniform vector
	PF_real_val = Uniform_Weights_Generate(size, Objs(), 3);
														
														//calculate the real PF points
	for (double i = 0; i < PF_real_val.size(); i++)
	{
		for (int j = 0; j < Objs(); j++)
		{
			PF_real_val[i][j] *= 0.5;
			if (indexr >= 0)
				PF_real << PF_real_val[i][j] << " ";
		}
		if (indexr >= 0)
			PF_real << std::endl;
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**************DTLZ2****************************/
std::vector<double> DTLZ2::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;				//fitness vector - output
	double sum_f = 0;							//temp value for the fitness calculation
	double g_x = 0;								//temp value for the fitness calculation
	short nobj = Objs();
												//calculate temp values
	for (int i = nobj - 1; i<Vars(); i++)
	{
		sum_f += pow((code[i] - 0.5), 2.0);
	}
	g_x = 1.0 + sum_f;
	sum_f = g_x;
	for (int j = 0; j<nobj - 1; j++)
	{
		sum_f = sum_f * cos(code[j] * pi / 2.0);
	}
	//Calculate and push fitness
	fitness.push_back(sum_f);	
	for (int i = 1; i<nobj; i++)
	{
		sum_f = g_x;
		for (int j = 0; j<nobj - 1 - i; j++)
		{
			sum_f = sum_f * cos(code[j] * pi / 2.0);
		}
		sum_f = sum_f * sin(code[nobj - 1 - i] * pi / 2.0);
		fitness.push_back(sum_f);
	}

																		//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}
void DTLZ2::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set boundaries for variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}
	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 1,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 2.5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> DTLZ2::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//Get the uniform vector
	PF_real_val = Uniform_Weights_Generate(size, Objs(),1);

	//calculate the real PF points
	if (indexr >= 0)
	{
		for (double i = 0; i < PF_real_val.size(); i++)
		{
			for (int j = 0; j < Objs(); j++)
			{
				PF_real << PF_real_val[i][j] << " ";
			}
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}
/**************DTLZ3****************************/
std::vector<double> DTLZ3::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;				//fitness vector - output
	double sum_f = 0;							//temp value for the fitness calculation
	double g_x = 0;								//temp value for the fitness calculation
	short nobj = Objs();

	for (int i = nobj - 1; i<Vars(); i++)
	{
		sum_f += pow((code[i] - 0.5), 2.0) - cos(20 * pi*(code[i] - 0.5));
	}
	g_x = 100 * (sum_f + Vars() - nobj + 1) + 1.0;
	sum_f = g_x;
	for (int j = 0; j<nobj - 1; j++)
	{
		sum_f = sum_f * cos(code[j] * pi / 2.0);
	}
	fitness.push_back(sum_f);

	for (int i = 1; i<nobj; i++)
	{
		sum_f = g_x;
		for (int j = 0; j<nobj - 1 - i; j++)
		{
			sum_f = sum_f * cos(code[j] * pi / 2.0);
		}
		sum_f = sum_f * sin(code[nobj - 1 - i] * pi / 2.0);
		fitness.push_back(sum_f);
	}

										//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}
void DTLZ3::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set boundaries for variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}
	//check if the amount of boudaries is correct

	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 100,-1 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> DTLZ3::Plot_PF(int indexr, int size)				//non continous
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//Get the uniform vector
	PF_real_val = Uniform_Weights_Generate(size, Objs(),1);

	//calculate the real PF points
	if (indexr >= 0)
	{
		for (double i = 0; i < PF_real_val.size(); i++)
		{
			for (int j = 0; j < Objs(); j++)
			{
				PF_real << PF_real_val[i][j] << " ";
			}
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}
/**************DTLZ4****************************/
std::vector<double> DTLZ4::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;			//fitness vector - output
	double sum_f = 0;						//temp values for the fitness calculation
	double g_x = 0;							//temp values for the fitness calculation
	short nobj = Objs();
		
	std::vector < double> temp_code = code;
	//calculate temp values
	for (int i = nobj - 1; i<Vars(); i++)
	{
		sum_f += pow((temp_code[i] - 0.5), 2.0);
	}
	for (int i = 0; i<nobj - 1; i++)
	{
		temp_code[i] = pow((float)temp_code[i], (float)100);
	}
	g_x = 1.0 + sum_f;
	sum_f = g_x;
	for (int j = 0; j<nobj - 1; j++)
	{
		sum_f = sum_f * cos(temp_code[j] * pi / 2.0);
	}
	fitness.push_back(sum_f);

	for (int i = 1; i<nobj; i++)
	{
		sum_f = g_x;
		for (int j = 0; j<nobj - 1 - i; j++)
		{
			sum_f = sum_f * cos(temp_code[j] * pi / 2.0);
		}
		sum_f = sum_f * sin(temp_code[nobj - 1 - i] * pi / 2.0);
		fitness.push_back(sum_f);
	}

															//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}
void DTLZ4::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes


	//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 1,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 2,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> DTLZ4::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//Get the uniform vector
	PF_real_val = Uniform_Weights_Generate(size, Objs(),1);

	//calculate the real PF points
	if (indexr >= 0)
	{
		for (double i = 0; i < PF_real_val.size(); i++)
		{
			for (int j = 0; j < Objs(); j++)
			{
				PF_real << PF_real_val[i][j] << " ";
			}
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}
/**************DTLZ5****************************/
std::vector<double> DTLZ5::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;					//fitness vector - output
	double sum_f = 0;								//temp value for the fitness calculation
	double g_x = 0;									//temp value for the fitness calculation
	short nobj = Objs();

	std::vector<double> x(nobj - 1, 0);

													//calculate the fitness and push it to the fitness vector
	for (int i = nobj - 1; i<Vars(); i++)
	{
		sum_f += pow((code[i] - 0.5), 2.0);
	}
	for (int i = 1; i<nobj - 1; i++)
	{
		x[i] = pi / (4 * (1 + sum_f))*(1 + 2 * sum_f*code[i]);
	}
	g_x = 1.0 + sum_f;
	sum_f = g_x;
	for (int j = 1; j<nobj - 1; j++)
	{
		sum_f = sum_f * cos(x[j]);
	}
	sum_f = sum_f * cos(code[0] * pi / 2.0);
	fitness.push_back(sum_f);

	for (int i = 1; i<nobj; i++)
	{
		sum_f = g_x;
		for (int j = 1; j<nobj - 1 - i; j++)
		{
			sum_f = sum_f * cos(x[j]);
		}
		if (i == nobj - 1)
		{
			sum_f = sum_f * sin(code[0] * pi / 2.0);
		}
		else
		{
			sum_f = sum_f * sin(x[nobj - 1 - i]);
			sum_f = sum_f * cos(code[0] * pi / 2.0);
		}
		fitness.push_back(sum_f);
	}

																	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}
void DTLZ5::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set boundaries for variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 1,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 2.5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> DTLZ5::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;				//ofstream file for the PF saving

										//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output
	short n_obj = Objs();

	//calculate the real PF points
	for (double f_m = 0; f_m <= 1; f_m += 1.0 / (double)(size / 2 - 1))
	{
		//calculate the PF value
		double f = sqrt(1 - pow(f_m, 2));			//PF value
		f /= sqrt(2);

		std::vector<double> temp_vect;		//temporary vector of fitness

		for (int i = 0; i < n_obj; i++)
		{
			double temp;
			if (i != n_obj - 1)
				temp = f;
			else
				temp = f_m;

			//push values to the fitness vector
			temp_vect.push_back(temp);

			if (indexr >= 0)
				PF_real << temp << " ";
		}
		if (indexr >= 0)
			PF_real << std::endl;

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

	}

	for (double f = 0; f <= 1; f += 1.0 / (double)(size / 2 - 1))
	{
		//calculate the PF value
		double f_m = sqrt(1 - pow(f, 2));			//PF value

		double f_val = f / sqrt(2);

		std::vector<double> temp_vect;		//temporary vector of fitness

		for (int i = 0; i < n_obj; i++)
		{
			double temp;
			if (i != n_obj - 1)
				temp = f_val;
			else
				temp = f_m;

			//push values to the fitness vector
			temp_vect.push_back(temp);

			if (indexr >= 0)
				PF_real << temp << " ";
		}
		if (indexr >= 0)
			PF_real << std::endl;

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

	}


	//close the PF file
	PF_real.close();

	return PF_real_val;
}
/**************DTLZ6****************************/
std::vector<double> DTLZ6::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;				//fitness vector - output
	double sum_f = 0;							//temp value for the fitness calculation
	double g_x = 0;								//temp value for the fitness calculation
	short nobj = Objs();
	std::vector<double> x(nobj - 1, 0);
												//calculate temp value
	for (int i = nobj - 1; i<Vars(); i++)
	{
		sum_f += pow((code[i]), 0.1);
	}
	for (int i = 1; i<nobj - 1; i++)
	{
		x[i] = pi / (4 * (1 + sum_f))*(1 + 2 * sum_f*code[i]);
	}
	g_x = 1.0 + sum_f;
	sum_f = g_x;
	for (int j = 1; j<nobj - 1; j++)
	{
		sum_f = sum_f * cos(x[j]);
	}
	sum_f = sum_f * cos(code[0] * pi / 2.0);
	fitness.push_back(sum_f);

	for (int i = 1; i<nobj; i++)
	{
		sum_f = g_x;
		for (int j = 1; j<nobj - 1 - i; j++)
		{
			sum_f = sum_f * cos(x[j]);
		}
		if (i == nobj - 1)
		{
			sum_f = sum_f * sin(code[0] * pi / 2.0);
		}
		else
		{
			sum_f = sum_f * sin(x[nobj - 1 - i]);
			sum_f = sum_f * cos(code[0] * pi / 2.0);
		}
		fitness.push_back(sum_f);
	}

																//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}
void DTLZ6::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set boundaries for variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 1,0 };								// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 15,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> DTLZ6::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;				//ofstream file for the PF saving

										//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output
	short n_obj = Objs();

	//calculate the real PF points
	for (double f_m = 0; f_m <= 1; f_m += 1.0 / (double)(size / 2 - 1))
	{
		//calculate the PF value
		double f = sqrt(1 - pow(f_m, 2));			//PF value
		f /= sqrt(2);

		std::vector<double> temp_vect;		//temporary vector of fitness

		for (int i = 0; i < n_obj; i++)
		{
			double temp;
			if (i != n_obj - 1)
				temp = f;
			else
				temp = f_m;

			//push values to the fitness vector
			temp_vect.push_back(temp);

			if (indexr >= 0)
				PF_real << temp << " ";
		}
		if (indexr >= 0)
			PF_real << std::endl;

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

	}

	for (double f = 0; f <= 1; f += 1.0 / (double)(size / 2 - 1))
	{
		//calculate the PF value
		double f_m = sqrt(1 - pow(f, 2));			//PF value

		double f_val = f / sqrt(2);

		std::vector<double> temp_vect;		//temporary vector of fitness

		for (int i = 0; i < n_obj; i++)
		{
			double temp;
			if (i != n_obj - 1)
				temp = f_val;
			else
				temp = f_m;

			//push values to the fitness vector
			temp_vect.push_back(temp);

			if (indexr >= 0)
				PF_real << temp << " ";
		}
		if (indexr >= 0)
			PF_real << std::endl;

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

	}


	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**************DTLZ7****************************/
std::vector<double> DTLZ7::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;				//fitness vector - output
	double sum_f = 0;							//temp value for the fitness calculation
	double f1, f2;								//temp values of fitness 
	double g_x = 0;								//temp value for the fitness calculation
	short nobj = Objs();
	double temp = 0;
												//calculate temp value
	for (int i = nobj - 1; i<Vars(); i++)
	{
		sum_f += code[i];
	}
	g_x = 1.0 + 9.0*sum_f / (Vars() - nobj + 1.0);
	sum_f = g_x;
	for (int i = 0; i<nobj - 1; i++)
	{
		fitness.push_back(code[i]);
	}
	for (int i = 0; i<nobj - 1; i++)
	{
		temp += (fitness[i] / (sum_f + 1))*(1 + sin(3 * pi*fitness[i]));
	}
	temp = nobj - temp;
	fitness.push_back((sum_f + 1)*temp);

																//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}
void DTLZ7::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set boundaries for variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 1,0 };								// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> DTLZ7::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;			//ofstream file for the PF saving

									//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output
	short n_obj = Objs();
	//Get the distribution of points
	std::vector<std::vector<double>> Var_val;

	double step = 1. / (floor(sqrt(size*2)) - 1);

	for (double i = 0; i <= 1; i += step)
	{
		for (double j = 0; j < 1; j += step)
		{
			std::vector<double>temp_var{ i, j };
			Var_val.push_back(temp_var);
		}

	}

	//calculate the real PF points
	for (double i = 0; i < Var_val.size(); i++)
	{
		std::vector<double> temp_var = Var_val[i];
		for (int j = 0; j < 20; j++)
			temp_var.push_back(0.);

		//Calculate the fitness
		std::vector<double> fit = this->Fitness_C(temp_var);


		//push the fitness vector to the PF vector
		PF_real_val.push_back(fit);
	}

	//remove dominated points
	for (int i = 0; i < PF_real_val.size(); i++)
	{
		for (int j = i + 1; j < PF_real_val.size(); j++)
		{
			short flag1 = 0, flag2 = 0;
			for (int k = 0; k < n_obj; k++)
			{
				//Copy ith fitness of each individual
				double fit_ind1;
				double fit_ind2;
				{
					fit_ind1 = PF_real_val[i][k];		//ith fitness of 1st individual
					fit_ind2 = PF_real_val[j][k];		//ith fitness of 2nd individual
				}
				//Check which fitness is greater
				if (fit_ind1 < fit_ind2)
					flag1 = 1;
				else if (fit_ind1 > fit_ind2)
					flag2 = 1;
			}
			//Check which individual dominate
			if (flag1 == 1 && flag2 == 0)
			{
				PF_real_val.erase(PF_real_val.begin() + j);
				j--;
			}
			else if (flag1 == 0 && flag2 == 1)
			{
				PF_real_val.erase(PF_real_val.begin() + i);
				i--;
				break;
			}
			else
				continue;
		}
	}


	if (indexr >= 0)
	{
		for (int i = 0; i < PF_real_val.size(); i++)
		{
			for (int j = 0; j < n_obj; j++)
				PF_real << PF_real_val[i][j] << " ";
			PF_real << std::endl;
		}
	}



	//close the PF file
	PF_real.close();

	return PF_real_val;
}


/**************MOP1****************************/
std::vector<double> MOP1::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;			//fitness vector - output
	double sum_f2 = 0;						//temp value for the fitness calculation
	double g_x = 0;							//temp values for the fitness calculation
	double f1, f2;		//Temporary fitness value
					//calculate temp values

	for (int i = 1; i < Vars(); i++)
	{
		double t;
		t = code[i] - sin(0.5*pi*code[0]);
		sum_f2 += -0.9*t*t + pow(abs(t), 0.6);
		//sum_f2 += sin(1.4*PI*pow(abs(t),0.5))+ pow(abs(t), 0.5);
	}
	//g_x = 1 + 2 * sin(PI*x_var[0])*sum_f2;  - from original paper (check)
	g_x = 1. + 2.* abs(sin(2 * pi*code[0]))*sum_f2;

	//calculate the fitness and push it to the fitness vector
	f1 = g_x*code[0];
	f2 = g_x*(1. - sqrt(code[0]));			//From paper
													//f2 = 1 - pow((f1 / g_x), 2);					//From website - ETH - wrong!!
	fitness.push_back(f1);										//fitness 1
	fitness.push_back(f2);										//fitness 2

																	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}
void MOP1::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set boundaries for variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 1,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> MOP1::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1 - sqrt(i));					//PF value
		std::vector<double> temp_vect;				//temporary vector of fitness

													//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}
/**************MOP2****************************/
std::vector<double> MOP2::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;				//fitness vector - output
	double sum_f2 = 0;							//temp value for the fitness calculation
	double g_x = 0;								//temp value for the fitness calculation
	double f1, f2;		//Temporary fitness value
	
	//calculate temp values
	for (int i = 1; i < Vars(); i++)
	{
		double t = code[i] - sin(0.5*pi*code[0]);
		sum_f2 += abs(t) / (1 + exp(5.*abs(t)));
	}
	g_x = 1. + 10. * sin(pi*code[0])*sum_f2;


	//calculate the fitness and push it to the fitness vector
	f1 = g_x*code[0];
	f2 = g_x*(1. - pow(code[0], 2));			//From paper
													//f2 = 1 - pow((f1 / g_x), 2);					//From website - ETH - wrong!!
	fitness.push_back(f1);										//fitness 1
	fitness.push_back(f2);										//fitness 2

																		//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}
void MOP2::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set boundaries for variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}
	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 1,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> MOP2::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;			//ofstream file for the PF saving

									//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1 - pow(i, 2));		//PF value
		std::vector<double> temp_vect;		//temporary vector of fitness

											//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}
/**************MOP3****************************/
std::vector<double> MOP3::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;				//fitness vector - output
	double sum_f2 = 0;							//temp value for the fitness calculation
	double g_x = 0;								//temp value for the fitness calculation
	double f1, f2;		//Temporary fitness value

	//calculate temp values
	for (int i = 1; i < Vars(); i++)
	{
		double t = code[i] - sin(0.5*pi*code[0]);
		sum_f2 += abs(t) / (1 + exp(5.*abs(t)));
	}
	g_x = 1. + 10. * sin(0.5*pi*code[0])*sum_f2;

	//calculate the fitness and push it to the fitness vector
	f1 = g_x*cos(0.5*pi*code[0]);
	f2 = g_x*sin(0.5*pi*code[0]);			//From paper
													//f2 = 1 - pow((f1 / g_x), 2);					//From website - ETH - wrong!!
	fitness.push_back(f1);										//fitness 1
	fitness.push_back(f2);										//fitness 2

										//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}
void MOP3::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set boundaries for variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}
	//check if the amount of boudaries is correct

	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 1,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> MOP3::Plot_PF(int indexr, int size)				//non continous
{
	std::ofstream PF_real;			//ofstream file for the PF saving

									//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = sqrt(1 - pow(i, 2));		//PF value
		std::vector<double> temp_vect;		//temporary vector of fitness

											//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}
/**************MOP4****************************/
std::vector<double> MOP4::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;				//fitness vector - output
	double sum_f2 = 0;							//temp value for the fitness calculation
	double g_x = 0;									//temp value for the fitness calculation
	double f1, f2;		//Temporary fitness value

	//calculate temp values
	for (int i = 1; i < Vars(); i++)
	{
		double t = code[i] - sin(0.5*pi*code[0]);
		sum_f2 += abs(t) / (1 + exp(5.*abs(t)));
	}
	g_x = 1. + 10. * sin(pi*code[0])*sum_f2;

	//calculate the fitness and push it to the fitness vector
	f1 = g_x*code[0];
	f2 = g_x*(1. - sqrt(code[0])*pow(cos(2.*pi*code[0]), 2));			
													//f2 = 1 - pow((f1 / g_x), 2);					//From website - ETH - wrong!!
	fitness.push_back(f1);										//fitness 1
	fitness.push_back(f2);										//fitness 2

															//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}
void MOP4::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 1,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> MOP4::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

	//calculate the real PF points
	double temp_min = 300; //storage of the previousely min value
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = 1-sqrt(i)*pow(cos(2*pi*i),2);			//PF value

		if (temp < temp_min)
			temp_min = temp;
		else
			continue;

		std::vector<double> temp_vect;		//temporary vector of fitness

											//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}
/**************MOP5****************************/
std::vector<double> MOP5::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;				//fitness vector - output
	double sum_f2 = 0;							//temp value for the fitness calculation
	double g_x = 0;								//temp value for the fitness calculation
	double f1, f2;			//Temporary fitness value

	//calculate temp values
	for (int i = 1; i < Vars(); i++)
	{
		double t = code[i] - sin(0.5*pi*code[0]);
		sum_f2 += -0.9*pow(t,2) + pow(abs(t), 0.6);
	}
	g_x = 1 + 2 * abs(cos(pi*code[0]))*sum_f2;

	//calculate the fitness and push it to the fitness vector
	f1 = g_x*code[0];
	f2 = g_x*(1 - sqrt(code[0]));			//From paper
													//f2 = 1 - pow((f1 / g_x), 2);					//From website - ETH - wrong!!
	fitness.push_back(f1);										//fitness 1
	fitness.push_back(f2);										//fitness 2

																	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}
void MOP5::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set boundaries for variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 1,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 8,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> MOP5::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;				//ofstream file for the PF saving

										//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1 - sqrt(i));			//PF value
		std::vector<double> temp_vect;		//temporary vector of fitness

											//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**************MOP6****************************/
std::vector<double> MOP6::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;				//fitness vector - output
	double sum_f2 = 0;							//temp value for the fitness calculation
	double g_x = 0;								//temp value for the fitness calculation
	double f1, f2, f3;			//Temporary fitness value

	//calculate temp values
	for (int i = 2; i < Vars(); i++)
	{
		double t = code[i] - code[0] * code[1];
		sum_f2 += -0.9*pow(t,2) + pow(abs(t), 0.6);
	}
	g_x = 1 + 2 * sin(pi*code[0])*sum_f2;

	//calculate the fitness and push it to the fitness vector
	f1 = g_x*code[0] * code[1];
	f2 = g_x*code[0] * (1 - code[1]);
	f3 = g_x*(1 - code[0]);

													//f2 = 1 - pow((f1 / g_x), 2);					//From website - ETH - wrong!!
	fitness.push_back(f1);										//fitness 1
	fitness.push_back(f2);										//fitness 2
	fitness.push_back(f3);										//fitness 3
	
	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}
void MOP6::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set boundaries for variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 1,0 };								// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> MOP6::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//Get the uniform vector
	PF_real_val = Uniform_Weights_Generate(size, Objs(), 3);

	//calculate the real PF points
	for (double i = 0; i < PF_real_val.size(); i++)
	{
		for (int j = 0; j < Objs(); j++)
		{
			if (indexr >= 0)
				PF_real << PF_real_val[i][j] << " ";
		}
		if (indexr >= 0)
			PF_real << std::endl;
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}


/**************MOP7****************************/
std::vector<double> MOP7::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;				//fitness vector - output
	double sum_f2 = 0;							//temp value for the fitness calculation
	double g_x = 0;								//temp value for the fitness calculation
	double f1, f2, f3;			//Temporary fitness value

	//calculate temp values
	for (int i = 2; i < Vars(); i++)
	{
		double t = code[i] - code[0] * code[1];
		sum_f2 += -0.9*pow(t,2) + pow(abs(t), 0.6);
	}
	g_x = 1 + 2 * sin(pi*code[0])*sum_f2;

	//calculate the fitness and push it to the fitness vector
	f1 = g_x*cos(0.5*pi*code[0])*cos(0.5*pi*code[1]);
	f2 = g_x*cos(0.5*pi*code[0])*sin(0.5*pi*code[1]);
	f3 = g_x*sin(0.5*pi*code[0]);

	fitness.push_back(f1);										//fitness 1
	fitness.push_back(f2);										//fitness 2
	fitness.push_back(f3);										//fitness 3

																//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}
void MOP7::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set boundaries for variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 1,0 };								// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> MOP7::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//Get the uniform vector
	PF_real_val = Uniform_Weights_Generate(size, Objs(), 1);

	//calculate the real PF points
	for (double i = 0; i < PF_real_val.size(); i++)
	{
		for (int j = 0; j < Objs(); j++)
		{
			if (indexr >= 0)
				PF_real << PF_real_val[i][j] << " ";
		}
		if (indexr >= 0)
			PF_real << std::endl;
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************UF1**************************/
std::vector<double> UF1::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;						//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;

	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - sin((6. * pi * code[0]) + ((double)i * pi) / (double)Vars());
		if (i % 2 == 1)
		{
			sizeJ1++;
			SumJ1 += pow(y, 2);
		}
		else
		{
			sizeJ2++;
			SumJ2 += pow(y, 2);
		}
	}

	std::vector<double> fitness;									//fitness vector - output

	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + (2.0 / sizeJ1)*SumJ1);				//fitness 1
	fitness.push_back(1. - sqrt(code[0]) + (2.0 / sizeJ2)*SumJ2);	//fitness 2
	
	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;

}

void UF1::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes
	
	//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };				
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}
	
	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 4,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 4,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> UF1::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	
	//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");
	
	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output
	
	//calculate the real PF points
	for (double i = 0; i <= 1.; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1. - sqrt(i));				//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness
		
		//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}
	
	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/*****************UF2**********************/

std::vector<double> UF2::Fitness_C(const std::vector<double> & code)
{
	double SumJ1, SumJ2;			//temp values for the fitness calculation
	int sizeJ1, sizeJ2;
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;

	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += pow(code[i-1] - (0.3*pow(code[0], 2)*cos(24 * pi*code[0] + 4. * i*pi / (double)Vars()) + 0.6*code[0])*cos(6 * pi*code[0] + i*pi / (double)Vars()), 2);
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += pow(code[i-1] - (0.3*pow(code[0], 2)*cos(24 * pi*code[0] + 4. * i*pi / (double)Vars()) + 0.6*code[0])*sin(6 * pi*code[0] + i*pi / (double)Vars()), 2);
		}
	}

	std::vector<double> fitness;				//fitness vector - output
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + (2. / sizeJ1)*SumJ1);				//fitness 1
	fitness.push_back(1. - sqrt(code[0]) + (2. / sizeJ2)*SumJ2);		//fitness 2
	
	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;

}

void UF2::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes
	
	//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);
	
	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}
	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 3,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 2,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> UF2::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	
	//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");
	
	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

	//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1. - sqrt(i));			//PF value
		std::vector<double> temp_vect;		//temporary vector of fitness

		//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}
	
	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/*****************UF3**********************/
std::vector<double> UF3::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2, MultipJ1, MultipJ2;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;
	MultipJ1 = MultipJ2 = 1;
	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i-1] - pow(code[0],(0.5*(1.0+3.0*(i-2)/(double)(Vars()-2.))));
		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += pow(y, 2);
			MultipJ1 *= cos((20 * y*pi) / sqrt(i));
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += pow(y, 2);
			MultipJ2 *= cos((20 * y*pi) / sqrt(i));
		}
	}

	std::vector<double> fitness;					//fitness vector - output
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + (2.0 / sizeJ1)*(4.0 * SumJ1 - 2.0 * MultipJ1 + 2.));				//fitness 1
	fitness.push_back(1.0 - sqrt(code[0]) + (2.0 / sizeJ2)*(4.0 * SumJ2 - 2.0 * MultipJ2 + 2.));		//fitness 2
	
	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;

}

void UF3::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes
	
	//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 6,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 6,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> UF3::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;			//ofstream file for the PF saving

	//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");
	
	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

	//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1.0 - sqrt(i));		//PF value
		std::vector<double> temp_vect;		//temporary vector of fitness

		//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/*****************UF4**********************/

std::vector<double> UF4::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;				//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;	

	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i-1] - sin(6 * pi*code[0] + (i*pi) / (double)Vars());
		double h = abs(y) / (1.0 + exp(2 * abs(y)));
		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += h;
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += h;
		}
	}

	std::vector<double> fitness;				//fitness vector - output
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + (2. / sizeJ1)*SumJ1);					//fitness 1
	fitness.push_back(1. - pow(code[0], 2) + (2. / sizeJ2)*SumJ2);		//fitness 2
	
	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;

}

void UF4::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes
	
	//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -2;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 1.5,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 1.5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> UF4::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving
	
	//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");
	
	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output
	
	//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1. - pow(i, 2));		//PF value
		std::vector<double> temp_vect;		//temporary vector of fitness

		//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}
	
	//close the PF file
	PF_real.close();
	return PF_real_val;
}

/*****************UF5**********************/
std::vector<double> UF5::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;
	double eps = 0.1;	// CEC '09
	short N = 10;		// CEC '09
	
	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i-1] - sin(6. * pi*code[0] + i*pi / (double)Vars()); //temp value for the fitness calculation
		double h = 2 * pow(y, 2) - cos(4. * pi*y) + 1.;	//temp value for the fitness calculation
		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += h;
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += h;
		}
	}

	std::vector<double> fitness;						//fitness vector - output
	//calculate the fitness and push it to the fitness vector
	double temp = (1.0 / (2.0 * N) + eps)*abs(sin(2. * N*pi*code[0]));	//temp value for the fitness calculation
	fitness.push_back(code[0] + temp + (2. / sizeJ1)*SumJ1);				//fitness 1
	fitness.push_back(1. - code[0] + temp + (2. / sizeJ2)*SumJ2);			//fitness 2
	
	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;

}

void UF5::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes
	
	//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);
	
	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> UF5::Plot_PF(int indexr, int size)
{
	short N = 10;		//Fucntion parameter - CEC '09
	std::ofstream PF_real;			//ofstream file for the PF saving
	
									//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");
	
	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output
	
	//calculate the real PF points
	for (double i = 0; i <= 2 * N; i++)
	{
		//calculate the PF value
		double temp = i / (2.0 * N);		//PF value

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << temp << " ";
			PF_real << (1. - temp);
			PF_real << std::endl;
		}

		std::vector<double> temp_vect;		//temporary vector of fitness

		//push values to the fitness vector
		temp_vect.push_back(temp);
		temp_vect.push_back(1. - temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

	}

	//close the PF file
	PF_real.close();
	
	return PF_real_val;
}

/*****************UF6**********************/
std::vector<double> UF6::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2, MultipJ1, MultipJ2;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0.;
	MultipJ1 = MultipJ2 = 1.;
	double eps = 0.1;	// CEC '09
	double N = 2;		// CEC '09
	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - sin((6.*pi*code[0])+(i*(pi/(double)Vars())));
		if (i % 2 == 1)
		{
			sizeJ1++;
			SumJ1 += y*y;
			MultipJ1 *= cos((20. * y*pi) / sqrt((double)i));
		}
		else
		{
			sizeJ2++;
			SumJ2 += y*y;
			MultipJ2 *= cos((20. * y*pi) / sqrt((double)i));
		}
	}

	std::vector<double> fitness;					//fitness vector - output
	double max, temp_max;							//value for max function

	//calculate max function value
	temp_max = 2.*( 0.5/ N + eps)*sin(2.*N*pi*code[0]);

	//calculate maximum
	if (temp_max >= 0.)
		max = temp_max;
	else
		max = 0;

	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + max + (2.0 / sizeJ1)*(4.0 * SumJ1 - 2.0 * MultipJ1 + 2.));				//fitness 1
	fitness.push_back(1.0 - code[0] + max + (2.0 / sizeJ2)*(4.0 * SumJ2 - 2.0 * MultipJ2 + 2.));		//fitness 2

																									//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;

}

void UF6::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1.;
		temp.lower = -1.;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> UF6::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;			//ofstream file for the PF saving

									//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

	//calculate the real PF points
	//Push 1st point
	std::vector<double> temp_vect;		//temporary vector of fitness

	//push values to the fitness vector
	temp_vect.push_back(0.);
	temp_vect.push_back(1.);
	//push the fitness vector to the PF vector
	PF_real_val.push_back(temp_vect);

	//save the fitness values to the file
	if (indexr >= 0)
	{
		PF_real << 0. << " ";
		PF_real << 1.;
		PF_real << std::endl;
	}

	//calculate the step of the i value for the loop
	double i_step = (0.5 / ((double)size - 2.0));

	//Calculate the other points
	for (double i = 0.25; i <= 0.5; i += i_step)
	{
		//calculate the PF value
		double temp = (1.0 - i);		//PF value
		std::vector<double> temp_vect;		//temporary vector of fitness

											//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}
	for (double i = 0.75; i <= 1; i += i_step)
	{
		//calculate the PF value
		double temp = (1.0 - i);		//PF value
		std::vector<double> temp_vect;		//temporary vector of fitness

											//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}
/*****************UF7**********************/
std::vector<double> UF7::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;

	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i-1] - sin(6. * pi*code[0] + i*pi / (double)Vars());	//temp value for the fitness calculation
		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += pow(y,2);
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += pow(y, 2);
		}
	}

	std::vector<double> fitness;			//fitness vector - output

	//calculate the fitness and push it to the fitness vector
	fitness.push_back(pow(code[0],0.2) + (2. / sizeJ1)*SumJ1);				//fitness 1
	fitness.push_back(1. - pow(code[0],0.2) + (2. / sizeJ2)*SumJ2);		//fitness 2
	
	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;

}

void UF7::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes
	
	//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);
	
	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 4.5,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 4.5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}

std::vector<std::vector<double>> UF7::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving
	
	//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");
	
	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

	//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1. - i);					//PF value
		
		std::vector<double> temp_vect;			//temporary vector of fitness

		//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl; 
		}
	}
	
	//close the PF file
	PF_real.close();

	return PF_real_val;
}
/**********************UF8**************************/
std::vector<double> UF8::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, sizeJ3, SumJ1, SumJ2, SumJ3;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = sizeJ3 = SumJ1 = SumJ2 = SumJ3 = 0;
	double f1, f2, f3;
	std::vector<double> fitness;					//fitness vector - output

													//calculate temp values
	for (int i = 3; i <= Vars(); i++)
	{
		double y = code[i - 1] - 2 * code[1] * sin(2.*pi*code[0] + ((double)i*pi / (double)Vars()));
		if (i % 3 == 1)
		{
			sizeJ1++;
			SumJ1 += pow(y, 2);
		}
		else if (i % 3 == 2)
		{
			sizeJ2++;
			SumJ2 += pow(y, 2);
		}
		else
		{
			sizeJ3++;
			SumJ3 += pow(y, 2);
		}
	}

	f1 = cos(0.5*code[0] * pi)*cos(0.5*code[1] * pi) + (2. / sizeJ1)*SumJ1;
	f2 = cos(0.5*code[0] * pi)*sin(0.5*code[1] * pi) + (2. / sizeJ2)*SumJ2;
	f3 = sin(0.5*code[0] * pi) + (2. / sizeJ3)*SumJ3;
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(f1);						//fitness 1
	fitness.push_back(f2);		//fitness 2
	fitness.push_back(f3);		//fitness 3
								//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void UF8::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);
	STRUCTURES::boundaries second{ 1,0 };
	bound.push_back(second);

	//set boundaries for the rest of variables
	for (int i = 2; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -2;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 40,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 40,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> UF8::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

	short nobj = Objs();
	//calculate the real PF points
	double step = 1. / (floor(sqrt(size*0.6)) - 1);
	for (double f1 = 0; f1 < 1. + step; f1 += step)
	{
		if (f1 >= 1.)
			f1 = 1.;

		for (double f2 = 0; pow(f2, 2) <= 1. - pow(f1, 2) + step; f2 += step)
		{
			double f3;
			if (pow(f1, 2) + pow(f2, 2) >= 1)
				f3 = 0;
			else
				f3 = sqrt(1. - pow(f1, 2) - pow(f2, 2));

			std::vector<double> temp_vect;			//temporary vector of fitness

													//push values to the fitness vector
			temp_vect.push_back(f1);
			temp_vect.push_back(f2);
			temp_vect.push_back(f3);

			//push the fitness vector to the PF vector
			PF_real_val.push_back(temp_vect);

			if (pow(f2, 2) != 1. - pow(f1, 2) && pow(f2 + step, 2) > 1. - pow(f1, 2))
			{
				double temp_f2 = sqrt(1. - pow(f1, 2));
				f3 = 0;
				std::vector<double> temp_vect;			//temporary vector of fitness

														//push values to the fitness vector
				temp_vect.push_back(f1);
				temp_vect.push_back(temp_f2);
				temp_vect.push_back(f3);

				//push the fitness vector to the PF vector
				PF_real_val.push_back(temp_vect);
			}
		}
	}
	for (double f2 = 0; f2 < 1. + step; f2 += step)
	{

		if (f2 >= 1.)
			f2 = 1.;
		for (double f3 = 0; pow(f3, 2) <= 1. - pow(f2, 2); f3 += step)
		{
			double f1;
			if (pow(f2, 2) + pow(f3, 2) >= 1)
				f1 = 0;
			else
				f1 = sqrt(1. - pow(f2, 2) - pow(f3, 2));

			std::vector<double> temp_vect;			//temporary vector of fitness

													//push values to the fitness vector
			temp_vect.push_back(f1);
			temp_vect.push_back(f2);
			temp_vect.push_back(f3);

			//push the fitness vector to the PF vector
			PF_real_val.push_back(temp_vect);

			if (pow(f3, 2) != 1. - pow(f2, 2) && pow(f3 + step, 2) > 1. - pow(f2, 2))
			{
				double temp_f3 = sqrt(1. - pow(f2, 2));
				f1 = 0;
				std::vector<double> temp_vect;			//temporary vector of fitness

														//push values to the fitness vector
				temp_vect.push_back(f1);
				temp_vect.push_back(f2);
				temp_vect.push_back(temp_f3);

				//push the fitness vector to the PF vector
				PF_real_val.push_back(temp_vect);
			}
		}
	}
	for (double f3 = 0; f3 < 1. + step; f3 += step)
	{
		if (f3 >= 1.)
			f3 = 1.;
		for (double f1 = 0; pow(f1, 2) <= 1. - pow(f3, 2); f1 += step)
		{
			double f2;
			if (pow(f3, 2) + pow(f1, 2) >= 1)
				f2 = 0;
			else
				f2 = sqrt(1. - pow(f3, 2) - pow(f1, 2));

			std::vector<double> temp_vect;			//temporary vector of fitness

													//push values to the fitness vector
			temp_vect.push_back(f1);
			temp_vect.push_back(f2);
			temp_vect.push_back(f3);

			//push the fitness vector to the PF vector
			PF_real_val.push_back(temp_vect);

			if (pow(f1, 2) != 1. - pow(f3, 2) && pow(f1 + step, 2) > 1. - pow(f3, 2))
			{
				double temp_f1 = sqrt(1. - pow(f3, 2));
				f2 = 0;
				std::vector<double> temp_vect;			//temporary vector of fitness

														//push values to the fitness vector
				temp_vect.push_back(temp_f1);
				temp_vect.push_back(f2);
				temp_vect.push_back(f3);

				//push the fitness vector to the PF vector
				PF_real_val.push_back(temp_vect);
			}
		}
	}
	//refine PF
	for (int i = 0; i < PF_real_val.size(); i++)
	{
		for (int j = i + 1; j < PF_real_val.size(); j++)
		{
			int eq = 0;
			for (int k = 0; k < nobj; k++)
			{
				if (abs(PF_real_val[i][k] - PF_real_val[j][k]) <= pow(10, -pow_eq_zero))
					eq++;
				else
					break;
			}
			if (eq == nobj)
			{
				PF_real_val.erase(PF_real_val.begin() + j);
				j--;
				continue;
			}
		}
		if (indexr >= 0)
		{
			for (int k = 0; k < nobj; k++)
			{
				PF_real << PF_real_val[i][k] << " ";
			}
			PF_real << std::endl;
		}
	}


	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************UF9**************************/
std::vector<double> UF9::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, sizeJ3, SumJ1, SumJ2, SumJ3;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = sizeJ3 = SumJ1 = SumJ2 = SumJ3 = 0;
	double f1, f2, f3;
	std::vector<double> fitness;					//fitness vector - output
	float eps = 0.1;
													//calculate temp values
	for (int i = 3; i <= Vars(); i++)
	{
		double y = code[i - 1] - 2 * code[1] * sin(2.*pi*code[0] + ((double)i*pi / (double)Vars()));
		if (i % 3 == 1)
		{
			sizeJ1++;
			SumJ1 += pow(y, 2);
		}
		else if (i % 3 == 2)
		{
			sizeJ2++;
			SumJ2 += pow(y, 2);
		}
		else
		{
			sizeJ3++;
			SumJ3 += pow(y, 2);
		}
	}
	double max = (1. + eps)*(1. - 4.*pow(2 * code[0] - 1, 2));
	if (max < 0)
		max = 0;
	f1 = 0.5*(max + 2 * code[0])*code[1] + (2. / sizeJ1)*SumJ1;
	f2 = 0.5*(max - 2 * code[0] + 2)*code[1] + +(2. / sizeJ2)*SumJ2;
	f3 = 1 - code[1] + (2. / sizeJ3)*SumJ3;
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(f1);						//fitness 1
	fitness.push_back(f2);		//fitness 2
	fitness.push_back(f3);		//fitness 3
								//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void UF9::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);
	STRUCTURES::boundaries second{ 1,0 };
	bound.push_back(second);

	//set boundaries for the rest of variables
	for (int i = 2; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -2;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 40,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 40,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> UF9::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

															//calculate the real PF points
	double step = 1. / (floor(sqrt(size * 5)) - 1);
	for (double f3 = 0; f3 <= 1.; f3 += step)
	{
		for (double f1 = 0; f1 <= 1. - f3; f1 += step)
		{
			if (f1 > 0.25*(1 - f3) && f1 < 0.75*(1 - f3))
				continue;
			double f2 = 1 - f1 - f3;

			std::vector<double> temp_vect;			//temporary vector of fitness

													//push values to the fitness vector
			temp_vect.push_back(f1);
			temp_vect.push_back(f2);
			temp_vect.push_back(f3);

			//push the fitness vector to the PF vector
			PF_real_val.push_back(temp_vect);

			//save the fitness values to the file
			if (indexr >= 0)
			{
				PF_real << f1 << " ";
				PF_real << f2 << " ";
				PF_real << f3;
				PF_real << std::endl;
			}
		}
	}
	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************UF10**************************/
std::vector<double> UF10::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, sizeJ3, SumJ1, SumJ2, SumJ3;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = sizeJ3 = SumJ1 = SumJ2 = SumJ3 = 0;
	double f1, f2, f3;
	std::vector<double> fitness;					//fitness vector - output

													//calculate temp values
	for (int i = 3; i <= Vars(); i++)
	{
		double y = code[i - 1] - 2 * code[1] * sin(2.*pi*code[0] + ((double)i*pi / (double)Vars()));
		if (i % 3 == 1)
		{
			sizeJ1++;
			SumJ1 += 4 * pow(y, 2) - cos(8 * pi*y) + 1.;
		}
		else if (i % 3 == 2)
		{
			sizeJ2++;
			SumJ2 += 4 * pow(y, 2) - cos(8 * pi*y) + 1.;
		}
		else
		{
			sizeJ3++;
			SumJ3 += 4 * pow(y, 2) - cos(8 * pi*y) + 1.;
		}
	}

	f1 = cos(0.5*code[0] * pi)*cos(0.5*code[1] * pi) + (2. / sizeJ1)*SumJ1;
	f2 = cos(0.5*code[0] * pi)*sin(0.5*code[1] * pi) + (2. / sizeJ2)*SumJ2;
	f3 = sin(0.5*code[0] * pi) + (2. / sizeJ3)*SumJ3;
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(f1);						//fitness 1
	fitness.push_back(f2);		//fitness 2
	fitness.push_back(f3);		//fitness 3
								//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void UF10::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);
	STRUCTURES::boundaries second{ 1,0 };
	bound.push_back(second);

	//set boundaries for the rest of variables
	for (int i = 2; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -2;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 40,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 40,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> UF10::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

	short nobj = Objs();
	//calculate the real PF points
	double step = 1. / (floor(sqrt(size*0.6)) - 1);
	for (double f1 = 0; f1 < 1. + step; f1 += step)
	{
		if (f1 >= 1.)
			f1 = 1.;

		for (double f2 = 0; pow(f2, 2) <= 1. - pow(f1, 2) + step; f2 += step)
		{
			double f3;
			if (pow(f1, 2) + pow(f2, 2) >= 1)
				f3 = 0;
			else
				f3 = sqrt(1. - pow(f1, 2) - pow(f2, 2));

			std::vector<double> temp_vect;			//temporary vector of fitness

													//push values to the fitness vector
			temp_vect.push_back(f1);
			temp_vect.push_back(f2);
			temp_vect.push_back(f3);

			//push the fitness vector to the PF vector
			PF_real_val.push_back(temp_vect);

			if (pow(f2, 2) != 1. - pow(f1, 2) && pow(f2 + step, 2) > 1. - pow(f1, 2))
			{
				double temp_f2 = sqrt(1. - pow(f1, 2));
				f3 = 0;
				std::vector<double> temp_vect;			//temporary vector of fitness

														//push values to the fitness vector
				temp_vect.push_back(f1);
				temp_vect.push_back(temp_f2);
				temp_vect.push_back(f3);

				//push the fitness vector to the PF vector
				PF_real_val.push_back(temp_vect);
			}
		}
	}
	for (double f2 = 0; f2 < 1. + step; f2 += step)
	{

		if (f2 >= 1.)
			f2 = 1.;
		for (double f3 = 0; pow(f3, 2) <= 1. - pow(f2, 2); f3 += step)
		{
			double f1;
			if (pow(f2, 2) + pow(f3, 2) >= 1)
				f1 = 0;
			else
				f1 = sqrt(1. - pow(f2, 2) - pow(f3, 2));

			std::vector<double> temp_vect;			//temporary vector of fitness

													//push values to the fitness vector
			temp_vect.push_back(f1);
			temp_vect.push_back(f2);
			temp_vect.push_back(f3);

			//push the fitness vector to the PF vector
			PF_real_val.push_back(temp_vect);

			if (pow(f3, 2) != 1. - pow(f2, 2) && pow(f3 + step, 2) > 1. - pow(f2, 2))
			{
				double temp_f3 = sqrt(1. - pow(f2, 2));
				f1 = 0;
				std::vector<double> temp_vect;			//temporary vector of fitness

														//push values to the fitness vector
				temp_vect.push_back(f1);
				temp_vect.push_back(f2);
				temp_vect.push_back(temp_f3);

				//push the fitness vector to the PF vector
				PF_real_val.push_back(temp_vect);
			}
		}
	}
	for (double f3 = 0; f3 < 1. + step; f3 += step)
	{
		if (f3 >= 1.)
			f3 = 1.;
		for (double f1 = 0; pow(f1, 2) <= 1. - pow(f3, 2); f1 += step)
		{
			double f2;
			if (pow(f3, 2) + pow(f1, 2) >= 1)
				f2 = 0;
			else
				f2 = sqrt(1. - pow(f3, 2) - pow(f1, 2));

			std::vector<double> temp_vect;			//temporary vector of fitness

													//push values to the fitness vector
			temp_vect.push_back(f1);
			temp_vect.push_back(f2);
			temp_vect.push_back(f3);

			//push the fitness vector to the PF vector
			PF_real_val.push_back(temp_vect);

			if (pow(f1, 2) != 1. - pow(f3, 2) && pow(f1 + step, 2) > 1. - pow(f3, 2))
			{
				double temp_f1 = sqrt(1. - pow(f3, 2));
				f2 = 0;
				std::vector<double> temp_vect;			//temporary vector of fitness

														//push values to the fitness vector
				temp_vect.push_back(temp_f1);
				temp_vect.push_back(f2);
				temp_vect.push_back(f3);

				//push the fitness vector to the PF vector
				PF_real_val.push_back(temp_vect);
			}
		}
	}
	//refine PF
	for (int i = 0; i < PF_real_val.size(); i++)
	{
		for (int j = i + 1; j < PF_real_val.size(); j++)
		{
			int eq = 0;
			for (int k = 0; k < nobj; k++)
			{
				if (abs(PF_real_val[i][k] - PF_real_val[j][k]) <= pow(10, -pow_eq_zero))
					eq++;
				else
					break;
			}
			if (eq == nobj)
			{
				PF_real_val.erase(PF_real_val.begin() + j);
				j--;
				continue;
			}
		}
		if (indexr >= 0)
		{
			for (int k = 0; k < nobj; k++)
			{
				PF_real << PF_real_val[i][k] << " ";
			}
			PF_real << std::endl;
		}
	}


	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************WFG1**************************/
std::vector<double> WFG1::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness(Objs(), 0.);					//fitness vector - output

	WFG_Calc(1, code, fitness);
								//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void WFG1::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes


	//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2 * (i + 1);
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 4,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> WFG1::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

	//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output
	
	int nobj = Objs();
	int k_num = 2 * (nobj - 1);
	std::vector<std::vector<double>> decision = Points_Generate(size, k_num);
	std::vector<double> fitness(nobj, 0.);
	for (int j = 0; j < decision.size(); j++)
	{
		std::vector<double> temp_decision = decision[j];

		for (int k = 1; k <= k_num; k++)
		{
			temp_decision[k - 1] *= 2 * k;
		}

		for (int k = k_num + 1; k <= k_num + 20; k++)
		{
			temp_decision.push_back(2 * k*0.35);
		}
		WFG_Calc(1, temp_decision, fitness);

		PF_real_val.push_back(fitness);

		for (int k = 0; k < nobj; k++)
		{
			PF_real << fitness[k] << " ";
		}
		PF_real << std::endl;

	}
	PF_real.close();

	return PF_real_val;
}


/**********************WFG2**************************/
std::vector<double> WFG2::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness(Objs(), 0.);					//fitness vector - output

	WFG_Calc(2, code, fitness);
	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void WFG2::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes


																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2 * (i + 1);
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 3,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> WFG2::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

	int nobj = Objs();
	int k_num = 2 * (nobj - 1);
	std::vector<std::vector<double>> decision = Points_Generate(size, k_num);
	std::vector<double> fitness(nobj, 0.);
	for (int j = 0; j < decision.size(); j++)
	{
		std::vector<double> temp_decision = decision[j];

		for (int k = 1; k <= k_num; k++)
		{
			temp_decision[k - 1] *= 2 * k;
		}

		for (int k = k_num + 1; k <= k_num + 20; k++)
		{
			temp_decision.push_back(2 * k*0.35);
		}
		WFG_Calc(2, temp_decision, fitness);

		PF_real_val.push_back(fitness);
	}


	//remove dominated points
	for (int i = 0; i < PF_real_val.size(); i++)
	{
		for (int j = i + 1; j < PF_real_val.size(); j++)
		{
			short flag1 = 0, flag2 = 0;
			for (int k = 0; k < nobj; k++)
			{
				//Copy ith fitness of each individual
				double fit_ind1;
				double fit_ind2;
				{
					fit_ind1 = PF_real_val[i][k];		//ith fitness of 1st individual
					fit_ind2 = PF_real_val[j][k];		//ith fitness of 2nd individual
				}
				//Check which fitness is greater
				if (fit_ind1 < fit_ind2)
					flag1 = 1;
				else if (fit_ind1 > fit_ind2)
					flag2 = 1;
			}
			//Check which individual dominate
			if (flag1 == 1 && flag2 == 0)
			{
				PF_real_val.erase(PF_real_val.begin() + j);
				j--;
			}
			else if (flag1 == 0 && flag2 == 1)
			{
				PF_real_val.erase(PF_real_val.begin() + i);
				i--;
				break;
			}
			else
				continue;
		}
	}


	if (indexr >= 0)
	{
		for (int i = 0; i < PF_real_val.size(); i++)
		{
			for (int j = 0; j < nobj; j++)
				PF_real << PF_real_val[i][j] << " ";
			PF_real << std::endl;
		}
	}

	//check dominance

	PF_real.close();

	return PF_real_val;
}

/**********************WFG3**************************/
std::vector<double> WFG3::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness(Objs(), 0.);					//fitness vector - output

	WFG_Calc(3, code, fitness);
	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void WFG3::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes


																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2 * (i + 1);
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 3,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> WFG3::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

	int nobj = Objs();
	int k_num = 2 * (nobj - 1);
	std::vector<std::vector<double>> decision = Points_Generate(size, k_num);
	std::vector<double> fitness(nobj, 0.);
	for (int j = 0; j < decision.size(); j++)
	{
		std::vector<double> temp_decision = decision[j];

		for (int k = 1; k <= k_num; k++)
		{
			temp_decision[k - 1] *= 2 * k;
		}

		for (int k = k_num + 1; k <= k_num + 20; k++)
		{
			temp_decision.push_back(2 * k*0.35);
		}
		WFG_Calc(3, temp_decision, fitness);

		PF_real_val.push_back(fitness);

		for (int k = 0; k < nobj; k++)
		{
			PF_real << fitness[k] << " ";
		}
		PF_real << std::endl;

	}
	PF_real.close();

	return PF_real_val;
}

/**********************WFG4**************************/
std::vector<double> WFG4::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness(Objs(), 0.);					//fitness vector - output

	WFG_Calc(4, code, fitness);
	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void WFG4::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes


																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2 * (i + 1);
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 3,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> WFG4::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//Get the uniform vector
	PF_real_val = Uniform_Weights_Generate(size, Objs(), 1);

	//calculate the real PF points
	if (indexr >= 0)
	{
		for (double i = 0; i < PF_real_val.size(); i++)
		{
			for (int j = 0; j < Objs(); j++)
			{
				PF_real_val[i][j] *= 2 * (j + 1);
				PF_real << PF_real_val[i][j] << " ";
			}
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************WFG5**************************/
std::vector<double> WFG5::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness(Objs(), 0.);					//fitness vector - output

	WFG_Calc(5, code, fitness);
	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void WFG5::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes


																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2 * (i + 1);
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 3,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> WFG5::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//Get the uniform vector
	PF_real_val = Uniform_Weights_Generate(size, Objs(), 1);

	//calculate the real PF points
	if (indexr >= 0)
	{
		for (double i = 0; i < PF_real_val.size(); i++)
		{
			for (int j = 0; j < Objs(); j++)
			{
				PF_real_val[i][j] *= 2 * (j + 1);
				PF_real << PF_real_val[i][j] << " ";
			}
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************WFG6**************************/
std::vector<double> WFG6::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness(Objs(), 0.);					//fitness vector - output

	WFG_Calc(6, code, fitness);
	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void WFG6::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes


																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2 * (i + 1);
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 3,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> WFG6::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//Get the uniform vector
	PF_real_val = Uniform_Weights_Generate(size, Objs(), 1);

	//calculate the real PF points
	if (indexr >= 0)
	{
		for (double i = 0; i < PF_real_val.size(); i++)
		{
			for (int j = 0; j < Objs(); j++)
			{
				PF_real_val[i][j] *= 2 * (j + 1);
				PF_real << PF_real_val[i][j] << " ";
			}
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************WFG7**************************/
std::vector<double> WFG7::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness(Objs(), 0.);					//fitness vector - output

	WFG_Calc(7, code, fitness);
	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void WFG7::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes


																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2 * (i + 1);
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 3,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> WFG7::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//Get the uniform vector
	PF_real_val = Uniform_Weights_Generate(size, Objs(), 1);

	//calculate the real PF points
	if (indexr >= 0)
	{
		for (double i = 0; i < PF_real_val.size(); i++)
		{
			for (int j = 0; j < Objs(); j++)
			{
				PF_real_val[i][j] *= 2 * (j + 1);
				PF_real << PF_real_val[i][j] << " ";
			}
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************WFG8**************************/
std::vector<double> WFG8::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness(Objs(), 0.);					//fitness vector - output

	WFG_Calc(8, code, fitness);
	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void WFG8::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes


																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2 * (i + 1);
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 3,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> WFG8::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//Get the uniform vector
	PF_real_val = Uniform_Weights_Generate(size, Objs(), 1);

	//calculate the real PF points
	if (indexr >= 0)
	{
		for (double i = 0; i < PF_real_val.size(); i++)
		{
			for (int j = 0; j < Objs(); j++)
			{
				PF_real_val[i][j] *= 2 * (j + 1);
				PF_real << PF_real_val[i][j] << " ";
			}
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************WFG9**************************/
std::vector<double> WFG9::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness(Objs(), 0.);					//fitness vector - output

	WFG_Calc(9, code, fitness);
	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void WFG9::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes


																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2 * (i + 1);
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 3,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> WFG9::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//Get the uniform vector
	PF_real_val = Uniform_Weights_Generate(size, Objs(), 1);

	//calculate the real PF points
	if (indexr >= 0)
	{
		for (double i = 0; i < PF_real_val.size(); i++)
		{
			for (int j = 0; j < Objs(); j++)
			{
				PF_real_val[i][j] *= 2 * (j + 1);
				PF_real << PF_real_val[i][j] << " ";
			}
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}
/**********************IMB1**************************/
std::vector<double> IMB1::Fitness_C(const std::vector<double> & code)
{
	double g_x = 0;			//temp values for the fitness calculation
	double f1, f2, f3;
	std::vector<double> fitness;					//fitness vector - output

	//calculate temp values
	if (code[0] <= 0.2 && code[0] >= 0.)
		;
	else
	{
		for (int i = 2; i <= Vars(); i++)
		{
			double t = code[i - 1] - sin(0.5*pi*code[0]);
			g_x += 0.5*(-0.9*pow(t, 2) + pow(abs(t), 0.6));
		}
	}

	f1 = code[0] * (1 + g_x);
	f2 = (1 - sqrt(code[0])) * (1 + g_x);
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(f1);						//fitness 1
	fitness.push_back(f2);		//fitness 2
								//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void IMB1::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes



	//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 3,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 2,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> IMB1::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = 1. - sqrt(i);					//PF value

		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************IMB2**************************/
std::vector<double> IMB2::Fitness_C(const std::vector<double> & code)
{
	double g_x = 0;			//temp values for the fitness calculation
	double f1, f2, f3;
	std::vector<double> fitness;					//fitness vector - output

													//calculate temp values
	if (code[0] >= 0.4 && code[0] <= 0.6)
		;
	else
	{
		for (int i = 2; i <= Vars(); i++)
		{
			double t = code[i - 1] - sin(0.5*pi*code[0]);
			g_x += 0.5*(-0.9*pow(t, 2) + pow(abs(t), 0.6));
		}
	}

	f1 = code[0] * (1 + g_x);
	f2 = (1 - code[0]) * (1 + g_x);
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(f1);						//fitness 1
	fitness.push_back(f2);		//fitness 2
								//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void IMB2::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes



																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 3,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 3,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> IMB2::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = 1. - i;					//PF value

		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************IMB3**************************/
std::vector<double> IMB3::Fitness_C(const std::vector<double> & code)
{
	double g_x = 0;			//temp values for the fitness calculation
	double f1, f2, f3;
	std::vector<double> fitness;					//fitness vector - output

													//calculate temp values
	if (code[0] >= 0.8 && code[0] <= 1)
		;
	else
	{
		for (int i = 2; i <= Vars(); i++)
		{
			double t = code[i - 1] - sin(0.5*pi*code[0]);
			g_x += 0.5*(-0.9*pow(t, 2) + pow(abs(t), 0.6));
		}
	}

	f1 = cos(0.5*pi*code[0]) * (1 + g_x);
	f2 = sin(0.5*pi*code[0]) * (1 + g_x);
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(f1);						//fitness 1
	fitness.push_back(f2);		//fitness 2
								//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void IMB3::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes



																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 3,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 3,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> IMB3::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = sqrt(1. - pow(i, 2));					//PF value

		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************IMB4**************************/
std::vector<double> IMB4::Fitness_C(const std::vector<double> & code)
{
	double g_x = 0;			//temp values for the fitness calculation
	double f1, f2, f3;
	std::vector<double> fitness;					//fitness vector - output

													//calculate temp values
	if (code[0] >= 2./3. && code[0] <= 1.)
		;
	else
	{
		for (int i = 3; i <= Vars(); i++)
		{
			double t = code[i - 1] - 0.5*(code[0]+code[1]);
			g_x +=-0.9*pow(t, 2) + pow(abs(t), 0.6);
		}
		g_x *= 2 * cos(0.5*pi*code[0]);
	}

	f1 = code[0] * code[1] * (1 + g_x);
	f2 = code[0] * (1 - code[1]) * (1 + g_x);
	f3 = (1 - code[0]) * (1 + g_x);
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(f1);						//fitness 1
	fitness.push_back(f2);						//fitness 2
	fitness.push_back(f3);						//fitness 3
								//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void IMB4::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes



																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 6,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 6,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> IMB4::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

	//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

	//Get the uniform vector
	PF_real_val = Uniform_Weights_Generate(size, Objs(), 3);

	//calculate the real PF points
	for (double i = 0; i < PF_real_val.size(); i++)
	{
		for (int j = 0; j < Objs(); j++)
		{
			if (indexr >= 0)
				PF_real << PF_real_val[i][j] << " ";
		}
		if (indexr >= 0)
			PF_real << std::endl;
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************IMB5**************************/
std::vector<double> IMB5::Fitness_C(const std::vector<double> & code)
{
	double g_x = 0;			//temp values for the fitness calculation
	double f1, f2, f3;
	std::vector<double> fitness;					//fitness vector - output

													//calculate temp values
	if (code[0] >= 0 && code[0] <= 0.5)
		;
	else
	{
		for (int i = 3; i <= Vars(); i++)
		{
			double t = code[i - 1] - 0.5*(code[0] + code[1]);
			g_x += -0.9*pow(t, 2) + pow(abs(t), 0.6);
		}
		g_x *= 2 * cos(0.5*pi*code[0]);
	}

	f1 = cos(0.5*pi*code[0])*cos(0.5*pi*code[1]) * (1 + g_x);
	f2 = cos(0.5*pi*code[0])*sin(0.5*pi*code[1]) * (1 + g_x);
	f3 = sin(0.5*pi*code[0]) * (1 + g_x);
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(f1);						//fitness 1
	fitness.push_back(f2);						//fitness 2
	fitness.push_back(f3);						//fitness 3
												//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void IMB5::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes



																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 40,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 40,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> IMB5::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//Get the uniform vector
	PF_real_val = Uniform_Weights_Generate(size, Objs(), 1);

	//calculate the real PF points
	for (double i = 0; i < PF_real_val.size(); i++)
	{
		for (int j = 0; j < Objs(); j++)
		{
			if (indexr >= 0)
				PF_real << PF_real_val[i][j] << " ";
		}
		if (indexr >= 0)
			PF_real << std::endl;
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************IMB6**************************/
std::vector<double> IMB6::Fitness_C(const std::vector<double> & code)
{
	double g_x = 0;			//temp values for the fitness calculation
	double f1, f2, f3;
	std::vector<double> fitness;					//fitness vector - output

													//calculate temp values
	if (code[0] >= 0 && code[0] <= 0.75)
		;
	else
	{
		for (int i = 3; i <= Vars(); i++)
		{
			double t = code[i - 1] - 0.5*(code[0] + code[1]);
			g_x += -0.9*pow(t, 2) + pow(abs(t), 0.6);
		}
		g_x *= 2 * cos(0.5*pi*code[0]);
	}

	f1 = code[0] * code[1] * (1 + g_x);
	f2 = code[0] * (1 - code[1]) * (1 + g_x);
	f3 = (1 - code[0]) * (1 + g_x);
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(f1);						//fitness 1
	fitness.push_back(f2);						//fitness 2
	fitness.push_back(f3);						//fitness 3
												//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void IMB6::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes



																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 40,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 40,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> IMB6::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//Get the uniform vector
	PF_real_val = Uniform_Weights_Generate(size, Objs(), 3);

	//calculate the real PF points
	for (double i = 0; i < PF_real_val.size(); i++)
	{
		for (int j = 0; j < Objs(); j++)
		{
			if (indexr >= 0)
				PF_real << PF_real_val[i][j] << " ";
		}
		if (indexr >= 0)
			PF_real << std::endl;
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************IMB7**************************/
std::vector<double> IMB7::Fitness_C(const std::vector<double> & code)
{
	double g_x = 0;			//temp values for the fitness calculation
	double f1, f2, f3;
	std::vector<double> fitness;					//fitness vector - output

													//calculate temp values
	if (code[0] <= 0.5 && code[0] >= 0.8)
	{
		for (int i = 2; i <= Vars(); i++)
		{
			double s = code[i - 1] - sin(0.5*pi*code[0]);
			g_x += -0.9*pow(s, 2) + pow(abs(s), 0.6);
		}
	}
	else
	{
		for (int i = 2; i <= Vars(); i++)
		{
			double t = code[i - 1] - 0.5;
			g_x += pow(abs(t), 0.6);
		}
	}

	f1 = code[0] * (1 + g_x);
	f2 = (1 - sqrt(code[0])) * (1 + g_x);
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(f1);						//fitness 1
	fitness.push_back(f2);		//fitness 2
								//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void IMB7::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes



																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 6,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> IMB7::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = 1. - sqrt(i);					//PF value

		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************IMB8**************************/
std::vector<double> IMB8::Fitness_C(const std::vector<double> & code)
{
	double g_x = 0;			//temp values for the fitness calculation
	double f1, f2, f3;
	std::vector<double> fitness;					//fitness vector - output

													//calculate temp values
	if (code[0] <= 0.5 && code[0] >= 0.8)
	{
		for (int i = 2; i <= Vars(); i++)
		{
			double s = code[i - 1] - sin(0.5*pi*code[0]);
			g_x += -0.9*pow(s, 2) + pow(abs(s), 0.6);
		}
	}
	else
	{
		for (int i = 2; i <= Vars(); i++)
		{
			double t = code[i - 1] - 0.5;
			g_x += pow(abs(t), 0.6);
		}
	}

	f1 = code[0] * (1 + g_x);
	f2 = (1 - code[0]) * (1 + g_x);
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(f1);						//fitness 1
	fitness.push_back(f2);		//fitness 2
								//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void IMB8::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes



																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 6,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> IMB8::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = 1. - i;					//PF value

		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************IMB9**************************/
std::vector<double> IMB9::Fitness_C(const std::vector<double> & code)
{
	double g_x = 0;			//temp values for the fitness calculation
	double f1, f2, f3;
	std::vector<double> fitness;					//fitness vector - output

													//calculate temp values
	if (code[0] <= 0.5 && code[0] >= 0.8)
	{
		for (int i = 2; i <= Vars(); i++)
		{
			double s = code[i - 1] - sin(0.5*pi*code[0]);
			g_x += -0.9*pow(s, 2) + pow(abs(s), 0.6);
		}
	}
	else
	{
		for (int i = 2; i <= Vars(); i++)
		{
			double t = code[i - 1] - 0.5;
			g_x += pow(abs(t), 0.6);
		}
	}

	f1 = cos(0.5*pi*code[0]) * (1 + g_x);
	f2 = sin(0.5*pi*code[0]) * (1 + g_x);
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(f1);						//fitness 1
	fitness.push_back(f2);		//fitness 2
								//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void IMB9::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes



																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 6,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> IMB9::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = sqrt(1. - pow(i, 2));					//PF value

		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************IMB10**************************/
std::vector<double> IMB10::Fitness_C(const std::vector<double> & code)
{
	double g_x = 0;			//temp values for the fitness calculation
	double f1, f2, f3;
	std::vector<double> fitness;					//fitness vector - output

													//calculate temp values
	if (code[0] <= 0.5 && code[0] >= 0.8 && code[1] <= 0.5 && code[1] >= 0.8)
	{
		for (int i = 3; i <= Vars(); i++)
		{
			double s = code[i - 1] - 0.5*(code[0] + code[1]);
			g_x += 2*(-0.9*pow(s, 2) + pow(abs(s), 0.6));
		}
	}
	else
	{
		for (int i = 3; i <= Vars(); i++)
		{
			double t = code[i - 1] - code[0] * code[1];
			g_x += pow(abs(t), 0.6);
		}
	}

	f1 = code[0] * code[1] * (1 + g_x);
	f2 = code[0] * (1 - code[1]) * (1 + g_x);
	f3 = (1 - code[0]) * (1 + g_x);
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(f1);						//fitness 1
	fitness.push_back(f2);						//fitness 2
	fitness.push_back(f3);						//fitness 3
												//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void IMB10::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes



																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 40,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 40,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> IMB10::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//Get the uniform vector
	PF_real_val = Uniform_Weights_Generate(size, Objs(), 3);

	//calculate the real PF points
	for (double i = 0; i < PF_real_val.size(); i++)
	{
		for (int j = 0; j < Objs(); j++)
		{
			if (indexr >= 0)
				PF_real << PF_real_val[i][j] << " ";
		}
		if (indexr >= 0)
			PF_real << std::endl;
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}


/********************************************
			CONSTRAINED FUNCTIONS
*********************************************/

/**********************CF1**************************/
std::vector<double> CF1::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;						//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;

	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - pow(code[0],(0.5*(1.+(3.*(i-2)/(double)(Vars()-2)))));
		if (i % 2 == 1)
		{
			sizeJ1++;
			SumJ1 += pow(y, 2);
		}
		else
		{
			sizeJ2++;
			SumJ2 += pow(y, 2);
		}
	}

	std::vector<double> fitness;									//fitness vector - output

																	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + (2.0 / sizeJ1)*SumJ1);				//fitness 1
	fitness.push_back(1. - code[0] + (2.0 / sizeJ2)*SumJ2);	//fitness 2

																	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;

}

void CF1::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 2.5f,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 2.5f,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> CF1::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	short N = 10;		// CEC '09
	
	//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

	//calculate the real PF points
	for (double i = 0; i <= 2*N; i++)
	{
		//calculate the PF value
		double temp = (1. - ((double)i / (double)(2 * N)));				//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back((double)i / (double)(2 * N));
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << (double)i / (double)(2 * N) << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> CF1::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	double cons_value;					//Value for constrain check
	short N = 10;		// CEC '09
	short a = 1;		// CEC '09
	std::vector<double> cons_val_temp;	//Values of the constrains - output
	
	//Calculate constrain
	cons_value = fit[0] + fit[1] - (double)a*abs(sin((double)N*pi*(fit[0] - fit[1] + 1.))) - 1.;

	//push value to the vector
	cons_val_temp.push_back(cons_value);

	//return vector
	return cons_val_temp;

}

/**********************CF2**************************/
std::vector<double> CF2::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;						//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;

	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y; 
		if (i % 2 == 1)
		{
			y = code[i - 1] - sin((6. * pi * code[0]) + ((double)i * pi) / (double)Vars());
			sizeJ1++;
			SumJ1 += pow(y, 2);
		}
		else
		{
			y = code[i - 1] - cos((6. * pi * code[0]) + ((double)i * pi) / (double)Vars());
			sizeJ2++;
			SumJ2 += pow(y, 2);
		}
	}

	std::vector<double> fitness;									//fitness vector - output

																	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + (2.0 / sizeJ1)*SumJ1);				//fitness 1
	fitness.push_back(1. - sqrt(code[0]) + (2.0 / sizeJ2)*SumJ2);	//fitness 2

																	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;

}

void CF2::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 6,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> CF2::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	double N = 2.;		// CEC '09

	//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

	

	//calculate the real PF points
	//Push 1st point
	std::vector<double> temp_vect;		//temporary vector of fitness

										//push values to the fitness vector
	temp_vect.push_back(0.);
	temp_vect.push_back(1.);
	//push the fitness vector to the PF vector
	PF_real_val.push_back(temp_vect);

	//save the fitness values to the file
	if (indexr >= 0)
	{
		PF_real << 0. << " ";
		PF_real << 1.;
		PF_real << std::endl;
	}

	//calculate the step of the i value for the loop
	double i_step = ((pow(2. / 4., 2) - pow(1. / 4., 2)) + (pow(4. / 4., 2) - pow(3. / 4., 2))) / ((double)size - 2.0);

	//Calculate the rest of the points
	for (double i = pow(1. / 4., 2); i <= pow(2. / 4., 2); i += i_step)
	{
		//calculate the PF value
		double temp = (1. - sqrt(i));				//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	for (double i = pow(3. / 4., 2); i <= pow(4. / 4., 2); i += i_step)
	{
		//calculate the PF value
		double temp = (1. - sqrt(i));				//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> CF2::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double time)
{
	double cons_value;					//Value for constrain check
	double t;							//Temporary value
	short N = 2;		// CEC '09
	short a = 1;		// CEC '09
	std::vector<double> cons_val_temp;	//Values of the constrains - output

	//Calculate t
	t = fit[1] + sqrt(fit[0]) - (double)a*sin((double)N*pi*(sqrt(fit[0]) - fit[1] + 1.)) - 1.;

	//Calculate constrain
	cons_value = t / (1. + exp(4.*abs(t)));

	//push value to the vector
	cons_val_temp.push_back(cons_value);

	//return vector
	return cons_val_temp;
}

/**********************CF3**************************/
std::vector<double> CF3::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2, MultipJ1, MultipJ2;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;
	MultipJ1 = MultipJ2 = 1;
	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - sin(6.*pi*code[0]+((double)i*pi/(double)Vars()));
		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += pow(y, 2);
			MultipJ1 *= cos((20. * y*pi) / sqrt(i));
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += pow(y, 2);
			MultipJ2 *= cos((20. * y*pi) / sqrt(i));
		}
	}

	std::vector<double> fitness;					//fitness vector - output

	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + (2.0 / sizeJ1)*(4.0 * SumJ1 - 2.0 * MultipJ1 + 2.));				//fitness 1
	fitness.push_back(1.0 - pow(code[0],2) + (2.0 / sizeJ2)*(4.0 * SumJ2 - 2.0 * MultipJ2 + 2.));		//fitness 2

																										//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void CF3::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -2;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 40,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 40,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> CF3::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	double N = 2.;		// CEC '09

						//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

	//calculate the real PF points
	//Push 1st point
	std::vector<double> temp_vect;		//temporary vector of fitness

	//push values to the fitness vector
	temp_vect.push_back(0.);
	temp_vect.push_back(1.);
	//push the fitness vector to the PF vector
	PF_real_val.push_back(temp_vect);


	//save the fitness values to the file
	if (indexr >= 0)
	{
		PF_real << 0. << " ";
		PF_real << 1.;
		PF_real << std::endl;
	}

	//calculate the step of the i value for the loop
	double i_step = ((sqrt(2./4.)-sqrt(1./4.))+ (sqrt(4. / 4.) - sqrt(3. / 4.))) / ((double)size - 2.0);

	//Calculate the rest of the points
	for (double i = sqrt(1. / 4.); i <= sqrt(2. / 4.); i += i_step)
	{
		//calculate the PF value
		double temp = (1. - pow(i, 2));				//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	for (double i = sqrt(3. / 4.); i <= sqrt(4. / 4.); i += i_step)
	{
		//calculate the PF value
		double temp = (1. - pow(i, 2));				//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> CF3::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	double cons_value;					//Value for constrain check
	short N = 2;		// CEC '09
	short a = 1;		// CEC '09
	std::vector<double> cons_val_temp;	//Values of the constrains - output
	
	//Calculate constrain
	cons_value = fit[1] + pow(fit[0], 2) - (double)a*sin((double)N*pi*(pow(fit[0], 2) - fit[1] + 1.)) - 1.;

	//push value to the vector
	cons_val_temp.push_back(cons_value);

	//return vector
	return cons_val_temp;
}

/**********************CF4**************************/
std::vector<double> CF4::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;
	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - sin(6.*pi*code[0] + ((double)i*pi / (double)Vars()));
		double h;

		if (i == 2)
		{
			if (y < (3. / 2.*(1. - (sqrt(2.) / 2.))))
				h = abs(y);
			else
				h = 0.125 + pow(y - 1., 2);
		}
		else
			h = pow(y, 2);

		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += h;
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += h;
		}
	}

	std::vector<double> fitness;					//fitness vector - output

													//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + SumJ1);				//fitness 1
	fitness.push_back(1.0 - code[0] + SumJ2);		//fitness 2

																										//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void CF4::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -2;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 25,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 25,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> CF4::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving

						//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

	//calculate the real PF points
	for (double i = 0; i <= 1.; i += 1.0 / (double)(size - 1))
	{
		double temp;				//PF value
		
		//calculate the PF value
		if (i <= 0.5)
			temp = 1. - i;
		else if (i <= 0.75)
			temp = -0.5*i + (3. / 4.);
		else
			temp = 1. - i + 0.125;
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}
	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> CF4::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double time)
{
	double cons_value;					//Value for constrain check
	double t;							//Temporary value
	std::vector<double> cons_val_temp;	//Values of the constrains - output

	//Calculate t
	t = code[1] - sin(6.*pi*code[0] + (2.*pi / (double)Vars()))-0.5*code[0]+0.25;

	//Calculate constrain
	cons_value = t / (1. + pow(e, 4.*abs(t)));

	//push value to the vector
	cons_val_temp.push_back(cons_value);

	//return vector
	return cons_val_temp;
}

/**********************CF5**************************/
std::vector<double> CF5::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;

	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y;
		if (i % 2 == 1)
			y = code[i - 1] - 0.8*code[0] * cos(6.*pi*code[0] + ((double)i*pi / (double)Vars()));
		else
			y = code[i - 1] - 0.8*code[0] * sin(6.*pi*code[0] + ((double)i*pi / (double)Vars()));
		double h;

		if (i == 2)
		{
			if (y < (3. / 2.*(1. - (sqrt(2.) / 2.))))
				h = abs(y);
			else
				h = 0.125 + pow(y - 1., 2);
		}
		else
			h = 2. * pow(y, 2) - cos(4.*pi*y) + 1.;

		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += h;
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += h;
		}
	}

	std::vector<double> fitness;					//fitness vector - output

	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + SumJ1);				//fitness 1
	fitness.push_back(1.0 - code[0] + SumJ2);		//fitness 2

													//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void CF5::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -2;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 35,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 35,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> CF5::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving

	//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

	//calculate the real PF points
	for (double i = 0; i <= 1.; i += 1.0 / (double)(size - 1))
	{
		double temp;				//PF value

		//calculate the PF value
		if (i <= 0.5)
			temp = 1. - i;
		else if (i <= 0.75)
			temp = -0.5*i + (3. / 4.);
		else
			temp = 1. - i + 0.125;
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}
	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> CF5::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	double cons_value;					//Value for constrain check
	std::vector<double> cons_val_temp;	//Values of the constrains - output


	//Calculate constrain value
	cons_value = code[1] - 0.8*code[0] * sin(6.*pi*code[0] + (2.*pi / (double)Vars())) - 0.5*code[0] + 0.25;

	//push value to the vector
	cons_val_temp.push_back(cons_value);

	//return vector
	return cons_val_temp;
}

/**********************CF6**************************/
std::vector<double> CF6::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;

	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y;
			
		if (i % 2 == 1)
		{
			y = code[i - 1] - 0.8*code[0] * cos(6.*pi*code[0] + ((double)i*pi / (double)Vars()));
			sizeJ1 += 1;
			SumJ1 += pow(y, 2);
		}
		else
		{
			y = code[i - 1] - 0.8*code[0] * sin(6.*pi*code[0] + ((double)i*pi / (double)Vars()));
			sizeJ2 += 1;
			SumJ2 += pow(y, 2);
		}
	}

	std::vector<double> fitness;					//fitness vector - output

	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + SumJ1);						//fitness 1
	fitness.push_back(pow(1.0 - code[0], 2) + SumJ2);		//fitness 2

													//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void CF6::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -2;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 19,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 1.5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> CF6::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving

											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

															//calculate the real PF points
	for (double i = 0.; i <= 1.; i += 1.0 / (double)(size - 1))
	{
		double temp;				//PF value

									//calculate the PF value
		if (i <= 0.5)
			temp = pow(1. - i, 2);
		else if (i <= 0.75)
			temp = 0.5*(1. - i);
		else
			temp = 0.25*sqrt(1. - i);
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}
	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> CF6::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	double cons_value;					//Value for the constrain check
	double cons_value2;					//2nd Value for the constrain check
	std::vector<double> cons_val_temp;	//Values of the constrains - output


	
	double temp, temp2;				//temporary values
	//Calculate temporary values
	temp = 0.5*(1. - code[0]) - pow(1. - code[0], 2);
	temp2 = 0.25 * sqrt(1. - code[0]) - 0.5*(1. - code[0]);

	//Calculate constrain values
	cons_value = code[1] - 0.8*code[0] * sin(6.*pi*code[0] + (2.*pi / (double)Vars())) - sgn(temp)*sqrt(abs(temp));
	cons_value2 = code[3] - 0.8*code[0] * sin(6.*pi*code[0] + (4.*pi / (double)Vars())) - sgn(temp2)*sqrt(abs(temp2));
	
	//push values to the vector
	cons_val_temp.push_back(cons_value);
	cons_val_temp.push_back(cons_value2);

	//return vector
	return cons_val_temp;
}

/**********************CF7**************************/
std::vector<double> CF7::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;

	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y;
		if (i % 2 == 1)
			y = code[i - 1] - cos(6.*pi*code[0] + ((double)i*pi / (double)Vars()));
		else
			y = code[i - 1] - sin(6.*pi*code[0] + ((double)i*pi / (double)Vars()));
		double h;

		if (i == 2 || i == 4)
			h = pow(y, 2);
		else
			h = 2. * pow(y, 2) - cos(4.*pi*y) + 1.;

		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += h;
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += h;
		}
	}

	std::vector<double> fitness;					//fitness vector - output

													//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + SumJ1);						//fitness 1
	fitness.push_back(pow(1.0 - code[0], 2) + SumJ2);		//fitness 2

															//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void CF7::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -2;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 40,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 40,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> CF7::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving

											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

															//calculate the real PF points
	for (double i = 0; i <= 1.; i += 1.0 / (double)(size - 1))
	{
		double temp;				//PF value

		//calculate the PF value
		if (i <= 0.5)
			temp = pow(1. - i, 2);
		else if (i <= 0.75)
			temp = 0.5*(1. - i);
		else
			temp = 0.25*sqrt(1. - i);
		std::vector<double> temp_vect;			//temporary vector of fitness

		//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}
	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> CF7::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	double cons_value;					//Value for the constrain check
	double cons_value2;					//2nd Value for the constrain check
	std::vector<double> cons_val_temp;	//Values of the constrains - output



	double temp, temp2;				//temporary values
									//Calculate temporary values
	temp = 0.5*(1 - code[0]) - pow(1 - code[0], 2);
	temp2 = 0.25 * sqrt(1 - code[0]) - 0.5*(1 - code[0]);

	//Calculate constrain values
	cons_value = code[1] - sin(6.*pi*code[0] + (2.*pi / (double)Vars())) - sgn(temp)*sqrt(abs(temp));
	cons_value2 = code[3] - sin(6.*pi*code[0] + (4.*pi / (double)Vars())) - sgn(temp2)*sqrt(abs(temp2));

	//push values to the vector
	cons_val_temp.push_back(cons_value);
	cons_val_temp.push_back(cons_value2);

	//return vector
	return cons_val_temp;
}

/**********************CF8**************************/
std::vector<double> CF8::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, sizeJ3, SumJ1, SumJ2, SumJ3;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = sizeJ3 = SumJ1 = SumJ2 = SumJ3 = 0;
	double f1, f2, f3;
	std::vector<double> fitness;					//fitness vector - output

	//calculate temp values
	for (int i = 3; i <= Vars(); i++)
	{
		double y = code[i - 1] - 2*code[1]*sin(2.*pi*code[0] + ((double)i*pi / (double)Vars()));
		if (i % 3 == 1)
		{
			sizeJ1++;
			SumJ1 += pow(y,2);
		}
		else if (i % 3 == 2)
		{
			sizeJ2++;
			SumJ2 += pow(y, 2);
		}
		else
		{
			sizeJ3++;
			SumJ3 += pow(y, 2);
		}
	}

	f1 = cos(0.5*code[0] * pi)*cos(0.5*code[1] * pi) + (2. / sizeJ1)*SumJ1;
	f2 = cos(0.5*code[0] * pi)*sin(0.5*code[1] * pi) + (2. / sizeJ2)*SumJ2;
	f3 = sin(0.5*code[0] * pi) + (2. / sizeJ3)*SumJ3;
													//calculate the fitness and push it to the fitness vector
	fitness.push_back(f1);						//fitness 1
	fitness.push_back(f2);		//fitness 2
	fitness.push_back(f3);		//fitness 3
															//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void CF8::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);
	STRUCTURES::boundaries second{ 1,0 };
	bound.push_back(second);

	//set boundaries for the rest of variables
	for (int i = 2; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 4;
		temp.lower = -4;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 40,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 40,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> CF8::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	int N = 2;
											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

															//calculate the real PF points
	for (double i = 0; i <= 1.; i += 1.0 / ((double)size / (2.*N + 1.) - 1.))
	{
		for (int j = 0; j <= 2 * N; j++)
		{
			double f1, f2, f3;

			f3 = i;
			f1 = sqrt((1 - pow(f3, 2))*j / (2 * N));

			if (pow(f1, 2) + pow(f3, 2) >= 1.)
				f2 = 0.;
			else
				f2 = sqrt(1 - pow(f1, 2) - pow(f3, 2));

			std::vector<double> temp_vect;			//temporary vector of fitness

													//push values to the fitness vector
			temp_vect.push_back(f1);
			temp_vect.push_back(f2);
			temp_vect.push_back(f3);

			//push the fitness vector to the PF vector
			PF_real_val.push_back(temp_vect);

			//save the fitness values to the file
			if (indexr >= 0)
			{
				PF_real << f1 << " ";
				PF_real << f2 << " ";
				PF_real << f3;
				PF_real << std::endl;
			}
		}
	}
	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> CF8::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	double cons_value;					//Value for the constrain check
	std::vector<double> cons_val_temp;	//Values of the constrains - output
	int a = 4;
	int N = 2;

	double temp1, temp2;				//temporary values
	//Calculate temporary values
	temp1 = (pow(fit[0], 2) + pow(fit[1], 2)) / (1. - pow(fit[2], 2));
	temp2 = a*abs(sin(N*pi* ((pow(fit[0], 2) - pow(fit[1], 2)) / (1. - pow(fit[2], 2)) + 1.)));

	//Calculate constrain values
	cons_value = temp1 - temp2 - 1;

	//push value to the vector
	cons_val_temp.push_back(cons_value);

	//return vector
	return cons_val_temp;
}

/**********************CF9**************************/
std::vector<double> CF9::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, sizeJ3, SumJ1, SumJ2, SumJ3;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = sizeJ3 = SumJ1 = SumJ2 = SumJ3 = 0;
	double f1, f2, f3;
	std::vector<double> fitness;					//fitness vector - output

													//calculate temp values
	for (int i = 3; i <= Vars(); i++)
	{
		double y = code[i - 1] - 2 * code[1] * sin(2.*pi*code[0] + ((double)i*pi / (double)Vars()));
		if (i % 3 == 1)
		{
			sizeJ1++;
			SumJ1 += pow(y, 2);
		}
		else if (i % 3 == 2)
		{
			sizeJ2++;
			SumJ2 += pow(y, 2);
		}
		else
		{
			sizeJ3++;
			SumJ3 += pow(y, 2);
		}
	}

	f1 = cos(0.5*code[0] * pi)*cos(0.5*code[1] * pi) + (2. / sizeJ1)*SumJ1;
	f2 = cos(0.5*code[0] * pi)*sin(0.5*code[1] * pi) + (2. / sizeJ2)*SumJ2;
	f3 = sin(0.5*code[0] * pi) + (2. / sizeJ3)*SumJ3;
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(f1);						//fitness 1
	fitness.push_back(f2);		//fitness 2
	fitness.push_back(f3);		//fitness 3
								//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void CF9::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);
	STRUCTURES::boundaries second{ 1,0 };
	bound.push_back(second);

	//set boundaries for the rest of variables
	for (int i = 2; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -2;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 40,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 40,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> CF9::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	int N = 2;
	//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

															//calculate the real PF points
															//line part
	for (double i = 0; i <= 1.; i += 1.0 / ((double)size / 40 - 1.))
	{

		double f1, f2, f3;

		f1 = 0;
		f2 = i;
		f3 = sqrt(1. - pow(f2, 2));

		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(f1);
		temp_vect.push_back(f2);
		temp_vect.push_back(f3);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		f1 = 0;
		f3 = i;
		f2 = sqrt(1. - pow(f3, 2));

		//push values to the fitness vector
		temp_vect[0] = f1;
		temp_vect[1] = f2;
		temp_vect[2] = f3;

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		f2 = 0;
		f1 = i;
		f3 = sqrt(1. - pow(f1, 2));


		//push values to the fitness vector
		temp_vect[0] = f1;
		temp_vect[1] = f2;
		temp_vect[2] = f3;

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		f2 = 0;
		f3 = i;
		f1 = sqrt(1. - pow(f3, 2));


		//push values to the fitness vector
		temp_vect[0] = f1;
		temp_vect[1] = f2;
		temp_vect[2] = f3;

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);
	}
	int size_f3, size_f1;
	int size_multip = 20 * N;
	size_f1 = sqrt((float)size / size_multip);
	size_f3 = size_f1 * size_multip;
	//surfaces part
	for (double i_f3 = 0; i_f3 <= 1.; i_f3 += 1.0 / (size_f3 - 1.))
	{
		double f3 = i_f3;
		for (int i_n = 1; i_n <= N; i_n++)
		{
			double min_f1, max_f1, step_f1;
			min_f1 = sqrt((2.*i_n - 1) / (2.*N)*(1 - pow(f3, 2)));
			max_f1 = sqrt((2.*i_n) / (2.*N)*(1 - pow(f3, 2)));
			if (1 - pow(f3, 2) < 0)
				abort();
			step_f1 = (max_f1 - min_f1) / (size_f1 - 1.);
			for (double i_f1 = min_f1; i_f1 <= max_f1; i_f1 += step_f1)
			{
				double f1, f2;
				f1 = i_f1;
				if (pow(f1, 2) + pow(f3, 2) >= 1)
					f2 = 0;
				else
					f2 = sqrt(1. - pow(f1, 2) - pow(f3, 2));

				if (f2 == 0)
					continue;

				std::vector<double> temp_vect;			//temporary vector of fitness


														//push values to the fitness vector
				temp_vect.push_back(f1);
				temp_vect.push_back(f2);
				temp_vect.push_back(f3);

				//push the fitness vector to the PF vector
				PF_real_val.push_back(temp_vect);

			}
		}
	}

	//refine PF
	short nobj = Objs();
	for (int i = 0; i < PF_real_val.size(); i++)
	{
		for (int j = i + 1; j < PF_real_val.size(); j++)
		{
			int eq = 0;
			for (int k = 0; k < nobj; k++)
			{
				if (abs(PF_real_val[i][k] - PF_real_val[j][k]) <= pow(10, -pow_eq_zero))
					eq++;
				else
					break;
			}
			if (eq == nobj)
			{
				PF_real_val.erase(PF_real_val.begin() + j);
				j--;
				continue;
			}
		}
		if (indexr >= 0)
		{
			for (int k = 0; k < nobj; k++)
			{
				PF_real << PF_real_val[i][k] << " ";
			}
			PF_real << std::endl;
		}
	}
	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> CF9::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	double cons_value;					//Value for the constrain check
	std::vector<double> cons_val_temp;	//Values of the constrains - output
	int a = 3;
	int N = 2;

	double temp1, temp2;				//temporary values
										//Calculate temporary values
	temp1 = (pow(fit[0], 2) + pow(fit[1], 2)) / (1. - pow(fit[2], 2));
	temp2 = a*sin(N*pi* ((pow(fit[0], 2) - pow(fit[1], 2)) / (1. - pow(fit[2], 2)) + 1.));

	//Calculate constrain values
	cons_value = temp1 - temp2 - 1;

	//push value to the vector
	cons_val_temp.push_back(cons_value);

	//return vector
	return cons_val_temp;
}

/**********************CF10**************************/
std::vector<double> CF10::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, sizeJ3, SumJ1, SumJ2, SumJ3;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = sizeJ3 = SumJ1 = SumJ2 = SumJ3 = 0;
	double f1, f2, f3;
	std::vector<double> fitness;					//fitness vector - output

													//calculate temp values
	for (int i = 3; i <= Vars(); i++)
	{
		double y = code[i - 1] - 2 * code[1] * sin(2.*pi*code[0] + ((double)i*pi / (double)Vars()));
		if (i % 3 == 1)
		{
			sizeJ1++;
			SumJ1 += 4 * pow(y, 2) - cos(8 * pi*y) + 1.;
		}
		else if (i % 3 == 2)
		{
			sizeJ2++;
			SumJ2 += 4 * pow(y, 2) - cos(8 * pi*y) + 1.;
		}
		else
		{
			sizeJ3++;
			SumJ3 += 4 * pow(y, 2) - cos(8 * pi*y) + 1.;
		}
	}

	f1 = cos(0.5*code[0] * pi)*cos(0.5*code[1] * pi) + (2. / sizeJ1)*SumJ1;
	f2 = cos(0.5*code[0] * pi)*sin(0.5*code[1] * pi) + (2. / sizeJ2)*SumJ2;
	f3 = sin(0.5*code[0] * pi) + (2. / sizeJ3)*SumJ3;
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(f1);						//fitness 1
	fitness.push_back(f2);		//fitness 2
	fitness.push_back(f3);		//fitness 3
								//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void CF10::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);
	STRUCTURES::boundaries second{ 1,0 };
	bound.push_back(second);

	//set boundaries for the rest of variables
	for (int i = 2; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -2;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 40,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 40,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> CF10::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	int N = 2;
	//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

															//calculate the real PF points
															//line part
	for (double i = 0; i <= 1.; i += 1.0 / ((double)size / 40 - 1.))
	{

		double f1, f2, f3;

		f1 = 0;
		f2 = i;
		f3 = sqrt(1. - pow(f2, 2));

		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(f1);
		temp_vect.push_back(f2);
		temp_vect.push_back(f3);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		f1 = 0;
		f3 = i;
		f2 = sqrt(1. - pow(f3, 2));

		//push values to the fitness vector
		temp_vect[0] = f1;
		temp_vect[1] = f2;
		temp_vect[2] = f3;

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		f2 = 0;
		f1 = i;
		f3 = sqrt(1. - pow(f1, 2));


		//push values to the fitness vector
		temp_vect[0] = f1;
		temp_vect[1] = f2;
		temp_vect[2] = f3;

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		f2 = 0;
		f3 = i;
		f1 = sqrt(1. - pow(f3, 2));


		//push values to the fitness vector
		temp_vect[0] = f1;
		temp_vect[1] = f2;
		temp_vect[2] = f3;

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);
	}
	int size_f3, size_f1;
	int size_multip = 20 * N;
	size_f1 = sqrt((float)size / size_multip);
	size_f3 = size_f1 * size_multip;
	//surfaces part
	for (double i_f3 = 0; i_f3 <= 1.; i_f3 += 1.0 / (size_f3 - 1.))
	{
		double f3 = i_f3;
		for (int i_n = 1; i_n <= N; i_n++)
		{
			double min_f1, max_f1, step_f1;
			min_f1 = sqrt((2.*i_n - 1) / (2.*N)*(1 - pow(f3, 2)));
			max_f1 = sqrt((2.*i_n) / (2.*N)*(1 - pow(f3, 2)));
			if (1 - pow(f3, 2) < 0)
				abort();
			step_f1 = (max_f1 - min_f1) / (size_f1 - 1.);
			for (double i_f1 = min_f1; i_f1 <= max_f1; i_f1 += step_f1)
			{
				double f1, f2;
				f1 = i_f1;
				if (pow(f1, 2) + pow(f3, 2) >= 1)
					f2 = 0;
				else
					f2 = sqrt(1. - pow(f1, 2) - pow(f3, 2));

				if (f2 == 0)
					continue;

				std::vector<double> temp_vect;			//temporary vector of fitness


														//push values to the fitness vector
				temp_vect.push_back(f1);
				temp_vect.push_back(f2);
				temp_vect.push_back(f3);

				//push the fitness vector to the PF vector
				PF_real_val.push_back(temp_vect);

			}
		}
	}

	//refine PF
	short nobj = Objs();
	for (int i = 0; i < PF_real_val.size(); i++)
	{
		for (int j = i + 1; j < PF_real_val.size(); j++)
		{
			int eq = 0;
			for (int k = 0; k < nobj; k++)
			{
				if (abs(PF_real_val[i][k] - PF_real_val[j][k]) <= pow(10, -pow_eq_zero))
					eq++;
				else
					break;
			}
			if (eq == nobj)
			{
				PF_real_val.erase(PF_real_val.begin() + j);
				j--;
				continue;
			}
		}
		if (indexr >= 0)
		{
			for (int k = 0; k < nobj; k++)
			{
				PF_real << PF_real_val[i][k] << " ";
			}
			PF_real << std::endl;
		}
	}
	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> CF10::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	double cons_value;					//Value for the constrain check
	std::vector<double> cons_val_temp;	//Values of the constrains - output
	int a = 1;
	int N = 2;

	double temp1, temp2;				//temporary values
										//Calculate temporary values
	temp1 = (pow(fit[0], 2) + pow(fit[1], 2)) / (1. - pow(fit[2], 2));
	temp2 = a*sin(N*pi* ((pow(fit[0], 2) - pow(fit[1], 2)) / (1. - pow(fit[2], 2)) + 1.));

	//Calculate constrain values
	cons_value = temp1 - temp2 - 1;

	//push value to the vector
	cons_val_temp.push_back(cons_value);

	//return vector
	return cons_val_temp;
}
/**********************DTLZ8**************************/
std::vector<double> DTLZ8::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;					//fitness vector - output

	short nobj = Objs();
	short nvar = Vars();
	//calculate temp values
	for (int i = 0; i < nobj; i++)
	{
		double temp_sum = 0;

		for (int j = floor(i*(double)nvar / nobj); j < floor((i + 1)*(double)nvar / nobj); j++)
		{
			temp_sum += code[j];
		}

		double temp_fit = (1. / floor((double)nvar / nobj))*temp_sum;
		fitness.push_back(temp_fit);
	}


															//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void DTLZ8::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

	//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 1.5,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 1.5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> DTLZ8::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving

											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output
	short n_obj = Objs();

	double step = 1. / (floor(sqrt(size * 11)) - 1);
	double f_m;
	//line part
	for (f_m = 1; f_m >= 0; f_m -= step)
	{
		//calculate the PF value
		double f = (1 - f_m) / 4.;			//PF value

		if (2 * f_m + 2 * f - 1 < 0)
			if (n_obj != 2)
				break;


		std::vector<double> temp_vect;		//temporary vector of fitness

		for (int i = 0; i < n_obj; i++)
		{
			double temp;
			if (i != n_obj - 1)
				temp = f;
			else
				temp = f_m;

			//push values to the fitness vector
			temp_vect.push_back(temp);
		}

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);


	}
	//hyperplane part
	for (; f_m > 0 - step; f_m -= step)
	{
		if (f_m < 0.)
			f_m = 0;
		//calculate the min function value
		double min_f = 1 - 2 * f_m;			//min f value


		for (double i_min_f = 0; i_min_f < min_f / 2; i_min_f += step)
		{
			double lower_f = i_min_f;
			double upper_f = min_f - i_min_f;

			if (f_m + 4 * lower_f - 1 < 0)
				continue;


			for (int i_obj = 0; i_obj < n_obj - 1; i_obj++)
			{
				std::vector<double> temp_vect;
				for (int i = 0; i < n_obj; i++)
				{
					double temp;
					if (i != n_obj - 1)
					{
						if (i == i_obj)
							temp = lower_f;
						else
							temp = upper_f;
					}
					else
						temp = f_m;

					//push values to the fitness vector
					temp_vect.push_back(temp);
				}
				//push the fitness vector to the PF vector
				PF_real_val.push_back(temp_vect);
			}
		}
		double lower_f = min_f / 2.;

		if (f_m + 4 * lower_f - 1 < 0)
			continue;

		std::vector<double> temp_vect;
		for (int i = 0; i < n_obj; i++)
		{
			double temp;
			if (i != n_obj - 1)
			{
				temp = lower_f;
			}
			else
				temp = f_m;

			//push values to the fitness vector
			temp_vect.push_back(temp);
		}
		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);




	}


	//save the PF to the file
	if (indexr >= 0)
	{
		for (int i = 0; i < PF_real_val.size(); i++)
		{
			for (int j = 0; j < n_obj; j++)
				PF_real << PF_real_val[i][j] << " ";
			PF_real << std::endl;
		}
	}

	return PF_real_val;
}

std::vector<double> DTLZ8::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{

	std::vector<double> cons_val_temp;	//Values of the constrains - output
	short nobj = Objs();

	for (int i = 0; i < nobj - 1; i++)
	{
		double cons_value = fit[nobj - 1] + 4*fit[i] - 1.;
		cons_val_temp.push_back(cons_value);
	}
	
	double min_fit_val = 1e30;
	for (int i = 0; i < nobj - 1; i++)
	{
		for (int j = 0; j < nobj - 1; j++)
		{
			if (i == j)
				continue;
			else
			{
				double temp_fit = fit[i] + fit[j];
				if (temp_fit < min_fit_val)
					min_fit_val = temp_fit;
			}
		}
	}


	double cons_value = 2.*fit[nobj - 1] + min_fit_val - 1.;
	cons_val_temp.push_back(cons_value);

	if (cons_val_temp.size() != Cons())
		abort();

	//return vector
	return cons_val_temp;
}

/**********************DTLZ9**************************/
std::vector<double> DTLZ9::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;					//fitness vector - output

	short nobj = Objs();
	short nvar = Vars();
	//calculate temp values
	for (int i = 0; i < nobj; i++)
	{
		double temp_sum = 0;

		for (int j = floor(i*(double)nvar / nobj); j < floor((i + 1)*(double)nvar / nobj); j++)
		{
			temp_sum += pow(code[j], 0.1);
		}

		double temp_fit = temp_sum;
		fitness.push_back(temp_fit);
	}

															//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void DTLZ9::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> DTLZ9::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;				//ofstream file for the PF saving

										//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output
	short n_obj = Objs();

	//calculate the real PF points
	for (double f_m = 0; f_m <= 1; f_m += 1.0 / (double)(size / 2 - 1))
	{
		//calculate the PF value
		double f = sqrt(1 - pow(f_m, 2));			//PF value

		std::vector<double> temp_vect;		//temporary vector of fitness

		for (int i = 0; i < n_obj; i++)
		{
			double temp;
			if (i != n_obj - 1)
				temp = f;
			else
				temp = f_m;

			//push values to the fitness vector
			temp_vect.push_back(temp);

			if (indexr >= 0)
				PF_real << temp << " ";
		}
		if (indexr >= 0)
			PF_real << std::endl;

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

	}

	for (double f = 0; f <= 1; f += 1.0 / (double)(size / 2 - 1))
	{
		//calculate the PF value
		double f_m = sqrt(1 - pow(f, 2));			//PF value

		std::vector<double> temp_vect;		//temporary vector of fitness

		for (int i = 0; i < n_obj; i++)
		{
			double temp;
			if (i != n_obj - 1)
				temp = f;
			else
				temp = f_m;

			//push values to the fitness vector
			temp_vect.push_back(temp);

			if (indexr >= 0)
				PF_real << temp << " ";
		}
		if (indexr >= 0)
			PF_real << std::endl;

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

	}


	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> DTLZ9::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	std::vector<double> cons_val_temp;	//Values of the constrains - output
	short nobj = Objs();

	for (int i = 0; i < nobj - 1; i++)
	{
		double cons_value = pow(fit[nobj - 1], 2) + pow(fit[i], 2) - 1.;
		cons_val_temp.push_back(cons_value);
	}

	if (cons_val_temp.size() != Cons())
		abort();

	//return vector
	return cons_val_temp;
}

/**********************IMB11**************************/
std::vector<double> IMB11::Fitness_C(const std::vector<double> & code)
{
	double g_x = 0;			//temp values for the fitness calculation
	double f1, f2, f3;
	std::vector<double> fitness;					//fitness vector - output

													//calculate temp values


	for (int i = 2; i <= Vars(); i++)
	{
		double t = code[i - 1] - code[0];
		g_x += 0.5*(-0.9*pow(t, 2) + pow(abs(t), 0.6));
	}


	f1 = code[0] * (1 + g_x);
	f2 = (1 - sqrt(code[0])) * (1 + g_x);
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(f1);						//fitness 1
	fitness.push_back(f2);		//fitness 2
								//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void IMB11::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes



																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 2,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 3,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> IMB11::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = 1. - sqrt(i);					//PF value

		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> IMB11::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	std::vector<double> cons_val_temp;	//Values of the constrains - output

	double g_x = 0;
	double G_x = 0;

	for (int i = 2; i <= Vars(); i++)
	{
		double t = code[i - 1] - code[0];
		g_x += 0.5*(-0.9*pow(t, 2) + pow(abs(t), 0.6));
	}
	if (code[0] > 0.6 && g_x > 0.001)
		G_x = 10 * g_x;

	cons_val_temp.push_back(-abs(G_x));

	if (cons_val_temp.size() != Cons())
		abort();

	//return vector
	return cons_val_temp;
}

/**********************IMB12**************************/
std::vector<double> IMB12::Fitness_C(const std::vector<double> & code)
{
	double g_x = 0;			//temp values for the fitness calculation
	double f1, f2, f3;
	std::vector<double> fitness;					//fitness vector - output

													//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double t = code[i - 1] - code[0];
		g_x += 0.5*(-0.9*pow(t, 2) + pow(abs(t), 0.6));
	}


	f1 = code[0] * (1 + g_x);
	f2 = (1 - code[0]) * (1 + g_x);
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(f1);						//fitness 1
	fitness.push_back(f2);		//fitness 2
								//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void IMB12::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes



																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 3,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 2,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> IMB12::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = 1. - i;					//PF value

		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> IMB12::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	std::vector<double> cons_val_temp;	//Values of the constrains - output

	double g_x = 0;
	double G_x = 0;

	for (int i = 2; i <= Vars(); i++)
	{
		double t = code[i - 1] - code[0];
		g_x += 0.5*(-0.9*pow(t, 2) + pow(abs(t), 0.6));
	}
	if ((code[0] < 0.2 || code[0] > 0.8) && g_x > 0.002)
		G_x = 10 * g_x;

	cons_val_temp.push_back(-abs(G_x));

	if (cons_val_temp.size() != Cons())
		abort();

	//return vector
	return cons_val_temp;
}
/**********************IMB13**************************/
std::vector<double> IMB13::Fitness_C(const std::vector<double> & code)
{
	double g_x = 0;			//temp values for the fitness calculation
	double f1, f2, f3;
	std::vector<double> fitness;					//fitness vector - output

													//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double t = code[i - 1] - code[0];
		g_x += 0.5*(-0.9*pow(t, 2) + pow(abs(t), 0.6));
	}


	f1 = code[0] * (1 + g_x);
	f2 = (1 - pow(code[0], 2)) * (1 + g_x);
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(f1);						//fitness 1
	fitness.push_back(f2);		//fitness 2
								//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void IMB13::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes



																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 3,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 3,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> IMB13::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//calculate the real PF points
	for (double i = 0; i <= 1; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = 1. - pow(i, 2);					//PF value

		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}
std::vector<double> IMB13::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	std::vector<double> cons_val_temp;	//Values of the constrains - output

	double g_x = 0;
	double G_x = 0;

	for (int i = 2; i <= Vars(); i++)
	{
		double t = code[i - 1] - code[0];
		g_x += 0.5*(-0.9*pow(t, 2) + pow(abs(t), 0.6));
	}
	if ((code[0] < 0.2 || code[0] > 0.8) && g_x > 0.001)
		G_x = 10 * g_x;

	cons_val_temp.push_back(-abs(G_x));

	if (cons_val_temp.size() != Cons())
		abort();

	//return vector
	return cons_val_temp;
}

/**********************IMB14**************************/
std::vector<double> IMB14::Fitness_C(const std::vector<double> & code)
{
	double g_x = 0;			//temp values for the fitness calculation
	double f1, f2, f3;
	std::vector<double> fitness;					//fitness vector - output

													//calculate temp values
	for (int i = 3; i <= Vars(); i++)
	{
		double t = code[i - 1] - 0.5*(code[0] + code[1]);
		g_x += 0.5*(-0.9*pow(t, 2) + pow(abs(t), 0.6));
	}


	f1 = code[0] * code[1] * (1 + g_x);
	f2 = code[0] * (1 - code[1]) * (1 + g_x);
	f3 = (1 - code[0]) * (1 + g_x);
	//calculate the fitness and push it to the fitness vector
	fitness.push_back(f1);						//fitness 1
	fitness.push_back(f2);						//fitness 2
	fitness.push_back(f3);						//fitness 3
												//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void IMB14::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes



																		//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 40,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 40,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> IMB14::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;		//ofstream file for the PF saving

								//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;		//Vector for PF data storage - output

														//Get the uniform vector
	PF_real_val = Uniform_Weights_Generate(size, Objs(), 3);

	//calculate the real PF points
	for (double i = 0; i < PF_real_val.size(); i++)
	{
		for (int j = 0; j < Objs(); j++)
		{
			if (indexr >= 0)
				PF_real << PF_real_val[i][j] << " ";
		}
		if (indexr >= 0)
			PF_real << std::endl;
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}
std::vector<double> IMB14::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	std::vector<double> cons_val_temp;	//Values of the constrains - output

	double g_x = 0;
	double G_x = 0;

	for (int i = 3; i <= Vars(); i++)
	{
		double t = code[i - 1] - 0.5*(code[0] + code[1]);
		g_x += 0.5*(-0.9*pow(t, 2) + pow(abs(t), 0.6));
	}
	if (code[0] > 0.75 && g_x > 0.001)
		G_x = 10 * g_x;

	cons_val_temp.push_back(-abs(G_x));

	if (cons_val_temp.size() != Cons())
		abort();
	//return vector
	return cons_val_temp;
}

/********************************************
DYNAMIC UNCONSTRAINED FUNCTIONS
*********************************************/

/**********************UDF1**************************/
std::vector<double> UDF1::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;						//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;

	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time

	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - sin((6. * pi * code[0]) + (((double)i * pi) / (double)Vars())) - Gt;
		if (i % 2 == 1)
		{
			sizeJ1++;
			SumJ1 += pow(y, 2);
		}
		else
		{
			sizeJ2++;
			SumJ2 += pow(y, 2);
		}
	}

	std::vector<double> fitness;									//fitness vector - output

																	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + abs(Gt) + (2.0 / sizeJ1)*SumJ1);			//fitness 1
	fitness.push_back(1. - code[0] + abs(Gt) + (2.0 / sizeJ2)*SumJ2);	//fitness 2

																	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;

}

void UDF1::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -2;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> UDF1::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

	double Gt = sin(0.5*pi*t);

	//calculate the real PF points
	for (double i = 0. + abs(Gt); i <= 1. + abs(Gt); i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1. - i + 2*abs(Gt));		//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************UDF2**************************/
std::vector<double> UDF2::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;						//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;

	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time

									//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - Gt - pow(code[0], (0.5*(2.0 + 3.0*(i - 2) / (double)(Vars() - 2)) + Gt));
		if (i % 2 == 1)
		{
			sizeJ1++;
			SumJ1 += pow(y, 2);
		}
		else
		{
			sizeJ2++;
			SumJ2 += pow(y, 2);
		}
	}

	std::vector<double> fitness;									//fitness vector - output

																	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + abs(Gt) + (2.0 / sizeJ1)*SumJ1);			//fitness 1
	fitness.push_back(1. - code[0] + abs(Gt) + (2.0 / sizeJ2)*SumJ2);	//fitness 2

																	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;

}

void UDF2::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> UDF2::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving

											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

	double Gt = sin(0.5*pi*t);

	//calculate the real PF points
	for (double i = 0. + abs(Gt); i <= 1. + abs(Gt); i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1. - i + 2 * abs(Gt));		//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************UDF3**************************/
std::vector<double> UDF3::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	int N = 2;
	float eps = 0.1f;
	double sizeJ1, sizeJ2, SumJ1, SumJ2, MultipJ1, MultipJ2;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;
	MultipJ1 = MultipJ2 = 1;
	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - sin(6*pi*code[0]+i*pi/(double)Vars());
		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += pow(y, 2);
			MultipJ1 *= cos((20 * y*pi) / sqrt(i));
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += pow(y, 2);
			MultipJ2 *= cos((20 * y*pi) / sqrt(i));
		}
	}

	std::vector<double> fitness;					//fitness vector - output
	//calculate the fitness and push it to the fitness vector
	double maxf = (1. / N + 2 * eps)*(sin(2 * N*pi*code[0]) - abs(2 * N*Gt));		//Max function
	//Check maxf
	if (maxf < 0.)
		maxf = 0.;
	double fitness1_temp = code[0] + (2.0 / sizeJ1)*pow((4.0 * SumJ1 - 2.0 * MultipJ1 + 2.), 2) + maxf;				//fitness 1
	double fitness2_temp = 1.0 - code[0] + (2.0 / sizeJ2)*pow((4.0 * SumJ2 - 2.0 * MultipJ2 + 2.), 2) + maxf;		//fitness 2

	fitness.push_back(fitness1_temp);				//fitness 1
	fitness.push_back(fitness2_temp);			//fitness 2

																										//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void UDF3::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 50,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 50,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> UDF3::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	short N = 2;
	float eps = 0.1f;

									//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

															//calculate the real PF points
															//Push 1st point
	std::vector<double> temp_vect;		//temporary vector of fitness

	for (double i = 0. ; i <= 1. ; i += 1.0 / (double)(size - 1))
	{
		if (sin(2 * N*pi*i) > abs(2 * N*Gt))
			continue;


		//calculate the PF value
		double temp = i;
		double temp2 = 1. - i;				//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(temp);
		temp_vect.push_back(temp2);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << temp << " ";
			PF_real << temp2;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}
/**********************UDF4**************************/
std::vector<double> UDF4::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;						//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;

	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double Mt = 0.5 + abs(Gt);		//Value of M(t) function - dynamic dependant on time
	double Kt = ceil(Vars()*Gt);					//Value of K(t) function - dynamic dependant on time	

									//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - sin((6. * pi * code[0]) + (((double)i + Kt) * pi) / (double)Vars());
		if (i % 2 == 1)
		{
			sizeJ1++;
			SumJ1 += pow(y, 2);
		}
		else
		{
			sizeJ2++;
			SumJ2 += pow(y, 2);
		}
	}

	std::vector<double> fitness;									//fitness vector - output

	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + (2.0 / sizeJ1)*SumJ1);						//fitness 1
	fitness.push_back(1. - (Mt*pow(code[0], Mt)) + (2.0 / sizeJ2)*SumJ2);	//fitness 2

																	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void UDF4::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> UDF4::Plot_PF(int indexr, double t, int size)
{
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double Mt = 0.5 + abs(Gt);					//Value of M(t) function - dynamic dependant on time
	double Ht = Mt;					//Value of H(t) function - dynamic dependant on time
	
	std::ofstream PF_real;					//ofstream file for the PF saving

											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output


	//calculate the real PF points
	for (double i = 0.; i <= 1.; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1. - (Mt*pow(i, Ht)));		//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************UDF5**************************/
std::vector<double> UDF5::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;						//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;

	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double Mt = 0.5 + abs(Gt);		//Value of M(t) function - dynamic dependant on time

	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - Gt - pow(code[0], (0.5*(2.0 + 3.0*(i - 2) / (double)(Vars() - 2)) + Gt));
		if (i % 2 == 1)
		{
			sizeJ1++;
			SumJ1 += pow(y, 2);
		}
		else
		{
			sizeJ2++;
			SumJ2 += pow(y, 2);
		}
	}

	std::vector<double> fitness;									//fitness vector - output

																	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + (2.0 / sizeJ1)*SumJ1);			//fitness 1
	fitness.push_back(1. - Mt*pow(code[0], Mt) + (2.0 / sizeJ2)*SumJ2);	//fitness 2

																	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;

}

void UDF5::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> UDF5::Plot_PF(int indexr, double t, int size)
{
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double Mt = 0.5 + abs(Gt);					//Value of M(t) function - dynamic dependant on time
	double Ht = Mt;					//Value of H(t) function - dynamic dependant on time

	std::ofstream PF_real;					//ofstream file for the PF saving

											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output


															//calculate the real PF points
	for (double i = 0.; i <= 1.; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1. - (Mt*pow(i, Ht)));		//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************UDF6**************************/
std::vector<double> UDF6::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double Mt = 0.5 + abs(Gt);					//Value of M(t) function - dynamic dependant on time
	int N = 10;
	float eps = 0.1f;
	double sizeJ1, sizeJ2, SumJ1, SumJ2;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;
	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - sin(6 * pi*code[0] + i*pi / Vars());
		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += (2 * pow(y, 2)) - cos(4 * pi*y) + 1.;
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += (2 * pow(y, 2)) - cos(4 * pi*y) + 1.;
		}
	}

	std::vector<double> fitness;					//fitness vector - output
													//calculate the fitness and push it to the fitness vector
	double h = (1. / (2. * N) + eps)*abs(sin(2 * N*pi*code[0]) - abs(2 * N*Gt));

	double fitness1_temp = code[0] + h + (2. / sizeJ1)*SumJ1;			//fitness 1
	double fitness2_temp = 1 - Mt * code[0] + h + (2. / sizeJ2)*SumJ2;		//fitness 2
	fitness.push_back(fitness1_temp);				//fitness 1
	fitness.push_back(fitness2_temp);			//fitness 2

												//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void UDF6::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> UDF6::Plot_PF(int indexr, double t, int size)
{


	std::ofstream PF_real;					//ofstream file for the PF saving
	short N = 10;
	float eps = 0.1f;
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double Mt = 0.5 + abs(Gt);		//Value of M(t) function - dynamic dependant on time
									//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

//calculate the real PF points


	if (abs(Gt) <= pow(10, -PF_res))
	{
		//save the fitness values to the file
		for (int i = 0; i <= 2 * N; i++)
		{
			//calculate the PF value
			double temp = (i / (2.*N)) ;				//PF x value
			double temp2 = 1. - (i / (2. * N))*0.5;			//PF y value
			std::vector<double> temp_vect;			//temporary vector of fitness

													//push values to the fitness vector

			temp_vect.push_back(temp);
			temp_vect.push_back(temp2);

			//push the fitness vector to the PF vector
			PF_real_val.push_back(temp_vect);

			//save the fitness values to the file
			if (indexr >= 0)
			{
				PF_real << temp << " ";
				PF_real << temp2;
				PF_real << std::endl;
			}
		}
	}
	else
	{
		//save the fitness values to the file
		for (int i = 0; i < 2 * N; i++)
		{
			//skip values where sin(x) == -1
			if (i % 2 == 1)
				continue;

			//calculate the PF value
			double temp = (i / (2.*N)) + ((1. / (2.*N) + eps)*(abs(2.*N*Gt) - 1.)) + (1. / (4.*N));		//PF x value
			double temp2 = 1. - (i / (2. * N) + (1. / (4.*N)))*Mt + ((1. / (2.*N) + eps)*(abs(2.*N*Gt) - 1.));			//PF y value
			std::vector<double> temp_vect;			//temporary vector of fitness

													//push values to the fitness vector

			temp_vect.push_back(temp);
			temp_vect.push_back(temp2);

			//push the fitness vector to the PF vector
			PF_real_val.push_back(temp_vect);

			//save the fitness values to the file
			if (indexr >= 0)
			{
				PF_real << temp << " ";
				PF_real << temp2;
				PF_real << std::endl;
			}
		}
	}


	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************UDF8**************************/
std::vector<double> UDF8::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double Gt1 = sin(0.5*pi*t_vector[0]);		//Value of G(t1) function - dynamic dependant on time
	double Gt2 = sin(0.5*pi*t_vector[1]);		//Value of G(t2) function - dynamic dependant on time
	double Gt3 = sin(0.5*pi*t_vector[2]);		//Value of G(t3) function - dynamic dependant on time
	double Gt4 = sin(0.5*pi*t_vector[3]);		//Value of G(t4) function - dynamic dependant on time
	double Gt5 = sin(0.5*pi*t_vector[4]);		//Value of G(t5) function - dynamic dependant on time
	double Kt1 = ceil(Vars()*Gt1);				//Value of K(t1) function - dynamic dependant on time
	double Ht4 = 0.5 + abs(Gt4);				//Value of H(t4) function - dynamic dependant on time
	double Ht5 = 0.5 + abs(Gt5);				//Value of H(t5) function - dynamic dependant on time

	double sizeJ1, sizeJ2, SumJ1, SumJ2;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;
	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - sin(6 * pi*code[0] + ((double)i + Kt1)*pi / Vars()) - Gt2;
		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += pow(y, 2);
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += pow(y, 2);
		}
	}

	std::vector<double> fitness;					//fitness vector - output

	double fitness1_temp = code[0] + abs(Gt3) + (2. / sizeJ1)*SumJ1;			//fitness 1
	double fitness2_temp = 1 - Ht4*pow(code[0], Ht5) + abs(Gt3) + (2. / sizeJ2)*SumJ2;		//fitness 2
	fitness.push_back(fitness1_temp);				//fitness 1
	fitness.push_back(fitness2_temp);			//fitness 2

												//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void UDF8::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> UDF8::Plot_PF(int indexr, double t, int size)
{

	std::ofstream PF_real;					//ofstream file for the PF saving

	double Gt3 = sin(0.5*pi*t_vector[2]);		//Value of G(t3) function - dynamic dependant on time
	double Gt4 = sin(0.5*pi*t_vector[3]);		//Value of G(t4) function - dynamic dependant on time
	double Gt5 = sin(0.5*pi*t_vector[4]);		//Value of G(t5) function - dynamic dependant on time
	double Ht4 = 0.5 + abs(Gt4);				//Value of H(t4) function - dynamic dependant on time
	double Ht5 = 0.5 + abs(Gt5);				//Value of H(t5) function - dynamic dependant on time

	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

															//calculate the real PF points


	//calculate the real PF points
	for (double i = 0.+ abs(Gt3); i <= 1. + abs(Gt3); i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = 1. - Ht4*pow(i - abs(Gt3), Ht5) + abs(Gt3);		//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}



	//close the PF file
	PF_real.close();

	return PF_real_val;
}
/**********************UDF9**************************/
std::vector<double> UDF9::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double Gt1 = sin(0.5*pi*t_vector[0]);		//Value of G(t1) function - dynamic dependant on time
	double Gt2 = sin(0.5*pi*t_vector[1]);		//Value of G(t2) function - dynamic dependant on time
	double Gt3 = sin(0.5*pi*t_vector[2]);		//Value of G(t3) function - dynamic dependant on time
	double Gt4 = sin(0.5*pi*t_vector[3]);		//Value of G(t4) function - dynamic dependant on time
	double Gt5 = sin(0.5*pi*t_vector[4]);		//Value of G(t5) function - dynamic dependant on time
	double Ht4 = 0.5 + abs(Gt4);				//Value of H(t4) function - dynamic dependant on time
	double Ht5 = 0.5 + abs(Gt5);				//Value of H(t5) function - dynamic dependant on time

	double sizeJ1, sizeJ2, SumJ1, SumJ2;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;
	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - pow(code[0], (0.5*(2. + 3.*(i - 2) / (Vars() - 2) + Gt1))) - Gt2;
		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += pow(y, 2);
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += pow(y, 2);
		}
	}

	std::vector<double> fitness;					//fitness vector - output

	double fitness1_temp = code[0] + abs(Gt3) + (2. / sizeJ1)*SumJ1;			//fitness 1
	double fitness2_temp = 1 - Ht4*pow(code[0], Ht5) + abs(Gt3) + (2. / sizeJ2)*SumJ2;		//fitness 2
	fitness.push_back(fitness1_temp);				//fitness 1
	fitness.push_back(fitness2_temp);			//fitness 2

												//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void UDF9::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -2;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> UDF9::Plot_PF(int indexr, double t, int size)
{

	std::ofstream PF_real;					//ofstream file for the PF saving

	double Gt3 = sin(0.5*pi*t_vector[2]);		//Value of G(t3) function - dynamic dependant on time
	double Gt4 = sin(0.5*pi*t_vector[3]);		//Value of G(t4) function - dynamic dependant on time
	double Gt5 = sin(0.5*pi*t_vector[4]);		//Value of G(t5) function - dynamic dependant on time
	double Ht4 = 0.5 + abs(Gt4);				//Value of H(t4) function - dynamic dependant on time
	double Ht5 = 0.5 + abs(Gt5);				//Value of H(t5) function - dynamic dependant on time

	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

															//calculate the real PF points


															//calculate the real PF points
	for (double i = 0. + abs(Gt3); i <= 1. + abs(Gt3); i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = 1. - Ht4*pow(i - abs(Gt3), Ht5) + abs(Gt3);		//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}



	//close the PF file
	PF_real.close();

	return PF_real_val;
}
/**********************JY1**************************/
std::vector<double> JY1::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	//Function parameters
	double At = 0.05;				//Value of At parameter - dynamic dependant on time, adjust the curvature
	double Wt = 6;					//Value of Wt parameter - dynamic dependant on time, control the number of mixed convex and concave segments
	double Alphat = 1;				//Value of Alpha t parameter - dynamic dependant on time, control the overall shape of POF
	double Betat = 1;				//Value of Beta t parameter - dynamic dependant on time, control the overall shape of POF
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time

	//Calculate temporary values
	double Sum_g_xII = 0.;

	for (int i = 1; i < Vars(); i++)
	{
		Sum_g_xII += pow((code[i] - Gt),2);
	}


	//calculate fitness and push to fitness vector
	std::vector<double> fitness;					//fitness vector - output
													//calculate the fitness and push it to the fitness vector
	//calculate
	double fitness1_temp = (1+Sum_g_xII)*(code[0] + At*sin(Wt*pi*code[0]));			//fitness 1
	double fitness2_temp = (1 + Sum_g_xII)*(1. - code[0] + At*sin(Wt*pi*code[0]));			//fitness 2
	//push
	fitness.push_back(fitness1_temp);				//fitness 1
	fitness.push_back(fitness2_temp);			//fitness 2

												//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void JY1::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> JY1::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	std::ifstream PF_real_input;			//ifstream file of real PF

	//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	double t_temp = 0;

	PF_real_input.open("Input/JY PF/" + this->Name_Show() + "_" + String_Prec(t_temp,1) + "t.dat");
	if (!PF_real_input.good())
	{
		std::cout << "Input/JY PF/" + this->Name_Show() + "_" + String_Prec(t_temp, 1) + "t.dat\n";
		std::cout << "ERROR#12: FUNCTION - PF_READ; FILE NOT OPEN";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}



	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output
	std::vector<double> PF_real_val_temp;	//Vector for storing temp data
	//Read values from file and save in vector and additional file
	double temp_val;			//Currently readed number
	//Read all values
	while (PF_real_input >> temp_val)		//will read as long as there are any values remained
	{
		//push value to temporary vector
		PF_real_val_temp.push_back(temp_val);
		//Push value to file
		PF_real << temp_val;
		//Send values to output vector
		if (PF_real_val_temp.size() == Objs())			//if size of temprary vector is equal to nubmer of objectives
		{
			//push to output vector
			PF_real_val.push_back(PF_real_val_temp);
			//clear the temporary vector
			PF_real_val_temp.clear();
			//end line in output file
			PF_real << std::endl;
		}
		else
			//push space to the file
			PF_real << " ";
	}
	
	//close the PF files
	PF_real.close();
	PF_real_input.close();

	//check the size of output vector
	if (PF_real_val.size() !=  PF_real_size + 1)
	{
		std::cout << "ERROR#12: FUNCTION - PF_READ; SIZE DO NOT MATCH\n";		//ERROR#12: FUNCTION - VECTOR SIZE
		std::cout << "Defined size: " << PF_real_size << "\nReaded values size: " << PF_real_val.size();
		system("pause");
		abort();
	}

	return PF_real_val;
}

/**********************JY2**************************/
std::vector<double> JY2::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	//Function parameters
	double At = 0.05;				//Value of At parameter - dynamic dependant on time, adjust the curvature
	double Wt = floor(6.*sin(0.5*pi*(t-1)));					//Value of Wt parameter - dynamic dependant on time, control the number of mixed convex and concave segments
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time

									//Calculate temporary values
	double Sum_g_xII = 0.;

	for (int i = 1; i < Vars(); i++)
	{
		Sum_g_xII += pow((code[i] - Gt), 2);
	}


	//calculate fitness and push to fitness vector
	std::vector<double> fitness;					//fitness vector - output
													//calculate the fitness and push it to the fitness vector
													//calculate
	double fitness1_temp = (1 + Sum_g_xII)*(code[0] + At*sin(Wt*pi*code[0]));			//fitness 1
	double fitness2_temp = (1 + Sum_g_xII)*(1. - code[0] + At*sin(Wt*pi*code[0]));			//fitness 2
																							//push
	fitness.push_back(fitness1_temp);				//fitness 1
	fitness.push_back(fitness2_temp);			//fitness 2

												//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void JY2::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> JY2::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	std::ifstream PF_real_input;			//ifstream file of real PF

											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	PF_real_input.open("Input/JY PF/" + this->Name_Show() + "_" + String_Prec(t, 1) + "t.dat");

	if (!PF_real_input.good())
	{
		std::cout << "ERROR#12: FUNCTION - PF_READ; FILE NOT OPEN";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}


	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output
	std::vector<double> PF_real_val_temp;	//Vector for storing temp data
											//Read values from file and save in vector and additional file
	double temp_val;			//Currently readed number
								//Read all values
	while (PF_real_input >> temp_val)		//will read as long as there are any values remained
	{
		//push value to temporary vector
		PF_real_val_temp.push_back(temp_val);
		//Push value to file
		PF_real << temp_val;
		//Send values to output vector
		if (PF_real_val_temp.size() == Objs())			//if size of temprary vector is equal to nubmer of objectives
		{
			//push to output vector
			PF_real_val.push_back(PF_real_val_temp);
			//clear the temporary vector
			PF_real_val_temp.clear();
			//end line in output file
			PF_real << std::endl;
		}
		else
			//push space to the file
			PF_real << " ";
	}

	//close the PF files
	PF_real.close();
	PF_real_input.close();

	//check the size of output vector
	if (PF_real_val.size() != PF_real_size + 1)
	{
		std::cout << "ERROR#12: FUNCTION - PF_READ; SIZE DO NOT MATCH\n";		//ERROR#12: FUNCTION - VECTOR SIZE
		std::cout << "Defined size: " << PF_real_size << "\nReaded values size: " << PF_real_val.size();
		system("pause");
		abort();
	}

	return PF_real_val;
}

/**********************JY3**************************/
std::vector<double> JY3::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	//Function parameters
	double At = 0.05;				//Value of At parameter - dynamic dependant on time, adjust the curvature
	double Wt = floor(6.*sin(0.5*pi*(t - 1)));	//Value of Wt parameter - dynamic dependant on time, control the number of mixed convex and concave segments

	//Calculate temporary values
	double Sum_g_xII = 0.;
	double alpha = floor(100 * pow(sin(0.5*pi*t), 2));
	double y_0 = abs(code[0] * sin((2 * alpha + 0.5)*pi*code[0]));


	for (int i = 1; i < Vars(); i++)
	{
		double y_i, y_j;

		y_i = code[i];

		if (i == i)
			y_j = y_0;
		else
			y_j = code[i - 1];
		
		Sum_g_xII += pow(pow(y_i, 2) - y_j, 2);
	}


	//calculate fitness and push to fitness vector
	std::vector<double> fitness;					//fitness vector - output
													//calculate the fitness and push it to the fitness vector
													//calculate
	double fitness1_temp = (1 + Sum_g_xII)*(y_0 + At*sin(Wt*pi*y_0));			//fitness 1
	double fitness2_temp = (1 + Sum_g_xII)*(1. - y_0 + At*sin(Wt*pi*y_0));			//fitness 2
																							//push
	fitness.push_back(fitness1_temp);				//fitness 1
	fitness.push_back(fitness2_temp);			//fitness 2

												//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void JY3::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> JY3::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	std::ifstream PF_real_input;			//ifstream file of real PF

											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	PF_real_input.open("Input/JY PF/" + this->Name_Show() + "_" + String_Prec(t, 1) + "t.dat");

	if (!PF_real_input.good())
	{
		std::cout << "ERROR#12: FUNCTION - PF_READ; FILE NOT OPEN";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}



	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output
	std::vector<double> PF_real_val_temp;	//Vector for storing temp data
											//Read values from file and save in vector and additional file
	double temp_val;			//Currently readed number
								//Read all values
	while (PF_real_input >> temp_val)		//will read as long as there are any values remained
	{
		//push value to temporary vector
		PF_real_val_temp.push_back(temp_val);
		//Push value to file
		PF_real << temp_val;
		//Send values to output vector
		if (PF_real_val_temp.size() == Objs())			//if size of temprary vector is equal to nubmer of objectives
		{
			//push to output vector
			PF_real_val.push_back(PF_real_val_temp);
			//clear the temporary vector
			PF_real_val_temp.clear();
			//end line in output file
			PF_real << std::endl;
		}
		else
			//push space to the file
			PF_real << " ";
	}

	//close the PF files
	PF_real.close();
	PF_real_input.close();

	//check the size of output vector
	if (PF_real_val.size() != PF_real_size + 1)
	{
		std::cout << "ERROR#12: FUNCTION - PF_READ; SIZE DO NOT MATCH\n";		//ERROR#12: FUNCTION - VECTOR SIZE
		std::cout << "Defined size: " << PF_real_size << "\nReaded values size: " << PF_real_val.size();
		system("pause");
		abort();
	}

	return PF_real_val;
}

/**********************JY4**************************/
std::vector<double> JY4::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	//Function parameters
	double At = 0.05;				//Value of At parameter - dynamic dependant on time, adjust the curvature
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double Wt = pow(10, 1 + abs(Gt));					//Value of Wt parameter - dynamic dependant on time, control the number of mixed convex and concave segments

									//Calculate temporary values
	double Sum_g_xII = 0.;

	for (int i = 1; i < Vars(); i++)
	{
		Sum_g_xII += pow((code[i] - Gt), 2);
	}


	//calculate fitness and push to fitness vector
	std::vector<double> fitness;					//fitness vector - output
													//calculate the fitness and push it to the fitness vector
													//calculate
	double fitness1_temp = (1 + Sum_g_xII)*(code[0] + At*sin(Wt*pi*code[0]));			//fitness 1
	double fitness2_temp = (1 + Sum_g_xII)*(1. - code[0] + At*sin(Wt*pi*code[0]));			//fitness 2
																							//push
	fitness.push_back(fitness1_temp);				//fitness 1
	fitness.push_back(fitness2_temp);			//fitness 2

												//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void JY4::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> JY4::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	std::ifstream PF_real_input;			//ifstream file of real PF

											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	PF_real_input.open("Input/JY PF/" + this->Name_Show() + "_" + String_Prec(t, 1) + "t.dat");

	if (!PF_real_input.good())
	{
		std::cout << "ERROR#12: FUNCTION - PF_READ; FILE NOT OPEN";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}



	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output
	std::vector<double> PF_real_val_temp;	//Vector for storing temp data
											//Read values from file and save in vector and additional file
	double temp_val;			//Currently readed number
								//Read all values
	while (PF_real_input >> temp_val)		//will read as long as there are any values remained
	{
		//push value to temporary vector
		PF_real_val_temp.push_back(temp_val);
		//Push value to file
		PF_real << temp_val;
		//Send values to output vector
		if (PF_real_val_temp.size() == Objs())			//if size of temprary vector is equal to nubmer of objectives
		{
			//push to output vector
			PF_real_val.push_back(PF_real_val_temp);
			//clear the temporary vector
			PF_real_val_temp.clear();
			//end line in output file
			PF_real << std::endl;
		}
		else
			//push space to the file
			PF_real << " ";
	}

	//close the PF files
	PF_real.close();
	PF_real_input.close();

	//check the size of output vector
	if (PF_real_val.size() != PF_real_size + 1)
	{
		std::cout << "ERROR#12: FUNCTION - PF_READ; SIZE DO NOT MATCH\n";		//ERROR#12: FUNCTION - VECTOR SIZE
		std::cout << "Defined size: " << PF_real_size << "\nReaded values size: " << PF_real_val.size();
		system("pause");
		abort();
	}

	return PF_real_val;
}

/**********************JY5**************************/
std::vector<double> JY5::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	//Function parameters
	double At = 0.3*sin(0.5*pi*(t-1));				//Value of At parameter - dynamic dependant on time, adjust the curvature
	double Wt = 1;					//Value of Wt parameter - dynamic dependant on time, control the number of mixed convex and concave segments

									//Calculate temporary values
	double Sum_g_xII = 0.;

	for (int i = 1; i < Vars(); i++)
	{
		Sum_g_xII += pow((code[i]), 2);
	}


	//calculate fitness and push to fitness vector
	std::vector<double> fitness;					//fitness vector - output
													//calculate the fitness and push it to the fitness vector
													//calculate
	double fitness1_temp = (1 + Sum_g_xII)*(code[0] + At*sin(Wt*pi*code[0]));			//fitness 1
	double fitness2_temp = (1 + Sum_g_xII)*(1. - code[0] + At*sin(Wt*pi*code[0]));			//fitness 2
																							//push
	fitness.push_back(fitness1_temp);				//fitness 1
	fitness.push_back(fitness2_temp);			//fitness 2

												//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void JY5::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> JY5::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	std::ifstream PF_real_input;			//ifstream file of real PF

											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	PF_real_input.open("Input/JY PF/" + this->Name_Show() + "_" + String_Prec(t, 1) + "t.dat");

	if (!PF_real_input.good())
	{
		std::cout << "ERROR#12: FUNCTION - PF_READ; FILE NOT OPEN";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}



	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output
	std::vector<double> PF_real_val_temp;	//Vector for storing temp data
											//Read values from file and save in vector and additional file
	double temp_val;			//Currently readed number
								//Read all values
	while (PF_real_input >> temp_val)		//will read as long as there are any values remained
	{
		//push value to temporary vector
		PF_real_val_temp.push_back(temp_val);
		//Push value to file
		PF_real << temp_val;
		//Send values to output vector
		if (PF_real_val_temp.size() == Objs())			//if size of temprary vector is equal to nubmer of objectives
		{
			//push to output vector
			PF_real_val.push_back(PF_real_val_temp);
			//clear the temporary vector
			PF_real_val_temp.clear();
			//end line in output file
			PF_real << std::endl;
		}
		else
			//push space to the file
			PF_real << " ";
	}

	//close the PF files
	PF_real.close();
	PF_real_input.close();

	//check the size of output vector
	if (PF_real_val.size() != PF_real_size + 1)
	{
		std::cout << "ERROR#12: FUNCTION - PF_READ; SIZE DO NOT MATCH\n";		//ERROR#12: FUNCTION - VECTOR SIZE
		std::cout << "Defined size: " << PF_real_size << "\nReaded values size: " << PF_real_val.size();
		system("pause");
		abort();
	}

	return PF_real_val;
}

/**********************JY6**************************/
std::vector<double> JY6::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	//Function parameters
	double At = 0.1;				//Value of At parameter - dynamic dependant on time, adjust the curvature
	double Wt = 3;					//Value of Wt parameter - dynamic dependant on time, control the number of mixed convex and concave segments
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double Kt = 2 * floor(10 * abs(Gt));

									//Calculate temporary values
	double Sum_g_xII = 0.;

	for (int i = 1; i < Vars(); i++)
	{
		double y_i = code[i] - Gt;
		Sum_g_xII += 4 * pow(y_i, 2) - cos(Kt*pi*y_i) + 1;
	}


	//calculate fitness and push to fitness vector
	std::vector<double> fitness;					//fitness vector - output
													//calculate the fitness and push it to the fitness vector
													//calculate
	double fitness1_temp = (1 + Sum_g_xII)*(code[0] + At*sin(Wt*pi*code[0]));			//fitness 1
	double fitness2_temp = (1 + Sum_g_xII)*(1. - code[0] + At*sin(Wt*pi*code[0]));			//fitness 2
																							//push
	fitness.push_back(fitness1_temp);				//fitness 1
	fitness.push_back(fitness2_temp);			//fitness 2

												//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void JY6::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> JY6::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	std::ifstream PF_real_input;			//ifstream file of real PF

											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	double t_temp = 0;

	PF_real_input.open("Input/JY PF/" + this->Name_Show() + "_" + String_Prec(t_temp, 1) + "t.dat");
	if (!PF_real_input.good())
	{
		std::cout << "Input/JY PF/" + this->Name_Show() + "_" + String_Prec(t_temp, 1) + "t.dat\n";
		std::cout << "ERROR#12: FUNCTION - PF_READ; FILE NOT OPEN";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}


	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output
	std::vector<double> PF_real_val_temp;	//Vector for storing temp data
											//Read values from file and save in vector and additional file
	double temp_val;			//Currently readed number
								//Read all values
	while (PF_real_input >> temp_val)		//will read as long as there are any values remained
	{
		//push value to temporary vector
		PF_real_val_temp.push_back(temp_val);
		//Push value to file
		PF_real << temp_val;
		//Send values to output vector
		if (PF_real_val_temp.size() == Objs())			//if size of temprary vector is equal to nubmer of objectives
		{
			//push to output vector
			PF_real_val.push_back(PF_real_val_temp);
			//clear the temporary vector
			PF_real_val_temp.clear();
			//end line in output file
			PF_real << std::endl;
		}
		else
			//push space to the file
			PF_real << " ";
	}

	//close the PF files
	PF_real.close();
	PF_real_input.close();

	//check the size of output vector
	if (PF_real_val.size() != PF_real_size + 1)
	{
		std::cout << "ERROR#12: FUNCTION - PF_READ; SIZE DO NOT MATCH\n";		//ERROR#12: FUNCTION - VECTOR SIZE
		std::cout << "Defined size: " << PF_real_size << "\nReaded values size: " << PF_real_val.size();
		system("pause");
		abort();
	}

	return PF_real_val;
}

/**********************JY7**************************/
std::vector<double> JY7::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	//Function parameters
	double At = 0.1;				//Value of At parameter - dynamic dependant on time, adjust the curvature
	double Wt = 3;					//Value of Wt parameter - dynamic dependant on time, control the number of mixed convex and concave segments
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double Alphat = 0.2+2.8*abs(Gt);				//Value of Alpha t parameter - dynamic dependant on time, control the overall shape of POF
	double Betat = Alphat;				//Value of Beta t parameter - dynamic dependant on time, control the overall shape of POF


									//Calculate temporary values
	double Sum_g_xII = 0.;

	for (int i = 1; i < Vars(); i++)
	{
		double y_i = code[i] - Gt;
		Sum_g_xII += pow(y_i, 2) - 10 * cos(2 * pi*y_i) + 10;
	}


	//calculate fitness and push to fitness vector
	std::vector<double> fitness;					//fitness vector - output
													//calculate the fitness and push it to the fitness vector
													//calculate
	double fitness1_temp = (1 + Sum_g_xII)*pow(code[0] + At*sin(Wt*pi*code[0]),Alphat);			//fitness 1
	double fitness2_temp = (1 + Sum_g_xII)*pow(1. - code[0] + At*sin(Wt*pi*code[0]),Betat);			//fitness 2
																							//push
	fitness.push_back(fitness1_temp);				//fitness 1
	fitness.push_back(fitness2_temp);			//fitness 2

												//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void JY7::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> JY7::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	std::ifstream PF_real_input;			//ifstream file of real PF

											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	PF_real_input.open("Input/JY PF/" + this->Name_Show() + "_" + String_Prec(t, 1) + "t.dat");

	if (!PF_real_input.good())
	{
		std::cout << "ERROR#12: FUNCTION - PF_READ; FILE NOT OPEN";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}



	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output
	std::vector<double> PF_real_val_temp;	//Vector for storing temp data
											//Read values from file and save in vector and additional file
	double temp_val;			//Currently readed number
								//Read all values
	while (PF_real_input >> temp_val)		//will read as long as there are any values remained
	{
		//push value to temporary vector
		PF_real_val_temp.push_back(temp_val);
		//Push value to file
		PF_real << temp_val;
		//Send values to output vector
		if (PF_real_val_temp.size() == Objs())			//if size of temprary vector is equal to nubmer of objectives
		{
			//push to output vector
			PF_real_val.push_back(PF_real_val_temp);
			//clear the temporary vector
			PF_real_val_temp.clear();
			//end line in output file
			PF_real << std::endl;
		}
		else
			//push space to the file
			PF_real << " ";
	}

	//close the PF files
	PF_real.close();
	PF_real_input.close();

	//check the size of output vector
	if (PF_real_val.size() != PF_real_size + 1)
	{
		std::cout << "ERROR#12: FUNCTION - PF_READ; SIZE DO NOT MATCH\n";		//ERROR#12: FUNCTION - VECTOR SIZE
		std::cout << "Defined size: " << PF_real_size << "\nReaded values size: " << PF_real_val.size();
		system("pause");
		abort();
	}

	return PF_real_val;
}

/**********************JY8**************************/
std::vector<double> JY8::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	//Function parameters
	double At = 0.05;				//Value of At parameter - dynamic dependant on time, adjust the curvature
	double Wt = 6;					//Value of Wt parameter - dynamic dependant on time, control the number of mixed convex and concave segments
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double Betat = 10-9.8*abs(Gt);				//Value of Beta t parameter - dynamic dependant on time, control the overall shape of POF
	double Alphat = 2./Betat;				//Value of Alpha t parameter - dynamic dependant on time, control the overall shape of POF


									//Calculate temporary values
	double Sum_g_xII = 0.;

	for (int i = 1; i < Vars(); i++)
	{
		Sum_g_xII += pow(code[i], 2);
	}


	//calculate fitness and push to fitness vector
	std::vector<double> fitness;					//fitness vector - output
													//calculate the fitness and push it to the fitness vector
													//calculate
	double fitness1_temp = (1 + Sum_g_xII)*pow(code[0] + At*sin(Wt*pi*code[0]), Alphat);			//fitness 1
	double fitness2_temp = (1 + Sum_g_xII)*pow(1. - code[0] + At*sin(Wt*pi*code[0]), Betat);			//fitness 2
																							//push
	fitness.push_back(fitness1_temp);				//fitness 1
	fitness.push_back(fitness2_temp);			//fitness 2

												//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void JY8::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> JY8::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	std::ifstream PF_real_input;			//ifstream file of real PF

											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	PF_real_input.open("Input/JY PF/" + this->Name_Show() + "_" + String_Prec(t, 1) + "t.dat");

	if (!PF_real_input.good())
	{
		std::cout << "ERROR#12: FUNCTION - PF_READ; FILE NOT OPEN";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}



	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output
	std::vector<double> PF_real_val_temp;	//Vector for storing temp data
											//Read values from file and save in vector and additional file
	double temp_val;			//Currently readed number
								//Read all values
	while (PF_real_input >> temp_val)		//will read as long as there are any values remained
	{
		//push value to temporary vector
		PF_real_val_temp.push_back(temp_val);
		//Push value to file
		PF_real << temp_val;
		//Send values to output vector
		if (PF_real_val_temp.size() == Objs())			//if size of temprary vector is equal to nubmer of objectives
		{
			//push to output vector
			PF_real_val.push_back(PF_real_val_temp);
			//clear the temporary vector
			PF_real_val_temp.clear();
			//end line in output file
			PF_real << std::endl;
		}
		else
			//push space to the file
			PF_real << " ";
	}

	//close the PF files
	PF_real.close();
	PF_real_input.close();

	//check the size of output vector
	if (PF_real_val.size() != PF_real_size + 1)
	{
		std::cout << "ERROR#12: FUNCTION - PF_READ; SIZE DO NOT MATCH\n";		//ERROR#12: FUNCTION - VECTOR SIZE
		std::cout << "Defined size: " << PF_real_size << "\nReaded values size: " << PF_real_val.size();
		system("pause");
		abort();
	}

	return PF_real_val;
}

/**********************JY9**************************/
int sigma;
extern int dyn_tau;
std::vector<double> JY9::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	//Function parameters
	double At = 0.05;				//Value of At parameter - dynamic dependant on time, adjust the curvature
	double Wt = floor(6*pow(sin(0.5*pi*(t-1)),sigma));					//Value of Wt parameter - dynamic dependant on time, control the number of mixed convex and concave segments
	double Alphat = 1;				//Value of Alpha t parameter - dynamic dependant on time, control the overall shape of POF
	double Betat = 1;				//Value of Beta t parameter - dynamic dependant on time, control the overall shape of POF
	double Gt = abs(sin(0.5*pi*t));		//Value of G(t) function - dynamic dependant on time

									//Calculate temporary values
	double Sum_g_xII = 0.;

	for (int i = 1; i < Vars(); i++)
	{
		Sum_g_xII += pow(code[i] - Gt + sigma, 2);
	}


	//calculate fitness and push to fitness vector
	std::vector<double> fitness;					//fitness vector - output
													//calculate the fitness and push it to the fitness vector
													//calculate
	double fitness1_temp = (1 + Sum_g_xII)*(code[0] + At*sin(Wt*pi*code[0]));			//fitness 1
	double fitness2_temp = (1 + Sum_g_xII)*(1. - code[0] + At*sin(Wt*pi*code[0]));			//fitness 2
																							//push
	fitness.push_back(fitness1_temp);				//fitness 1
	fitness.push_back(fitness2_temp);			//fitness 2

												//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void JY9::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> JY9::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	std::ifstream PF_real_input;			//ifstream file of real PF
	int ro_t = 5;					//severity of change
	sigma = (int)floor((double)dyn_tau / (T_dyn*ro_t)) % 3;
											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	PF_real_input.open("Input/JY PF/" + this->Name_Show() + "_" + std::to_string(sigma) + "sigma_" + String_Prec(t, 1) + "t.dat");

	if (!PF_real_input.good())
	{
		std::cout << "ERROR#12: FUNCTION - PF_READ; FILE NOT OPEN";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}



	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output
	std::vector<double> PF_real_val_temp;	//Vector for storing temp data
											//Read values from file and save in vector and additional file
	double temp_val;			//Currently readed number
								//Read all values
	while (PF_real_input >> temp_val)		//will read as long as there are any values remained
	{
		//push value to temporary vector
		PF_real_val_temp.push_back(temp_val);
		//Push value to file
		PF_real << temp_val;
		//Send values to output vector
		if (PF_real_val_temp.size() == Objs())			//if size of temprary vector is equal to nubmer of objectives
		{
			//push to output vector
			PF_real_val.push_back(PF_real_val_temp);
			//clear the temporary vector
			PF_real_val_temp.clear();
			//end line in output file
			PF_real << std::endl;
		}
		else
			//push space to the file
			PF_real << " ";
	}

	//close the PF files
	PF_real.close();
	PF_real_input.close();

	//check the size of output vector
	if (PF_real_val.size() != PF_real_size + 1)
	{
		std::cout << "ERROR#12: FUNCTION - PF_READ; SIZE DO NOT MATCH\n";		//ERROR#12: FUNCTION - VECTOR SIZE
		std::cout << "Defined size: " << PF_real_size << "\nReaded values size: " << PF_real_val.size();
		system("pause");
		abort();
	}

	return PF_real_val;
}

/**********************JY10**************************/
std::vector<double> JY10::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	//Function parameters
	double At = 0.05;				//Value of At parameter - dynamic dependant on time, adjust the curvature
	double Wt = 6;					//Value of Wt parameter - dynamic dependant on time, control the number of mixed convex and concave segments
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	
	double Alphat = 1+sigma*Gt;				//Value of Alpha t parameter - dynamic dependant on time, control the overall shape of POF
	double Betat = Alphat;				//Value of Beta t parameter - dynamic dependant on time, control the overall shape of POF

									//Calculate temporary values
	double Sum_g_xII = 0.;

	for (int i = 1; i < Vars(); i++)
	{
		Sum_g_xII += pow(code[i] - Gt + sigma, 2);
	}


	//calculate fitness and push to fitness vector
	std::vector<double> fitness;					//fitness vector - output
													//calculate the fitness and push it to the fitness vector
													//calculate
	double fitness1_temp = (1 + Sum_g_xII)*pow(code[0] + At*sin(Wt*pi*code[0]), Alphat);			//fitness 1
	double fitness2_temp = (1 + Sum_g_xII)*pow(1. - code[0] + At*sin(Wt*pi*code[0]), Betat);			//fitness 2
																							//push
	fitness.push_back(fitness1_temp);				//fitness 1
	fitness.push_back(fitness2_temp);			//fitness 2

												//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void JY10::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> JY10::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	std::ifstream PF_real_input;			//ifstream file of real PF
	int ro_t = 5;					//severity of change
	sigma = (int)floor((double)dyn_tau / (T_dyn*ro_t) + Random_I(1, 3)) % 3;
											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	PF_real_input.open("Input/JY PF/" + this->Name_Show() + "_" + std::to_string(sigma) + "sigma_" + String_Prec(t, 1) + "t.dat");

	if (!PF_real_input.good())
	{
		std::cout << "ERROR#12: FUNCTION - PF_READ; FILE NOT OPEN";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}



	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output
	std::vector<double> PF_real_val_temp;	//Vector for storing temp data
											//Read values from file and save in vector and additional file
	double temp_val;			//Currently readed number
								//Read all values
	while (PF_real_input >> temp_val)		//will read as long as there are any values remained
	{
		//push value to temporary vector
		PF_real_val_temp.push_back(temp_val);
		//Push value to file
		PF_real << temp_val;
		//Send values to output vector
		if (PF_real_val_temp.size() == Objs())			//if size of temprary vector is equal to nubmer of objectives
		{
			//push to output vector
			PF_real_val.push_back(PF_real_val_temp);
			//clear the temporary vector
			PF_real_val_temp.clear();
			//end line in output file
			PF_real << std::endl;
		}
		else
			//push space to the file
			PF_real << " ";
	}

	//close the PF files
	PF_real.close();
	PF_real_input.close();

	//check the size of output vector
	if (PF_real_val.size() != PF_real_size + 1)
	{
		std::cout << "ERROR#12: FUNCTION - PF_READ; SIZE DO NOT MATCH\n";		//ERROR#12: FUNCTION - VECTOR SIZE
		std::cout << "Defined size: " << PF_real_size << "\nReaded values size: " << PF_real_val.size();
		system("pause");
		abort();
	}

	return PF_real_val;
}

/**********************FDA1**************************/
std::vector<double> FDA1::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double sum_x = 0;						//temp values for the fitness calculation
	

	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time

									//calculate temp values
	for (int i = 1; i < Vars(); i++)
	{
		sum_x += pow(code[i] - Gt, 2);
	}

	double g_x = 1 + sum_x;
	double h = 1 - sqrt(code[0] / g_x);

	std::vector<double> fitness;									//fitness vector - output
	//calculate the fitness
	double fit1 = code[0];
	double fit2 = g_x * h;
	//push to the fitness vector
	fitness.push_back(fit1);			//fitness 1
	fitness.push_back(fit2);	//fitness 2

																		//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;

}

void FDA1::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> FDA1::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

	

	//calculate the real PF points
	for (double i = 0. ; i <= 1.; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1. - sqrt(i));		//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************FDA2**************************/
std::vector<double> FDA2::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double sum_x2 = 0;						//temp values for the fitness calculation
	double sum_x3 = 0;

	double Ht = 0.75 + 0.7*sin(0.5*pi*t); // Value of H(t) function - dynamic dependant on time
									//calculate temp values
	for (int i = 1; i < 16; i++)
	{
		sum_x2 += pow(code[i], 2);
	}

	for (int i = 16; i < Vars(); i++)
	{
		sum_x3 += pow(code[i] - Ht, 2);
	}


	double g_x = 1 + sum_x2;
	double h = 1 - pow((code[0] / g_x),(Ht+sum_x3));

	std::vector<double> fitness;									//fitness vector - output
																	//calculate the fitness
	double fit1 = code[0];
	double fit2 = g_x * h;
	//push to the fitness vector
	fitness.push_back(fit1);			//fitness 1
	fitness.push_back(fit2);	//fitness 2

								//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;

}

void FDA2::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> FDA2::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

	double Ht = 0.75 + 0.7*sin(0.5*pi*t); // Value of H(t) function - dynamic dependant on time

	//calculate the real PF points
	for (double i = 0.; i <= 1.; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double mod;

		if (Ht <= 1)
			mod = Ht;
		else
			mod = Ht + 15 * (pow(Ht - 1, 2));
		double temp = (1. - pow(i,mod));		//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************FDA3**************************/
std::vector<double> FDA3::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double sum_x1 = 0;						//temp values for the fitness calculation
	double sum_x2 = 0;

	double Gt = abs(sin(0.5*pi*t));		//Value of G(t) function - dynamic dependant on time
	double Ft = pow(10,2*(sin(0.5*pi*t)));		//Value of F(t) function - dynamic dependant on time
									//calculate temp values
	for (int i = 0; i < 1; i++)
	{
		sum_x1 += pow(code[i], Ft);
	}

	for (int i = 1; i < Vars(); i++)
	{
		sum_x2 += pow(code[i] - Gt, 2);
	}

	double g_x = 1 + Gt + sum_x2;
	double h = 1 - sqrt(sum_x1 / g_x);

	std::vector<double> fitness;									//fitness vector - output
																	//calculate the fitness
	double fit1 = sum_x1;
	double fit2 = g_x * h;
	//push to the fitness vector
	fitness.push_back(fit1);			//fitness 1
	fitness.push_back(fit2);	//fitness 2

								//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;

}

void FDA3::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

	

	//set boundaries for the rest of variables
	for (int i = 0; i < 1; i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> FDA3::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

	double Gt = abs(sin(0.5*pi*t));
	double Ft = pow(10, 2 * (sin(0.5*pi*t)));		//Value of F(t) function - dynamic dependant on time

										  //calculate the real PF points
	for (double i = 0.; i <= 1.; i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp1 = i;

		
		double temp2 = (1. -sqrt(temp1/(1.+Gt)))*(1+Gt);		//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(temp1);
		temp_vect.push_back(temp2);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << temp1 << " ";
			PF_real << temp2;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************FDA4**************************/
std::vector<double> FDA4::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double sum_x2 = 0;						//temp values for the fitness calculation
	double Multip_x = 1;

	double Gt = abs(sin(0.5*pi*t));		//Value of G(t) function - dynamic dependant on time

									//calculate temp values
	for (int i = 1; i < Vars(); i++)
	{
		sum_x2 += pow(code[i] - Gt, 2);
	}

	Multip_x = cos(code[0] * pi / 2);

	double g_x2 = sum_x2;

	std::vector<double> fitness;									//fitness vector - output
																	//calculate the fitness
	double fit1 = (1+g_x2)*Multip_x;
	double fit2 = (1 + g_x2)*sin(code[0]*pi/2);
	//push to the fitness vector
	fitness.push_back(fit1);			//fitness 1
	fitness.push_back(fit2);	//fitness 2

								//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;

}

void FDA4::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> FDA4::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

	double Gt = sin(0.5*pi*t);

	//calculate the real PF points
	for (double i = 0. + abs(Gt); i <= 1. + abs(Gt); i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1. - i + 2 * abs(Gt));		//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************FDA5**************************/
std::vector<double> FDA5::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double sum_x2 = 0;						//temp values for the fitness calculation
	double Multip_y = 1;

	double Gt = abs(sin(0.5*pi*t));		//Value of G(t) function - dynamic dependant on time
	double Ft = 1+pow(sin(0.5*pi*t),4);		//Value of F(t) function - dynamic dependant on time
										//calculate temp values
	for (int i = 1; i < Vars(); i++)
	{
		sum_x2 += pow(code[i] - Gt, 2);
	}

	double y1 = pow(code[0], Ft);

	Multip_y = cos(y1 * pi / 2);

	double g_x2 = Gt + sum_x2;

	std::vector<double> fitness;									//fitness vector - output
																	//calculate the fitness
	double fit1 = (1 + g_x2)*Multip_y;
	double fit2 = (1 + g_x2)*sin(y1 * pi / 2);
	//push to the fitness vector
	fitness.push_back(fit1);			//fitness 1
	fitness.push_back(fit2);	//fitness 2

								//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;

}

void FDA5::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> FDA5::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

	double Gt = sin(0.5*pi*t);

	//calculate the real PF points
	for (double i = 0. + abs(Gt); i <= 1. + abs(Gt); i += 1.0 / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1. - i + 2 * abs(Gt));		//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/********************************************
DYNAMIC CONSTRAINED FUNCTIONS
*********************************************/

/**********************CDF1**************************/
std::vector<double> CDF1::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;						//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;

	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time

									//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - Gt - pow(code[0], (0.5*(1.0 + 3.0*(i - 2) / (double)(Vars() - 2))));
		if (i % 2 == 1)
		{
			sizeJ1++;
			SumJ1 += pow(y, 2);
		}
		else
		{
			sizeJ2++;
			SumJ2 += pow(y, 2);
		}
	}

	std::vector<double> fitness;									//fitness vector - output

																	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + abs(Gt) + (2.0 / sizeJ1)*SumJ1);			//fitness 1
	fitness.push_back(1. - code[0] + abs(Gt) + (2.0 / sizeJ2)*SumJ2);	//fitness 2

																		//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;

}

void CDF1::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -2;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> CDF1::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	short N = 10;		//CEC 09
											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

	double Gt = sin(0.5*pi*t);
	/*
	//calculate the real PF points
	for (double i = 0; i <= 2 * N; i++)
	{
		//calculate the PF value
		double temp = (double)i / (double)(2 * N) + abs(Gt);
		double temp2 = (1. - ((double)i / (double)(2 * N)) + abs(Gt));				//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(temp);
		temp_vect.push_back(temp2);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << temp << " ";
			PF_real << temp2;
			PF_real << std::endl;
		}
	}*/
	//Calculate the rest of the points
	for (double i = 0 + abs(Gt); i <= 1. + abs(Gt); i += 1. / (double)(size))
	{
		//calculate the PF valued
		double temp = 1. - i + 2* abs(Gt);		//PF value

		double cons_value = i + temp - 2 * abs(Gt) - abs(sin((double)N*pi*(i - temp + 1.))) - 1.;

		if (cons_value < cons_check_param)
		{
			//double i221 = 21;
			//	std::cout << i << " " << cons_value << std::endl;
			continue;

		}
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}
std::vector<double> CDF1::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	double cons_value;					//Value for constrain check
	short N = 10;		// CEC '09
	short a = 1;		// CEC '09
	std::vector<double> cons_val_temp;	//Values of the constrains - output
	double Gt = sin(0.5*pi*t);
										//Calculate constrain
	cons_value = fit[0] + fit[1] -2*abs(Gt)- (double)a*abs(sin((double)N*pi*(fit[0] - fit[1] + 1.))) - 1.;

	//push value to the vector
	cons_val_temp.push_back(cons_value);

	//return vector
	return cons_val_temp;

}

/**********************CDF2**************************/
std::vector<double> CDF2::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;
	double Gt = sin(0.5*pi*t);

	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y;

		if (i % 2 == 1)
		{
			y = code[i - 1] - 0.8*code[0] * cos(6.*pi*code[0] + ((double)i*pi / (double)Vars()));
			sizeJ1 += 1;
			SumJ1 += pow(y - abs(Gt), 2);
		}
		else
		{
			y = code[i - 1] - 0.8*code[0] * sin(6.*pi*code[0] + ((double)i*pi / (double)Vars()));
			sizeJ2 += 1;
			SumJ2 += pow(y - abs(Gt), 2);
		}
	}

	std::vector<double> fitness;					//fitness vector - output

													//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + SumJ1 + abs(Gt));						//fitness 1
	fitness.push_back(pow(1.0 - code[0], 2) + SumJ2 + abs(Gt));		//fitness 2

															//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void CDF2::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -2;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 40,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 1.5,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> CDF2::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	double Gt = sin(0.5*pi*t);
											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

															//calculate the real PF points
	for (double i = 0.; i <= 1.; i += 1.0 / (double)(size - 1))
	{
		double temp;				//PF value

									//calculate the PF value
		if (i <= 0.5 )
			temp = pow(1. - i, 2);
		else if (i <= 0.75 )
			temp = 0.5*(1. - i);
		else
			temp = 0.25*sqrt(1. - i);
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i + abs(Gt));
		temp_vect.push_back(temp + abs(Gt));

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i + abs(Gt) << " ";
			PF_real << temp + abs(Gt);
			PF_real << std::endl;
		}
	}
	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> CDF2::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	double cons_value;					//Value for the constrain check
	double cons_value2;					//2nd Value for the constrain check
	std::vector<double> cons_val_temp;	//Values of the constrains - output
	double Gt = sin(0.5*pi*t);


	double temp, temp2;				//temporary values
									//Calculate temporary values
	temp = 0.5*(1. - code[0]) - pow(1. - code[0], 2);
	temp2 = 0.25 * sqrt(1. - code[0]) - 0.5*(1. - code[0]);

	//Calculate constrain values
	cons_value = code[1] - 0.8*code[0] * sin(6.*pi*code[0] + (2.*pi / (double)Vars())) - sgn(temp)*sqrt(abs(temp)) - abs(Gt);
	cons_value2 = code[3] - 0.8*code[0] * sin(6.*pi*code[0] + (4.*pi / (double)Vars())) - sgn(temp2)*sqrt(abs(temp2)) - abs(Gt);

	//push values to the vector
	cons_val_temp.push_back(cons_value);
	cons_val_temp.push_back(cons_value2);

	//return vector
	return cons_val_temp;
}

/**********************CDF3**************************/
std::vector<double> CDF3::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;						//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;

	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double Mt = 0.5 + abs(Gt);		//Value of M(t) function - dynamic dependant on time
	double Kt = ceil(Vars()*Gt);					//Value of K(t) function - dynamic dependant on time	

													//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - sin((6. * pi * code[0]) + (((double)i + Kt) * pi) / (double)Vars());
		if (i % 2 == 1)
		{
			sizeJ1++;
			SumJ1 += pow(y, 2);
		}
		else
		{
			sizeJ2++;
			SumJ2 += pow(y, 2);
		}
	}

	std::vector<double> fitness;									//fitness vector - output

																	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + (2.0 / sizeJ1)*SumJ1);						//fitness 1
	fitness.push_back(1. - (Mt*pow(code[0], Mt)) + (2.0 / sizeJ2)*SumJ2);	//fitness 2

																			//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void CDF3::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> CDF3::Plot_PF(int indexr, double t, int size)
{
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double Mt = 0.5 + abs(Gt);					//Value of M(t) function - dynamic dependant on time
	double Ht = Mt;					//Value of H(t) function - dynamic dependant on time
	short N = 2;
	std::ofstream PF_real;					//ofstream file for the PF saving

											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output


	std::vector<double> temp_vect1;	
	temp_vect1.push_back(0.);
	temp_vect1.push_back(1.);
	//push the fitness vector to the PF vector
	PF_real_val.push_back(temp_vect1);

	//save the fitness values to the file
	if (indexr >= 0)
	{
		PF_real << 0. << " ";
		PF_real << 1.;
		PF_real << std::endl;
	}

	//Calculate the rest of the points
	for (double i = 0; i <= 1.; i += 1. / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1. - (Mt*pow(i, Ht)));		//PF value

		if (((temp + sqrt(i)) - sin((double)N*pi*(sqrt(i) - temp + 1.)) - 1.) < cons_check_param)
			continue;
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}


	/*double i1min, i1max, i2min, i2max;

	i1min = pow((2. / (2. * N)), 1 / Ht);
	i1max = pow((3. / (2. * N)), 1 / Ht);
	i2min = pow((3. / (2. * N)), 1 / Ht);
	i2max = pow((4. / (2. * N)), 1 / Ht);

	//calculate the step of the i value for the loop
	double i_step = ((i1max - i1min) + (i2max - i2min)) / ((double)size - 2.0);

	//Calculate the rest of the points
	for (double i = i1min; i <= i1max; i += i_step)
	{
		//calculate the PF value
		double temp = (1. - (Mt*pow(i, Ht)));		//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	for (double i = i2min; i <= i2max; i += i_step)
	{
		//calculate the PF value
		double temp = (1. - (Mt*pow(i, Ht)));		//PF value
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}*/


	//close the PF file
	PF_real.close();


	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> CDF3::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	double cons_value;					//Value for constrain check
	double temp;							//Temporary value
	short N = 2;		// CEC '09
	short a = 1;		// CEC '09
	std::vector<double> cons_val_temp;	//Values of the constrains - output



										//Calculate t
	temp = fit[1] + sqrt(fit[0]) - (double)a*sin((double)N*pi*(sqrt(fit[0]) - fit[1] + 1.)) - 1.;
	
	//temp = fit[1] + Mt*pow(fit[0],Ht) - (double)a*sin((double)N*pi*(Mt*pow(fit[0], Ht) - fit[1] + 1.)) - 1.;

	//Calculate constrain
	cons_value = temp / (1. + exp(4.*abs(temp)));

	//push value to the vector
	cons_val_temp.push_back(cons_value);

	//return vector
	return cons_val_temp;
}

/**********************CDF4**************************/
std::vector<double> CDF4::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;						//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;

	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double Mt = 0.5 + abs(Gt);		//Value of M(t) function - dynamic dependant on time

									//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - Gt - pow(code[0], (0.5*(2.0 + 3.0*(i - 2) / (double)(Vars() - 2)) + Gt));
		if (i % 2 == 1)
		{
			sizeJ1++;
			SumJ1 += pow(y, 2);
		}
		else
		{
			sizeJ2++;
			SumJ2 += pow(y, 2);
		}
	}

	std::vector<double> fitness;									//fitness vector - output

																	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + (2.0 / sizeJ1)*SumJ1);			//fitness 1
	fitness.push_back(1. - Mt*pow(code[0], Mt) + (2.0 / sizeJ2)*SumJ2);	//fitness 2

																		//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;

}

void CDF4::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> CDF4::Plot_PF(int indexr, double t, int size)
{
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double Mt = 0.5 + abs(Gt);					//Value of M(t) function - dynamic dependant on time
	double Ht = Mt;					//Value of H(t) function - dynamic dependant on time
	short N = 2;
	std::ofstream PF_real;					//ofstream file for the PF saving

											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output


	std::vector<double> temp_vect1;
	temp_vect1.push_back(0.);
	temp_vect1.push_back(1.);
	//push the fitness vector to the PF vector
	PF_real_val.push_back(temp_vect1);

	//save the fitness values to the file
	if (indexr >= 0)
	{
		PF_real << 0. << " ";
		PF_real << 1.;
		PF_real << std::endl;
	}
	
	//Calculate the rest of the points
	for (double i = 0; i <= 1.; i += 1. / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1. - (Mt*pow(i, Ht)));		//PF value
		
		if ((temp + pow(i, 0.5) - (double)1 * sin((double)N*pi*(pow(i, 0.5) - temp + 1.)) - 1.) < cons_check_param)
			continue;
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();


	//close the PF file
	PF_real.close();

	return PF_real_val;
}
std::vector<double> CDF4::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	double cons_value;					//Value for constrain check
	short N = 2;		// CEC '09
	short a = 1;		// CEC '09
	std::vector<double> cons_val_temp;	//Values of the constrains - output

	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double Mt = 0.5 + abs(Gt);					//Value of M(t) function - dynamic dependant on time
	double Ht = Mt;					//Value of H(t) function - dynamic dependant on time

										//Calculate constrain
	cons_value = fit[1] + pow(fit[0], 0.5) - (double)a*sin((double)N*pi*(pow(fit[0], 0.5) - fit[1] + 1.)) - 1.;

	//push value to the vector
	cons_val_temp.push_back(cons_value);

	//return vector
	return cons_val_temp;
}


/**********************CDF5**************************/
std::vector<double> CDF5::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time

	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y;
		if (i % 2 == 1)
			y = code[i - 1] - 0.8*code[0] * cos(6.*pi*code[0] + ((double)i*pi / (double)Vars()))- Gt;
		else
			y = code[i - 1] - 0.8*code[0] * sin(6.*pi*code[0] + ((double)i*pi / (double)Vars())) - Gt;
		double h;

		if (i == 2)
		{
			if (y < (3. / 2.*(1. - (sqrt(2.) / 2.))))
				h = abs(y);
			else
				h = 0.125 + pow(y - 1., 2);
		}
		else
			h = 2. * pow(y, 2) - cos(4.*pi*y) + 1.;

		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += h;
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += h;
		}
	}

	std::vector<double> fitness;					//fitness vector - output

													//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + SumJ1 + abs(Gt));				//fitness 1
	fitness.push_back(1.0 - code[0] + SumJ2 + abs(Gt));		//fitness 2

													//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void CDF5::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -2;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 35,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 35,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> CDF5::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time

											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

															//calculate the real PF points
	for (double i = 0; i <= 1.; i += 1.0 / (double)(size - 1))
	{
		double temp;				//PF value

									//calculate the PF value
		if (i <= 0.5 )
			temp = 1. - i;
		else if (i <= 0.75)
			temp = -0.5*i + (3. / 4.);
		else
			temp = 1. - i + 0.125;
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i + abs(Gt));
		temp_vect.push_back(temp + abs(Gt));

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i + abs(Gt) << " ";
			PF_real << temp + abs(Gt);
			PF_real << std::endl;
		}
	}
	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> CDF5::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	double cons_value;					//Value for constrain check
	std::vector<double> cons_val_temp;	//Values of the constrains - output
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time

										//Calculate constrain value
	cons_value = code[1] - 0.8*code[0] * sin(6.*pi*code[0] + (2.*pi / (double)Vars())) - 0.5*code[0] + 0.25 - Gt;

	//push value to the vector
	cons_val_temp.push_back(cons_value);

	//return vector
	return cons_val_temp;
}
/**********************CDF6**************************/
std::vector<double> CDF6::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double Gt1 = sin(0.5*pi*t_vector[0]);		//Value of G(t1) function - dynamic dependant on time
	double Gt2 = sin(0.5*pi*t_vector[1]);		//Value of G(t2) function - dynamic dependant on time
	double Gt3 = sin(0.5*pi*t_vector[2]);		//Value of G(t3) function - dynamic dependant on time
	double Gt4 = sin(0.5*pi*t_vector[3]);		//Value of G(t4) function - dynamic dependant on time
	double Gt5 = sin(0.5*pi*t_vector[4]);		//Value of G(t5) function - dynamic dependant on time
	double Kt1 = ceil(Vars()*Gt1);				//Value of K(t1) function - dynamic dependant on time
	double Ht4 = 0.5 + abs(Gt4);				//Value of H(t4) function - dynamic dependant on time
	double Ht5 = 0.5 + abs(Gt5);				//Value of H(t5) function - dynamic dependant on time

	double sizeJ1, sizeJ2, SumJ1, SumJ2, MultipJ1, MultipJ2;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;
	MultipJ1 = MultipJ2 = 1;
	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - sin(6 * pi*code[0] + ((double)i + Kt1)*pi / Vars()) - Gt2;
		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += pow(y, 2);
			//MultipJ1 *= cos((20. * y*pi) / sqrt(i));
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += pow(y, 2);
			//MultipJ2 *= cos((20. * y*pi) / sqrt(i));
		}
	}

	/*for (int i = 2; i <= Vars(); i++)
	{
	double y = code[i - 1] - sin(6 * pi*code[0] + ((double)i + Kt1)*pi / Vars()) - Gt2;
	if (i % 2 == 1)
	{
	sizeJ1 += 1;
	SumJ1 += pow(y, 2);
	}
	else
	{
	sizeJ2 += 1;
	SumJ2 += pow(y, 2);
	}
	}*/

	std::vector<double> fitness;					//fitness vector - output

	double fitness1_temp = code[0] + abs(Gt3) + (2.0 / sizeJ1)*(SumJ1);			//fitness 1
	double fitness2_temp = 1 - Ht4*pow(code[0], Ht5) + abs(Gt3) + (2.0 / sizeJ2)*(SumJ2);		//fitness 2
	fitness.push_back(fitness1_temp);				//fitness 1
	fitness.push_back(fitness2_temp);			//fitness 2

												//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void CDF6::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -2;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 10,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> CDF6::Plot_PF(int indexr, double t, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving

	double Gt2 = sin(0.5*pi*t_vector[1]);		//Value of G(t2) function - dynamic dependant on time
	double Gt3 = sin(0.5*pi*t_vector[2]);		//Value of G(t3) function - dynamic dependant on time
	double Gt4 = sin(0.5*pi*t_vector[3]);		//Value of G(t4) function - dynamic dependant on time
	double Gt5 = sin(0.5*pi*t_vector[4]);		//Value of G(t5) function - dynamic dependant on time
	double Ht4 = 0.5 + abs(Gt4);				//Value of H(t4) function - dynamic dependant on time
	double Ht5 = 0.5 + abs(Gt5);				//Value of H(t5) function - dynamic dependant on time

	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

															//calculate the real PF points

															//calculate the real PF points
	for (double i = 0.; i <= 1.; i += 1.0 / (double)(size - 1))
	{

		//calculate the PF value
		double temp = 1. - Ht4*pow(i, Ht5) + abs(Gt3);		//PF value
															//temp = 1. - Ht4*pow(i - abs(Gt3), Ht5) + abs(Gt3);

		if ((temp + Ht4*pow(i + abs(Gt3), Ht5) - sin(2.*pi*(Ht4*pow(i + abs(Gt3), Ht5) - temp + 1.)) - 1.) < cons_check_param)
			continue;

		std::vector<double> temp_vect;			//temporary vector of fitness
												//push values to the fitness vector
		temp_vect.push_back(i + abs(Gt3));
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i + abs(Gt3) << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> CDF6::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double time)
{
	double cons_value;					//Value for constrain check
	double t;							//Temporary value
	short N = 2;

	double Gt3 = sin(0.5*pi*t_vector[2]);		//Value of G(t3) function - dynamic dependant on time
	double Gt4 = sin(0.5*pi*t_vector[3]);		//Value of G(t4) function - dynamic dependant on time
	double Gt5 = sin(0.5*pi*t_vector[4]);		//Value of G(t5) function - dynamic dependant on time
	double Ht4 = 0.5 + abs(Gt4);				//Value of H(t4) function - dynamic dependant on time
	double Ht5 = 0.5 + abs(Gt5);				//Value of H(t5) function - dynamic dependant on time

	std::vector<double> cons_val_temp;	//Values of the constrains - output

										//Calculate t
	cons_value = fit[1] + Ht4*pow(fit[0], Ht5) - sin((double)N*pi*(Ht4*pow(fit[0], Ht5) - fit[1] + 1.)) - 1.;

	//push value to the vector
	cons_val_temp.push_back(cons_value);

	//return vector
	return cons_val_temp;
}

/**********************CDF7**************************/
std::vector<double> CDF7::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double Ht = 0.5 + abs(Gt);					//Value of H(t) function - dynamic dependant on time
	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y;
		if (i % 2 == 1)
		{
			y = code[i - 1] - sin((6. * pi * code[0]) + ((double)i * pi) / (double)Vars());
			sizeJ1++;
			SumJ1 += pow(y, 2);
		}
		else
		{
			y = code[i - 1] - cos((6. * pi * code[0]) + ((double)i * pi) / (double)Vars());
			sizeJ2++;
			SumJ2 += pow(y, 2);
		}
	}

	std::vector<double> fitness;									//fitness vector - output

																	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + (2.0 / sizeJ1)*SumJ1);				//fitness 1
	fitness.push_back(1. - pow(code[0],Ht) + (2.0 / sizeJ2)*SumJ2);	//fitness 2

																	//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void CDF7::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{6,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 6,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> CDF7::Plot_PF(int indexr, double t, int size)
{
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double Ht = 0.5 + abs(Gt);					//Value of M(t) function - dynamic dependant on time
	short N = 2;
	std::ofstream PF_real;					//ofstream file for the PF saving

											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

	//Calculate the rest of the points
	for (double i = 0; i <= 1.; i += 1. / (double)(size - 1))
	{
		//calculate the PF value
		double temp = (1. - (pow(i, Ht)));		//PF value
		
		double temp_const = (temp + sqrt(i) - sin((double)N*pi*(sqrt(i) - temp + 1.)) - 1.);
		double cons_value = temp_const / (1. + exp(4.*abs(temp_const)));
		
		if (cons_value < cons_check_param)
			continue;
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();


	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> CDF7::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	double cons_value;					//Value for constrain check
	double temp;							//Temporary value
	short N = 2;		// CEC '09
	short a = 1;		// CEC '09
	std::vector<double> cons_val_temp;	//Values of the constrains - output

										//Calculate t
	temp = fit[1] + sqrt(fit[0]) - (double)a*sin((double)N*pi*(sqrt(fit[0]) - fit[1] + 1.)) - 1.;

	//Calculate constrain
	cons_value = temp / (1. + exp(4.*abs(temp)));

	//push value to the vector
	cons_val_temp.push_back(cons_value);

	//return vector
	return cons_val_temp;
}

/**********************CDF8**************************/
std::vector<double> CDF8::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
												//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - sin(6.*pi*code[0] + ((double)i*pi / (double)Vars())) - Gt;
		double h;

		if (i == 2)
		{
			if (y < (3. / 2.*(1. - (sqrt(2.) / 2.))))
				h = abs(y);
			else
				h = 0.125 + pow(y - 1., 2);
		}
		else
			h = pow(y, 2);

		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += h;
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += h;
		}
	}

	std::vector<double> fitness;					//fitness vector - output

													//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + SumJ1);				//fitness 1
	fitness.push_back(1.0 - code[0] + SumJ2);		//fitness 2

													//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void CDF8::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -2;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 25,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 25,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> CDF8::Plot_PF(int indexr, double t, int size)
{
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double Ht = 0.5 + abs(Gt);					//Value of M(t) function - dynamic dependant on time
	short N = 2;
	std::ofstream PF_real;					//ofstream file for the PF saving

											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

	//calculate the real PF points
	for (double i = 0; i <= 1.; i += 1.0 / (double)(size - 1))
	{
		double temp;				//PF value

									//calculate the PF value
		if (i <= 0.5)
			temp = 1. - i;
		else if (i <= 0.75)
			temp = -0.5*i + (3. / 4.);
		else
			temp = 1. - i + 0.125;
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);

		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}
	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> CDF8::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	double cons_value;					//Value for constrain check
	double temp;							//Temporary value
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	std::vector<double> cons_val_temp;	//Values of the constrains - output

										//Calculate temp
	temp = code[1] - sin(6.*pi*code[0] + (2.*pi / (double)Vars())) - 0.5*code[0] + 0.25 - Gt;

	//Calculate constrain
	cons_value = temp / (1. + pow(e, 4.*abs(temp)));

	//push value to the vector
	cons_val_temp.push_back(cons_value);

	//return vector
	return cons_val_temp;
}


/**********************CDF9**************************/
std::vector<double> CDF9::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;						//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;

	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - pow(code[0], (0.5*(1. + (3.*(i - 2) / (double)(Vars() - 2)))));
		if (i % 2 == 1)
		{
			sizeJ1++;
			SumJ1 += pow(y, 2);
		}
		else
		{
			sizeJ2++;
			SumJ2 += pow(y, 2);
		}
	}

	std::vector<double> fitness;									//fitness vector - output

																	//calculate the fitness and push it to the fitness vector
	fitness.push_back(code[0] + (2.0 / sizeJ1)*SumJ1);				//fitness 1
	fitness.push_back(1. - code[0] + (2.0 / sizeJ2)*SumJ2);	//fitness 2

															//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void CDF9::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 6,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 6,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> CDF9::Plot_PF(int indexr, double t, int size)
{
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	short N = 10;
	std::ofstream PF_real;					//ofstream file for the PF saving

											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output


	//Calculate the rest of the points
	for (double i = 0; i <= 1.; i += 1. / (double)(size))
	{
		//calculate the PF value
		double temp = 1. - i;		//PF value

		double cons_value = i + temp - abs(sin((double)N*pi*(i - temp + 1.))) - 1. + abs(Gt);

		if (cons_value < cons_check_param)
		{
			//double i221 = 21;
		//	std::cout << i << " " << cons_value << std::endl;
			continue;
			
		}
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();


	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> CDF9::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double cons_value;					//Value for constrain check
	short N = 10;		// CEC '09
	short a = 1;		// CEC '09
	std::vector<double> cons_val_temp;	//Values of the constrains - output

										//Calculate constrain
	cons_value = fit[0] + fit[1] - abs(sin((double)N*pi*(fit[0] - fit[1] + 1.))) - 1. + abs(Gt);

	//push value to the vector
	cons_val_temp.push_back(cons_value);

	//return vector
	return cons_val_temp;
}

/**********************CDF10**************************/
std::vector<double> CDF10::Fitness_C(const std::vector<double> & code, int tau, double t)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2, MultipJ1, MultipJ2;			//temp values for the fitness calculation
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;
	MultipJ1 = MultipJ2 = 1;
	//calculate temp values
	for (int i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - sin(6.*pi*code[0] + ((double)i*pi / (double)Vars()));

		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += pow(y, 2);
			//MultipJ1 *= cos((20. * y*pi) / sqrt(i));
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += pow(y, 2);
			//MultipJ2 *= cos((20. * y*pi) / sqrt(i));
		}
	}

	std::vector<double> fitness;					//fitness vector - output

													//calculate the fitness and push it to the fitness vector
	//fitness.push_back(code[0] + (2.0 / sizeJ1)*(4.0 * SumJ1 - 2.0 * MultipJ1 + 2.));				//fitness 1
//	fitness.push_back(1.0 - pow(code[0], 2) + (2.0 / sizeJ2)*(4.0 * SumJ2 - 2.0 * MultipJ2 + 2.));		//fitness 2
	fitness.push_back(code[0] + (2.0 / sizeJ1)*(SumJ1));				//fitness 1
	fitness.push_back(1.0 - pow(code[0], 2) + (2.0 / sizeJ2)*(SumJ2));		//fitness 2
																										//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void CDF10::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

																		//set the boundaries for the first variable - if different than others
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);

	//set boundaries for the rest of variables
	for (int i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -2;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 40,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 40,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> CDF10::Plot_PF(int indexr, double t, int size)
{
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double Ht = 0.5 + abs(Gt);					//Value of M(t) function - dynamic dependant on time
	short N = 2;
	std::ofstream PF_real;					//ofstream file for the PF saving

											//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output


															//Calculate the rest of the points
	for (double i = 0; i <= 1.; i += 1. / (double)(size - 1))
	{
		//calculate the PF value
		double temp = 1. - pow(i,2);		//PF value

		double cons_value = pow(i, 2) + temp - sin((double)N*pi*(pow(i, 2) - temp + 1. + Gt)) - 1. ;

		if (cons_value < cons_check_param)
			continue;
		std::vector<double> temp_vect;			//temporary vector of fitness

												//push values to the fitness vector
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		//push the fitness vector to the PF vector
		PF_real_val.push_back(temp_vect);

		//save the fitness values to the file
		if (indexr >= 0)
		{
			PF_real << i << " ";
			PF_real << temp;
			PF_real << std::endl;
		}
	}

	//close the PF file
	PF_real.close();


	//close the PF file
	PF_real.close();

	return PF_real_val;
}

std::vector<double> CDF10::Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t)
{
	double Gt = sin(0.5*pi*t);		//Value of G(t) function - dynamic dependant on time
	double cons_value;					//Value for constrain check
	short N = 2;		// CEC '09
	short a = 1;		// CEC '09
	std::vector<double> cons_val_temp;	//Values of the constrains - output

										//Calculate constrain
	cons_value = fit[1] + pow(fit[0], 2) - (double)a*sin((double)N*pi*(pow(fit[0], 2) - fit[1] + 1. + Gt)) - 1.;

	//push value to the vector
	cons_val_temp.push_back(cons_value);

	//return vector
	return cons_val_temp;
}




std::vector<double> key_score;
std::vector<double> key_std;
/**********************GEN**************************/
std::vector<double> GEN::Fitness_C(const std::vector<double> & code)
{
	double fit1 = 0, fit2 = 0, fit3 = 0;
	short vars = Vars();
	for (int ix = 0; ix < vars; ix++)
	{
		double temp_code = code[ix];
		if (temp_code >= 0.5)
		{
			fit1 += key_score[ix];
			fit2 += key_std[ix];
			fit3++;
		}
	}
	if (fit2 != 0)
	{
		fit1 /= 2.3955;
		fit2 /= 6.0567;
	}

	fit1 = 1 - fit1;
	fit3 = (vars - fit3) / vars;

	




	std::vector<double> fitness;					//fitness vector - output

	//calculate the fitness and push it to the fitness vector
	fitness.push_back(fit1);				//fitness 1
	fitness.push_back(fit2);		//fitness 2
	fitness.push_back(fit3);		//fitness 3
													//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void GEN::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes

		

	//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 2,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 1000,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> GEN::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving

	std::ifstream key_file;
	std::string line;
	key_score.clear();
	key_std.clear();
	//open the key file
	key_file.open("input/gene/CADD_norm_gdi_proc_score_refined_0.01.csv");
	int ix = 0;
	while (std::getline(key_file, line))
	{
		if (ix == 0)
		{
			ix++;
			continue;
		}
		int index = 0;
		std::stringstream ss(line);

		std::string val;
		while (std::getline(ss, val, ','))
		{
			if (index == 11)
			{
				std::stringstream ss2(val);
				double value;
				ss2 >> value;

				key_std.push_back(value);
			}
			else if (index == 13)
			{
				std::stringstream ss2(val);
				double value;
				ss2 >> value;

				key_score.push_back(value);
			}
			index++;
		}
	}

	if (key_score.size() != Vars() || key_std.size() != Vars())
		abort();


	//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

	//calculate the real PF points

	double i = 0;
	double temp = 0;				//PF value

	std::vector<double> temp_vect;			//temporary vector of fitness

											//push values to the fitness vector
	temp_vect.push_back(i);
	temp_vect.push_back(temp);

	//push the fitness vector to the PF vector
	PF_real_val.push_back(temp_vect);

	//save the fitness values to the file
	if (indexr >= 0)
	{
		PF_real << i << " ";
		PF_real << temp;
		PF_real << std::endl;
	}
	//close the PF file
	PF_real.close();

	return PF_real_val;
}

/**********************GEN_pat**************************/
std::vector<double> GEN_pat::Fitness_C(const std::vector<double> & code)
{
	double fit1 = 0, fit2 = 0, fit3 = 0;
	short vars = Vars();
	for (int ix = 0; ix < vars; ix++)
	{
		double temp_code = code[ix];
		if (temp_code >= 0.5)
		{
			fit1 += key_score[ix];
			fit2 += key_std[ix];
			fit3++;
		}
	}
	if (fit2 != 0)
	{
		fit1 /= 2.3955;
		fit2 /= 6.4159;
	}

	fit1 = 1 - fit1;
	fit3 = (vars - fit3) / vars;






	std::vector<double> fitness;					//fitness vector - output

	//calculate the fitness and push it to the fitness vector
	fitness.push_back(fit1);				//fitness 1
	fitness.push_back(fit2);		//fitness 2
	fitness.push_back(fit3);		//fitness 3
													//check if the fitness has a proper size
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	return fitness;
}

void GEN_pat::Bound_Set()
{
	std::vector<STRUCTURES::boundaries> bound;							//vector of boundaries structure for the variables
	std::vector<STRUCTURES::boundaries> max_min;						//vector of boundaries structure for the fitness - only for the plotting purposes



	//set boundaries for the rest of variables
	for (int i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}

	//check if the amount of boudaries is correct
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set boundaries for the fitness
	STRUCTURES::boundaries t1{ 2,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 1000,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	//assign boundary vectors to class vectors
	function::Bound_Set(bound);
	function::Max_Min_Fit_Set(max_min);
}
std::vector<std::vector<double>> GEN_pat::Plot_PF(int indexr, int size)
{
	std::ofstream PF_real;					//ofstream file for the PF saving

	std::ifstream key_file;
	std::string line;
	key_score.clear();
	key_std.clear();
	//open the key file
	key_file.open("input/gene/CADD_norm_gdi_proc_score_refined_0.01.csv");
	int ix = 0;
	while (std::getline(key_file, line))
	{
		if (ix == 0)
		{
			ix++;
			continue;
		}
		int index = 0;
		std::stringstream ss(line);

		std::string val;
		while (std::getline(ss, val, ','))
		{
			if (index == 11)
			{
				std::stringstream ss2(val);
				double value;
				ss2 >> value;

				key_std.push_back(value);
			}
			else if (index == 13)
			{
				std::stringstream ss2(val);
				double value;
				ss2 >> value;

				key_score.push_back(value);
			}
			index++;
		}
	}

	if (key_score.size() != Vars() || key_std.size() != Vars())
		abort();


	//check if the index is correct
	if (indexr >= 0)
		//open a new PF file
		PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(0., 1) + ".x1");

	std::vector<std::vector<double>> PF_real_val;			//Vector for PF data storage - output

	//calculate the real PF points

	double i = 0;
	double temp = 0;				//PF value

	std::vector<double> temp_vect;			//temporary vector of fitness

											//push values to the fitness vector
	temp_vect.push_back(i);
	temp_vect.push_back(temp);

	//push the fitness vector to the PF vector
	PF_real_val.push_back(temp_vect);

	//save the fitness values to the file
	if (indexr >= 0)
	{
		PF_real << i << " ";
		PF_real << temp;
		PF_real << std::endl;
	}
	//close the PF file
	PF_real.close();

	return PF_real_val;
}

