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

/*The SBX crossover method is based on K. Deb and R. Bhushan Agrawal, “Simulated Binary Crossover for Continuous Search Space,” Complex Syst., vol. 9, pp. 115–148, 1995.*/


#pragma once
//****************************************
//				CROSSOVER header
//	Storage of different crossover types
//			and basic crossover class
//****************************************
#ifndef CROSSOVER_H
#define CROSSOVER_H

#include "GA_Functions.h"
#include "Fit_Functions.h"
#include <cmath>
#include "Define.h"

/**Crossover class template**/
template <typename tname>
class crossover
{
private:
	short ccode;					//code of the crossover type
	std::string name;				//name of the crossover type
public:
	/*
	*Default normal constructor*
	@param cc crossover type code
	@param na crossover type name
	*/
	crossover(short cc, std::string na) { ccode = cc; name = na; };
	/**Default constuctor, throw error - crossover cannot be empty**/
	crossover() { std::cout << "ERROR#08: CROSSOVER - CREATION"; system("pause"); abort(); }
	~crossover() {};
	/**Returning name of the current crossover type**/
	std::string Name_Show() { return name; };
	/*
	*Empty crossover - make crossover for the given population*
	@param selected vector containing selected individuals
	@param gapara GA_Parameters class object
	@param fcode function class object 
	*/
	virtual std::vector<tname>  Crossover(const std::vector<tname> & selected, GA_parameters & gapara, function & fcode) const { std::vector<tname> a; abort(); return a; };				//Take 2 individuals and return 2 indi
};

/****************************************
		Crossover #1 - REAL VALUES 
****************************************/
template <typename tname>
class crossover_1 : public crossover<tname>
{
public: 
	/**Default constructor**/
	crossover_1() : crossover(1, "Crossover_SBX_Changed") {};
	~crossover_1() {};
	/*
	*make crossover for the given population*
	@param selected vector containing selected individuals
	@param gapara GA_Parameters class object
	@param fcode function class object
	*/
	std::vector<tname> Crossover(const std::vector<tname> & selected, GA_parameters & gapara, function & fcode) const;
};

/*
*make crossover for the given population*
@param selected vector containing selected individuals
@param gapara GA_Parameters class object
@param fcode function class object
*/
template <typename tname>
std::vector<tname> crossover_1<tname>::Crossover(const std::vector<tname> & selected, GA_parameters & gapara, function & fcode) const
{
	std::vector<tname> cross_indi;		//Output crossover vector
	int psize = selected.size();		//size of the given selection vector
	bool psize_odd = false;

	//Check if given selection vector is odd - crossover can be done only on the even vector
	if (psize % 2 == 1)
	{
		psize_odd = true;
		//decreasing size of selection vector to make it even
		psize -= 1;
	}
	std::vector<double> upper_b_v,lower_b_v;	//boundaries storage vectors
	//assigning boundaries to storage vectors
	for (int i = 0; i < fcode.Vars(); i++)
	{
		upper_b_v.push_back(fcode.Bound(i, "upper"));
		lower_b_v.push_back(fcode.Bound(i, "lower"));
	}
	//actual crossover loop
	for (int j = 0; j < psize; j += 2)
	{
		std::vector<double> parent1_v = selected[j].Code_Show(); //parent 1 storage vector
		std::vector<double> parent2_v = selected[j+1].Code_Show(); //parent 2 storage vector
		float cross_r = Random_F();					//random 
		
		//Checking if crossover occur
		if (cross_r <= gapara.Cross_Prob())
		{
			std::vector<double> code1;	//code of 1 parent storage vector
			std::vector<double> code2;	//code of 2 parent storage vector

			//Loop for number of variables
			for (int i = 0; i < fcode.Vars(); i++)
			{
				double parent1; //code variable of 1 parent - safe storage
				double parent2;	//code variable of 2 parent - safe storage
				double y2;		//code variable of 2, bigger value parent 
				double y1;		//code variable of 1, lower value parent
				double beta, alpha, betaq; //crossover parameters
				
				//copying variables from storage vector
				parent1 = parent1_v[i];
				parent2 = parent2_v[i];

				//copying boundaries from storage vector
				double upper_b = upper_b_v[i];
				double lower_b = lower_b_v[i];

				//checking if variables are not the same
				if (abs(parent1 - parent2) > pow(10, -14))
				{
					//checking which one is greater - for beta calculaton
					if (parent2 > parent1)					
					{
						y2 = parent2;
						y1 = parent1;
					}
					else
					{
						y2 = parent1;
						y1 = parent2;
					}

					//beta and alpha calculation
					if ((y1 - lower_b) > (upper_b - y2))
						beta = 1.0 + 2.0 * ((upper_b - y2) / (y2 - y1));
					else
						beta = 1.0 + 2.0 * ((y1 - lower_b) / (y2 - y1));
					beta = 1.0 / beta;
					alpha = 2.0 - std::pow(beta, -(gapara.Di_C() + 1.0));

					//randomly choosing degree of change of the values 
					double rnd = Random();
					
					//betaq calculation - in what extent values will change
					if (rnd < 1. / alpha)
					{
						alpha *= rnd;
						betaq = std::pow(alpha, 1.0 / (gapara.Di_C() + 1));
					}
					else
					{
						alpha *= rnd;
						alpha = 1. / (2. - alpha);
						betaq = std::pow(alpha, 1.0 / (gapara.Di_C() + 1));
					}
					double c1, c2;
					c1 = 0.5*((y1 + y2) - betaq*(y2 - y1));
					c2 = 0.5*((y1 + y2) + betaq*(y2 - y1));

					if (c1 < lower_b)
						c1 = lower_b;
					else if (c1 > upper_b)
						c1 = upper_b;
					if (c2 < lower_b)
						c2 = lower_b;
					else if (c2 > upper_b)
						c2 = upper_b;
					
					//final code calculation and sending to storage vector
					code1.push_back(c1);
					code2.push_back(c2);
				}
				//if are the same, crossover do not occur
				else
				{	/*In original was:*/
					//cross_indi[1] = parent2;
					//cross_indi[0] = parent1;
					code1.push_back(parent2);
					code2.push_back(parent1);
				}
			}
			//creation of new individuals with calculated code
			tname c1{ code1, fcode };
			tname c2{ code2, fcode };

			//sending new individuals to output vector
			cross_indi.push_back(c1);
			cross_indi.push_back(c2);
		}
		//if not offsprings are genetically identical to parents
		else
		{
			cross_indi.push_back(selected[j]);
			cross_indi.push_back(selected[j+1]);
		}

		if (TGM == true)
		{
			//Random roll
			float roll = Random_F();
			//Create the output vectors
			std::vector<std::vector<double>> out_1 = selected[j].TGM_fitness;
			std::vector<std::vector<double>> out_2 = selected[j+1].TGM_fitness;

			//perform changes (remove the last position, and insert the new one at the beginning)
			out_1.insert(out_1.begin(), std::vector<double>(fcode.Objs(), 0));
			out_2.insert(out_2.begin(), std::vector<double>(fcode.Objs(), 0));
			out_1.pop_back();
			out_2.pop_back();
			
			//assign the vectors to the offsprings
			if (roll < 0.5f)
			{
				cross_indi[j].TGM_fitness = out_1;
				cross_indi[j+1].TGM_fitness = out_2;
			}
			else
			{
				cross_indi[j].TGM_fitness = out_2;
				cross_indi[j + 1].TGM_fitness = out_1;
			}
			//check the correct sizes
			if (cross_indi[j].TGM_fitness.size() != TGM_size + 1 || cross_indi[j + 1].TGM_fitness.size() != TGM_size + 1)
				abort();
		}
	} //end of crossover loop
	//add the last individual (in the case of odd size vector)
	if (psize_odd == true)
	{
		//if is, last individual is added to new vector without change
		cross_indi.push_back(selected.back());
		if (TGM == true)
		{
			cross_indi.back().TGM_fitness.insert(cross_indi.back().TGM_fitness.begin(), std::vector<double>(fcode.Objs(), 0));
			cross_indi.back().TGM_fitness.pop_back();

			if (cross_indi.back().TGM_fitness.size() != TGM_size + 1)
				abort();
		}

	}
	//Checking if output vector have the right size
	if (cross_indi.size() != selected.size())
	{
		std::cout << "ERROR#10: CROSSOVER - VECTOR SIZE"; //ERROR#10: CROSSOVER - VECTOR SIZE
		system("pause");
		abort();
	}

	//end of function
	return cross_indi;
};


/****************************************
Crossover #1b - REAL VALUES SBX
****************************************/
template <typename tname>
class crossover_1b : public crossover<tname>
{
public:
	/**Default constructor**/
	crossover_1b() : crossover(1, "Crossover_SBX") {};
	~crossover_1b() {};
	/*
	*make crossover for the given population*
	@param selected vector containing selected individuals
	@param gapara GA_Parameters class object
	@param fcode function class object
	*/
	std::vector<tname> Crossover(const std::vector<tname> & selected, GA_parameters & gapara, function & fcode) const;
};

/*
*make crossover for the given population*
@param selected vector containing selected individuals
@param gapara GA_Parameters class object
@param fcode function class object
*/
template <typename tname>
std::vector<tname> crossover_1b<tname>::Crossover(const std::vector<tname> & selected, GA_parameters & gapara, function & fcode) const
{
	std::vector<tname> cross_indi;		//Output crossover vector
	int psize = selected.size();		//size of the given selection vector
	bool psize_odd = false;

	//Check if given selection vector is odd - crossover can be done only on the even vector
	if (psize % 2 == 1)
	{
		psize_odd = true;
		//decreasing size of selection vector to make it even
		psize -= 1;
	}
	std::vector<double> upper_b_v, lower_b_v;	//boundaries storage vectors
												//assigning boundaries to storage vectors
	for (int i = 0; i < fcode.Vars(); i++)
	{
		upper_b_v.push_back(fcode.Bound(i, "upper"));
		lower_b_v.push_back(fcode.Bound(i, "lower"));
	}
	//actual crossover loop
	for (int j = 0; j < psize; j += 2)
	{
		std::vector<double> parent1_v = selected[j].Code_Show(); //parent 1 storage vector
		std::vector<double> parent2_v = selected[j + 1].Code_Show(); //parent 2 storage vector
		float cross_r = Random_F();					//random 

													//Checking if crossover occur
		if (cross_r <= gapara.Cross_Prob())
		{
			std::vector<double> code1;	//code of 1 parent storage vector
			std::vector<double> code2;	//code of 2 parent storage vector

										//Loop for number of variables
			for (int i = 0; i < fcode.Vars(); i++)
			{
				if (Random() <= 0.5)
				{
					double parent1; //code variable of 1 parent - safe storage
					double parent2;	//code variable of 2 parent - safe storage
					double y2;		//code variable of 2, bigger value parent 
					double y1;		//code variable of 1, lower value parent
					double beta, alpha, betaq; //crossover parameters
					double c1, c2;
					//copying variables from storage vector
					parent1 = parent1_v[i];
					parent2 = parent2_v[i];

					//copying boundaries from storage vector
					double upper_b = upper_b_v[i];
					double lower_b = lower_b_v[i];

					//checking if variables are not the same
					if (abs(parent1 - parent2) > pow(10, -14))
					{
						//checking which one is greater - for beta calculaton
						if (parent2 > parent1)
						{
							y2 = parent2;
							y1 = parent1;
						}
						else
						{
							y2 = parent1;
							y1 = parent2;
						}
						double rnd = Random();


						//beta and alpha calculation
						beta = 1.0 + 2.0 * ((y1 - lower_b) / (y2 - y1));
						alpha = 2.0 - std::pow(beta, -(gapara.Di_C() + 1.0));
						if (rnd < (1.0 / alpha))
						{
							betaq = std::pow((rnd*alpha), (1.0 / (gapara.Di_C() + 1)));
						}
						else
						{
							betaq = std::pow((1.0 / (2.0 - rnd*alpha)), (1.0 / (gapara.Di_C() + 1)));
						}
						c1 = 0.5*((y1 + y2) - betaq*(y2 - y1));

						beta = 1.0 + 2.0 * ((upper_b - y2) / (y2 - y1));
						alpha = 2.0 - std::pow(beta, -(gapara.Di_C() + 1.0));
						if (rnd < (1.0 / alpha))
						{
							betaq = std::pow((rnd*alpha), (1.0 / (gapara.Di_C() + 1)));
						}
						else
						{
							betaq = std::pow((1.0 / (2.0 - rnd*alpha)), (1.0 / (gapara.Di_C() + 1)));
						}
						c2 = 0.5*((y1 + y2) + betaq*(y2 - y1));

						if (c1 < lower_b)
							c1 = lower_b;
						else if (c1 > upper_b)
							c1 = upper_b;
						if (c2 < lower_b)
							c2 = lower_b;
						else if (c2 > upper_b)
							c2 = upper_b;

						//final code calculation and sending to storage vector
						if (Random() <= 0.5)
						{
							code1.push_back(c2);
							code2.push_back(c1);
						}
						else
						{
							code2.push_back(c2);
							code1.push_back(c1);
						}
					}
					//if are the same, crossover do not occur
					else
					{
						code1.push_back(parent1);
						code2.push_back(parent2);
					}
				}
				else
				{
					code1.push_back(parent1_v[i]);
					code2.push_back(parent2_v[i]);
				}
			}
			//creation of new individuals with calculated code
			tname c1{ code1, fcode };
			tname c2{ code2, fcode };

			//sending new individuals to output vector
			cross_indi.push_back(c1);
			cross_indi.push_back(c2);
		}
		//if not offsprings are genetically identical to parents
		else
		{
			cross_indi.push_back(selected[j]);
			cross_indi[j].Crowd_Dist_Set(0.0);
			cross_indi[j].Rank_Set(0);
			cross_indi.push_back(selected[j + 1]);
			cross_indi[j + 1].Crowd_Dist_Set(0.0);
			cross_indi[j + 1].Rank_Set(0);
		}
		if (TGM == true)
		{
			//Random roll
			float roll = Random_F();
			//Create the output vectors
			std::vector<std::vector<double>> out_1 = selected[j].TGM_fitness;
			std::vector<std::vector<double>> out_2 = selected[j + 1].TGM_fitness;

			//perform changes (remove the last position, and insert the new one at the beginning)
			out_1.insert(out_1.begin(), std::vector<double>(fcode.Objs(), 0));
			out_2.insert(out_2.begin(), std::vector<double>(fcode.Objs(), 0));
			out_1.pop_back();
			out_2.pop_back();

			//assign the vectors to the offsprings
			if (roll < 0.5f)
			{
				cross_indi[j].TGM_fitness = out_1;
				cross_indi[j + 1].TGM_fitness = out_2;
			}
			else
			{
				cross_indi[j].TGM_fitness = out_2;
				cross_indi[j + 1].TGM_fitness = out_1;
			}
			//check the correct sizes
			if (cross_indi[j].TGM_fitness.size() != TGM_size + 1 || cross_indi[j + 1].TGM_fitness.size() != TGM_size + 1)
				abort();
		}
	} //end of crossover loop

	  //add the last individual (in the case of odd size vector)
	if (psize_odd == true)
	{
		//if is, last individual is added to new vector without change
		cross_indi.push_back(selected.back());
		if (TGM == true)
		{
			cross_indi.back().TGM_fitness.insert(cross_indi.back().TGM_fitness.begin(), std::vector<double>(fcode.Objs(), 0));
			cross_indi.back().TGM_fitness.pop_back();

			if (cross_indi.back().TGM_fitness.size() != TGM_size + 1)
				abort();
		}

	}
	  //Checking if output vector have the right size
	if (cross_indi.size() != selected.size())
	{
		std::cout << "ERROR#10: CROSSOVER - VECTOR SIZE"; //ERROR#10: CROSSOVER - VECTOR SIZE
		system("pause");
		abort();
	}

	//end of function
	return cross_indi;
};

/****************************************
Crossover #B1 - Uniform Crossover
****************************************/
template <typename tname>
class crossover_B1 : public crossover<tname>
{
public:
	/**Default constructor**/
	crossover_B1() : crossover(1, "Uniform Binary Crossover") {};
	~crossover_B1() {};
	/*
	*make crossover for the given population*
	@param selected vector containing selected individuals
	@param gapara GA_Parameters class object
	@param fcode function class object
	*/
	std::vector<tname> Crossover(const std::vector<tname> & selected, GA_parameters & gapara, function & fcode) const;
};

/*
*make crossover for the given population*
@param selected vector containing selected individuals
@param gapara GA_Parameters class object
@param fcode function class object
*/
template <typename tname>
std::vector<tname> crossover_B1<tname>::Crossover(const std::vector<tname> & selected, GA_parameters & gapara, function & fcode) const
{
	std::vector<tname> cross_indi;		//Output crossover vector
	int psize = selected.size();		//size of the given selection vector
	bool psize_odd = false;

	//Check if given selection vector is odd - crossover can be done only on the even vector
	if (psize % 2 == 1)
	{
		psize_odd = true;
		//decreasing size of selection vector to make it even
		psize -= 1;
	}
													//actual crossover loop
	for (int j = 0; j < psize; j += 2)
	{
		std::vector<bool> parent1_v = selected[j].Code_Bin_Show(); //parent 1 storage vector
		std::vector<bool> parent2_v = selected[j + 1].Code_Bin_Show(); //parent 2 storage vector
		float cross_r = Random_F();					//random 
		int bin_string_size = parent1_v.size();
		//Checking if crossover occur
		if (cross_r <= gapara.Cross_Prob())
		{
			std::vector<bool> code1;	//code of 1 parent storage vector
			std::vector<bool> code2;	//code of 2 parent storage vector

										//Loop for the binary string
			for (int i = 0; i < bin_string_size; i++)
			{
				//get random value
				float rand = Random_F();

				if (rand <= 0.5)
				{
					//first offspring gets one gene from the first parent and second from the second one
					code1.push_back(parent1_v[i]);
					code2.push_back(parent2_v[i]);
				}
				else
				{
					//first offspring gets one gene from the second parent and second from the first one
					code1.push_back(parent2_v[i]);
					code2.push_back(parent1_v[i]);
				}
			}
			//create the empty real values vector
			std::vector<double> code_temp = std::vector<double>(fcode.Vars(), 0.);

			//creation of new individuals with calculated code
			tname c1{ code1,code_temp , fcode };
			tname c2{ code2,code_temp , fcode };

			//sending new individuals to output vector
			cross_indi.push_back(c1);
			cross_indi.push_back(c2);
		}
		//if not offsprings are genetically identical to parents
		else
		{
			cross_indi.push_back(selected[j]);
			cross_indi[j].Crowd_Dist_Set(0.0);
			cross_indi[j].Rank_Set(0);
			cross_indi.push_back(selected[j + 1]);
			cross_indi[j + 1].Crowd_Dist_Set(0.0);
			cross_indi[j + 1].Rank_Set(0);
		}

		if (TGM == true)
		{
			//Random roll
			float roll = Random_F();
			//Create the output vectors
			std::vector<std::vector<double>> out_1 = selected[j].TGM_fitness;
			std::vector<std::vector<double>> out_2 = selected[j + 1].TGM_fitness;

			//perform changes (remove the last position, and insert the new one at the beginning)
			out_1.insert(out_1.begin(), std::vector<double>(fcode.Objs(), 0));
			out_2.insert(out_2.begin(), std::vector<double>(fcode.Objs(), 0));
			out_1.pop_back();
			out_2.pop_back();

			//assign the vectors to the offsprings
			if (roll < 0.5f)
			{
				cross_indi[j].TGM_fitness = out_1;
				cross_indi[j + 1].TGM_fitness = out_2;
			}
			else
			{
				cross_indi[j].TGM_fitness = out_2;
				cross_indi[j + 1].TGM_fitness = out_1;
			}
			//check the correct sizes
			if (cross_indi[j].TGM_fitness.size() != TGM_size + 1 || cross_indi[j + 1].TGM_fitness.size() != TGM_size + 1)
				abort();
		}
	} //end of crossover loop

	  //add the last individual (in the case of odd size vector)
	if (psize_odd == true)
	{
		//if is, last individual is added to new vector without change
		cross_indi.push_back(selected.back());
		if (TGM == true)
		{
			cross_indi.back().TGM_fitness.insert(cross_indi.back().TGM_fitness.begin(), std::vector<double>(fcode.Objs(), 0));
			cross_indi.back().TGM_fitness.pop_back();

			if (cross_indi.back().TGM_fitness.size() != TGM_size + 1)
				abort();
		}

	}
	  //Checking if output vector have the right size
	if (cross_indi.size() != selected.size())
	{
		std::cout << "ERROR#10: CROSSOVER - VECTOR SIZE"; //ERROR#10: CROSSOVER - VECTOR SIZE
		system("pause");
		abort();
	}

	//end of function
	return cross_indi;
};

/****************************************
Crossover #B2 - Mult-Point Binary Crossover
****************************************/
template <typename tname>
class crossover_B2 : public crossover<tname>
{
public:
	/**Default constructor**/
	crossover_B2() : crossover(2, "Mult-Point Binary Crossover") {};
	~crossover_B2() {};
	/*
	*make crossover for the given population*
	@param selected vector containing selected individuals
	@param gapara GA_Parameters class object
	@param fcode function class object
	*/
	std::vector<tname> Crossover(const std::vector<tname> & selected, GA_parameters & gapara, function & fcode) const;
};

/*
*make crossover for the given population*
@param selected vector containing selected individuals
@param gapara GA_Parameters class object
@param fcode function class object
*/
template <typename tname>
std::vector<tname> crossover_B2<tname>::Crossover(const std::vector<tname> & selected, GA_parameters & gapara, function & fcode) const
{
	std::vector<tname> cross_indi;		//Output crossover vector
	int psize = selected.size();		//size of the given selection vector
	bool psize_odd = false;

	//Check if given selection vector is odd - crossover can be done only on the even vector
	if (psize % 2 == 1)
	{
		psize_odd = true;
		//decreasing size of selection vector to make it even
		psize -= 1;
	}
	//actual crossover loop
	for (int j = 0; j < psize; j += 2)
	{
		std::vector<bool> parent1_v = selected[j].Code_Bin_Show(); //parent 1 storage vector
		std::vector<bool> parent2_v = selected[j + 1].Code_Bin_Show(); //parent 2 storage vector
		float cross_r = Random_F();					//random 
		int bin_string_size = parent1_v.size();
													//Checking if crossover occur
		if (cross_r <= gapara.Cross_Prob())
		{
			std::vector<bool> code1;	//code of 1 parent storage vector
			std::vector<bool> code2;	//code of 2 parent storage vector

			//Get the vector of points for crossover
			std::vector<int> p_vector;
			for (int i = 0; i < MPCross_size; i++)
			{
				if (i == 0)
					p_vector.push_back(Random_I(1, (bin_string_size - 1) / MPCross_size));
				else if (i < MPCross_size - 1)
					p_vector.push_back(p_vector[i - 1] + Random_I(1, (bin_string_size - 1) / MPCross_size));
				else
					p_vector.push_back(Random_I(p_vector[i - 1], (bin_string_size - 1)));
			}
			p_vector.push_back(bin_string_size);
			//Loop for the binary string
			int i = 0;
			for (int k = 0; k < p_vector.size(); k++)
			{
				
				for (; i < p_vector[k]; i++)
				{
					if (k % 2 == 0)
					{
						//first offspring gets one gene from the first parent and second from the second one
						code1.push_back(parent1_v[i]);
						code2.push_back(parent2_v[i]);
					}
					else
					{
						//first offspring gets one gene from the second parent and second from the first one
						code1.push_back(parent2_v[i]);
						code2.push_back(parent1_v[i]);
					}
				}
			}
			//create the empty real values vector
			std::vector<double> code_temp = std::vector<double>(fcode.Vars(),0.);

			//creation of new individuals with calculated code
			tname c1{ code1,code_temp , fcode };
			tname c2{ code2,code_temp , fcode };

			//sending new individuals to output vector
			cross_indi.push_back(c1);
			cross_indi.push_back(c2);
		}
		//if not offsprings are genetically identical to parents
		else
		{
			cross_indi.push_back(selected[j]);
			cross_indi[j].Crowd_Dist_Set(0.0);
			cross_indi[j].Rank_Set(0);
			cross_indi.push_back(selected[j + 1]);
			cross_indi[j + 1].Crowd_Dist_Set(0.0);
			cross_indi[j + 1].Rank_Set(0);
		}
		if (TGM == true)
		{
			//Random roll
			float roll = Random_F();
			//Create the output vectors
			std::vector<std::vector<double>> out_1 = selected[j].TGM_fitness;
			std::vector<std::vector<double>> out_2 = selected[j + 1].TGM_fitness;

			//perform changes (remove the last position, and insert the new one at the beginning)
			out_1.insert(out_1.begin(), std::vector<double>(fcode.Objs(), 0));
			out_2.insert(out_2.begin(), std::vector<double>(fcode.Objs(), 0));
			out_1.pop_back();
			out_2.pop_back();

			//assign the vectors to the offsprings
			if (roll < 0.5f)
			{
				cross_indi[j].TGM_fitness = out_1;
				cross_indi[j + 1].TGM_fitness = out_2;
			}
			else
			{
				cross_indi[j].TGM_fitness = out_2;
				cross_indi[j + 1].TGM_fitness = out_1;
			}
			//check the correct sizes
			if (cross_indi[j].TGM_fitness.size() != TGM_size + 1 || cross_indi[j + 1].TGM_fitness.size() != TGM_size + 1)
				abort();
		}
	} //end of crossover loop

	  //add the last individual (in the case of odd size vector)
	if (psize_odd == true)
	{
		//if is, last individual is added to new vector without change
		cross_indi.push_back(selected.back());
		if (TGM == true)
		{

			cross_indi.back().TGM_fitness.insert(cross_indi.back().TGM_fitness.begin(), std::vector<double>(fcode.Objs(), 0));
			cross_indi.back().TGM_fitness.pop_back();

			if (cross_indi.back().TGM_fitness.size() != TGM_size + 1)
				abort();
		}
	}

	  //Checking if output vector have the right size
	if (cross_indi.size() != selected.size())
	{
		std::cout << "ERROR#10: CROSSOVER - VECTOR SIZE"; //ERROR#10: CROSSOVER - VECTOR SIZE
		system("pause");
		abort();
	}

	//end of function
	return cross_indi;
};
#endif // !CROSSOVER_H

/*template <typename tname>
class crossover
{
private:
	int ccode;
	std::string name;
public:
	crossover(int cc, std::string na) { ccode = cc; name = na; };
	crossover() { std::cout << "ERROR#08: CROSSOVER - CREATION"; abort(); }
	~crossover() {};
	std::string Name_Show() { return name; };
	virtual std::vector<tname>  Crossover(const tname & type1, const tname & type2, GA_parameters & gapara, function & fcode) const { std::vector<tname> a; return a; };				//Take 2 individuals and return 2 indi
};

/*Crossover 1
template <typename tname>
class crossover_1 : public crossover<tname>
{
public:
	crossover_1() : crossover(1, "Crossover_1") {};
	~crossover_1() {};
	std::vector<tname>  Crossover(const tname & type1, const tname & type2, GA_parameters & gapara, function & fcode) const;
};

template <typename tname>
std::vector<tname> crossover_1<tname>::Crossover(const tname & type1, const tname & type2, GA_parameters & gapara, function & fcode) const
{
	std::vector<tname> cross_indi;

	double cross_r = rand() % 101;   //random seed
	if (cross_r / 100 <= gapara.Cross_Prob())
	{
		std::vector<double> code1;
		std::vector<double> code2;
		for (int i = 0; i < fcode.Vars(); i++)
		{
			double parent1, parent2, y2, y1, beta, alpha, betaq;
			parent1 = type1.Code_Show(i);
			parent2 = type2.Code_Show(i);
			double upper_b = fcode.Bound(i, "upper");
			double lower_b = fcode.Bound(i, "lower");
			if (abs(parent1 - parent2) > 0.000001)
			{
				if (parent2 > parent1)			//checking which one is greater
				{
					y2 = parent2;
					y1 = parent1;
				}
				else
				{
					y2 = parent1;
					y1 = parent2;
				}
				if ((y1 - lower_b) > (upper_b - y2))			//find beta
					beta = 1 + 2 * ((upper_b - y2) / (y2 - y1));
				else
					beta = 1 + 2 * ((y1 - lower_b) / (y2 - y1));
				beta = 1 / beta;
				alpha = 2 - std::pow(beta, (gapara.Di_C() + 1));
				double rnd = rand();
				if (rnd / 100 < 1 / alpha)
				{
					alpha *= rnd / 100;
					betaq = std::pow(alpha, 1 / (gapara.Di_C() + 1));
				}
				else
				{
					alpha *= rnd / 100;
					alpha = 1 / (2 - alpha);
					betaq = std::pow(alpha, 1 / (gapara.Di_C() + 1));
				}
				code1.push_back(0.5*((y1 + y2) - betaq*(y2 - y1)));
				code2.push_back(0.5*((y1 + y2) + betaq*(y2 - y1)));
			}
			else
			{	/*In original was:
				//cross_indi[1] = parent2;
				//cross_indi[0] = parent1;
				code1.push_back(parent2);
				code2.push_back(parent1);
			}
		}
		tname c1{ code1, fcode };
		tname c2{ code2, fcode };
		cross_indi.push_back(c1);
		cross_indi.push_back(c2);
	}
	else
	{
		cross_indi.push_back(type1);
		cross_indi.push_back(type2);
	}
	if (cross_indi.size() != 2)
	{
		abort();
		std::cout << "ERROR#10: CROSSOVER - VECTOR SIZE"; //ERROR#10: CROSSOVER - VECTOR SIZE
	}
	return cross_indi;
}; */