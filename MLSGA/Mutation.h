#pragma once
//****************************************
//				MUTATION header
//		Different mutation types storage 
//			and basic mutation class
//****************************************
#ifndef MUTATION_H
#define MUTATION_H
#include <vector>
#include <string>
#include "Random.h"
#include "GA_Functions.h"

template <typename tname>
class mutation
{
private:
	short mcode;							//Code of the mutation type
	std::string name;						//Name of the mutation type
public:	
	/*
	*Default normal constructor*
	@param mc code of the mutation type
	@param na name of the mutation type
	*/
	mutation(short mc, std::string nm) { mcode = mc; name = nm; };
	/**Default constructor, throw ERROR - mutation class cannot be empty**/
	mutation() { std::cout << "ERROR#07: MUTATION - CREATION"; system("pause"); abort(); }
	~mutation() {};

	/*
	*Mutation algorithm*
	@param code source code which will be mutated
	@param fcode function class type as boundaries source
	@param gapara GA_parameter class type
	*/
	virtual void Mutate(std::vector<double> & code, function & fcode, GA_parameters & gapara) const {};
	/*
	*Mutation algorithm*
	@param code source code which will be mutated
	@param fcode function class type as boundaries source
	@param gapara GA_parameter class type
	*/
	virtual void Mutate_Bin(std::vector<bool> & code, function & fcode, GA_parameters & gapara) const {};
	/**Return the name of the current mutation type**/
	std::string Name_Show() { return name; };
};

/****************************************
			Mutation #1 
	 real value distribution - suttipong
****************************************/
template <typename tname>
class mutation_1 : public mutation<tname>
{
public:
	/**Default constructor**/
	mutation_1() : mutation(1, "Polynomial - changed") {};
	~mutation_1() {};

	/*
	*Mutation algorithm*
	@param code source code which will be mutated
	@param fcode function class type as boundaries source
	@param gapara GA_parameter class type
	*/
	void Mutate(std::vector<double> & code, function & fcode, GA_parameters & gapara) const;
};

/*
*Mutation algorithm*
@param code source code which will be mutated
@param fcode function class type as boundaries source
@param gapara GA_parameter class type
*/
template <typename tname>
void mutation_1<tname>::Mutate(std::vector<double> & code, function & fcode, GA_parameters & gapara) const
{

	//loop for the variables/code
	for (int i = 0; i < code.size(); i++)
	{
		//Check if mutation occur
		float mut_rand = Random_F();
		if (mut_rand <= gapara.Mut_Prob())	//check if mutation occur
		{
			//assign values
			double code_val = code[i];					//variable storage
			double delta;								//mutation variable delta
			double upper_b = fcode.Bound(i, "upper");	//upper boundary
			double lower_b = fcode.Bound(i, "lower");	//lower boundary

			//check if the variable is not equal to boundaries
			if (code_val > lower_b && code_val < upper_b)					
			{
				//check if the value is closer to lower or upper boundary
				if ((code_val - lower_b) < (upper_b - code_val))
					delta = (code_val - lower_b) / (upper_b - lower_b);
				else
					delta = (upper_b - code_val) / (upper_b - lower_b);
				//delta cannot be too small
				if (delta < pow(10, -pow_eq_zero))
					delta = pow(10, -pow_eq_zero)/2;

				double indi = 1.0 / (gapara.Di_M() + 1.0);	//mutation variable indi - distribution index
				double rnd1 = Random();						//Random seed
				double deltaq;								//real change of the code value

				//Check if variable will have lower or higher value than original
				if (rnd1 <= 0.5)
				{
					double val;								//temporary variable

					//calculate the real change of the code value
					val = 2 * rnd1 + (1.0 - 2.0 * rnd1)*pow((1.0 - delta), (gapara.Di_M() + 1.));
					deltaq = pow(val, indi) - 1.0;
				}
				else
				{
					double val;								//temporary variable

					//calculating the real change of the code value
					val = 2 * (1.0 - rnd1) + 2.0 * (rnd1 - 0.5)*pow((1.0 - delta), (gapara.Di_M() + 1.));
					deltaq = 1.0 - pow(val, indi);
				}

				//Calculae a new code varaible value 
				code[i] += (deltaq*(upper_b - lower_b));

				//Check if is not out of boundary
				if (code[i] < lower_b)
					code[i] = lower_b;
				else if (code[i] > upper_b)
					code[i] = upper_b;
			}

			//if code[i] == lower_b || code[i] == upper_b;
			else							
			{
				//Assigning the random value
				double x = Random();
				code[i] = lower_b + x * (upper_b - lower_b);
				
			}
		}
	}//end of the loop for the variables/code
};

/*End of Mutation 1*/

/****************************************
Mutation #1b
real value distribution
****************************************/
template <typename tname>
class mutation_1b : public mutation<tname>
{
public:
	/**Default constructor**/
	mutation_1b() : mutation(1, "Real value distribution") {};
	~mutation_1b() {};

	/*
	*Mutation algorithm*
	@param code source code which will be mutated
	@param fcode function class type as boundaries source
	@param gapara GA_parameter class type
	*/
	void Mutate(std::vector<double> & code, function & fcode, GA_parameters & gapara) const;
};

/*
*Mutation algorithm*
@param code source code which will be mutated
@param fcode function class type as boundaries source
@param gapara GA_parameter class type
*/
template <typename tname>
void mutation_1b<tname>::Mutate(std::vector<double> & code, function & fcode, GA_parameters & gapara) const
{

	//loop for the variables/code
	for (int i = 0; i < code.size(); i++)
	{
		//Check if mutation occur
		float mut_rand = Random_F();
		if (mut_rand <= gapara.Mut_Prob())	//check if mutation occur
		{
			//assign values
			double code_val = code[i];					//variable storage
			double delta1, delta2;								//mutation variable delta
			double upper_b = fcode.Bound(i, "upper");	//upper boundary
			double lower_b = fcode.Bound(i, "lower");	//lower boundary

														//check if the variable is not equal to boundaries

			delta1 = (code_val - lower_b) / (upper_b - lower_b);
			delta2 = (upper_b - code_val) / (upper_b - lower_b);

			double indi = 1.0 / (gapara.Di_M() + 1.0);	//mutation variable indi - distribution index
			double rnd1 = Random();						//Random seed
			double deltaq;								//real change of the code value

														//Check if variable will have lower or higher value than original
			if (rnd1 <= 0.5)
			{
				double val;								//temporary variable
				//calculate the real change of the code value
				val = 2 * rnd1 + (1.0 - 2.0 * rnd1)*pow((1.0 - delta1), (gapara.Di_M() + 1.));
				deltaq = pow(val, indi) - 1.0;
			}
			else
			{
				double val;								//temporary variable
				//calculating the real change of the code value
				val = 2 * (1.0 - rnd1) + 2.0 * (rnd1 - 0.5)*pow((1.0 - delta2), (gapara.Di_M() + 1.));
				deltaq = 1.0 - pow(val, indi);
			}

			//Calculae a new code varaible value 
			code[i] += (deltaq*(upper_b - lower_b));

			//Check if is not out of boundary
			if (code[i] < lower_b)
				code[i] = lower_b;
			else if (code[i] > upper_b)
				code[i] = upper_b;
		}

		//if code[i] == lower_b || code[i] == upper_b;

	}//end of the loop for the variables/code
};

/*End of Mutation 1b*/

/****************************************
			Mutation #2
 real value random - Uniform mutation
****************************************/
template <typename tname>
class mutation_2 : public mutation<tname>
{
public:
	/**Default constructor**/
	mutation_2() : mutation(2, "Uniform") {};
	~mutation_2() {};

	/*
	*Mutation algorithm*
	@param code source code which will be mutated
	@param fcode function class type as boundaries source
	@param gapara GA_parameter class type
	*/
	void Mutate(std::vector<double> & code, function & fcode, GA_parameters & gapara) const;
};

/*
*Mutation algorithm*
@param code source code which will be mutated
@param fcode function class type as boundaries source
@param gapara GA_parameter class type
*/
template <typename tname>
void mutation_2<tname>::Mutate(std::vector<double> & code, function & fcode, GA_parameters & gapara) const
{
	//loop for the variables/code
	for (int i = 0; i < code.size(); i++)
	{
		//Check if mutation occur
		float mut_rand = Random_F();
		float mut_prob = 1.f/(float)fcode.Vars();
		if (mut_rand <= mut_prob)
		{
			//Get boundaries for the varaible
			double upper_b = fcode.Bound(i, "upper");
			double lower_b = fcode.Bound(i, "lower");
			
			//Change value of the varaible to the random value between boundaries
			double x = Random();
			code[i] = lower_b + x * (upper_b - lower_b);
		}
	} //end of the loop for the variables/code
};

/*End of Mutation 2*/

/****************************************
			Mutation #3
		real value distribution - PG - not working!!!!!!!!!!!!
****************************************/
template <typename tname>
class mutation_3 : public mutation<tname>
{
public:
	/**Default constructor**/
	mutation_3() : mutation(3, "Real value distribution - PG") {};
	~mutation_3() {};

	/*
	*Mutation algorithm*
	@param code source code which will be mutated
	@param fcode function class type as boundaries source
	@param gapara GA_parameter class type
	*/
	void Mutate(std::vector<double> & code, function & fcode, GA_parameters & gapara) const;
};

/*
*Mutation algorithm*
@param code source code which will be mutated
@param fcode function class type as boundaries source
@param gapara GA_parameter class type
*/
template <typename tname>
void mutation_3<tname>::Mutate(std::vector<double> & code, function & fcode, GA_parameters & gapara) const
{
	//loop for the variables/code
	for (int i = 0; i < code.size(); i++)
	{
		//Check if mutation occur
		float mut_rand = Random_F();
		if (mut_rand <= gapara.Mut_Prob())	//check if mutation occur
		{
			//assign values
			double code_val = code[i];					//variable storage
			
			//copy boundaries
			double upper_b = fcode.Bound(i, "upper");	//upper boundary
			double lower_b = fcode.Bound(i, "lower");	//lower boundary

			double delta = (upper_b - lower_b) * Di_m2; 								//mutation variable delta - change of variable
			if (code_val != upper_b && code_val != lower_b)
			{
				//calculate the real change
				double delta_r = Random(1.0, 10.0);					//real variable change
				delta_r = log(delta_r) * delta;

				//check if new value will be bigger or lower
				float rand = Random_F();
				if (rand > 0.5f)
				{
					code[i] += delta_r;
				}
				else
				{
					code[i] -= delta_r;
				}
				//Check if is not out of boundary
				if (code[i] < lower_b)
					code[i] = lower_b;
				else if (code[i] > upper_b)
					code[i] = upper_b;
			}

			//if code[i] == lower_b || code[i] == upper_b;
			else
			{
				//Assign a random value
				double x = Random();
				code[i] = lower_b + x * (upper_b - lower_b);
			}
		}
	}//end of the loop for the variables/code
}

//end of Mutation #3

/****************************************
Mutation #B1
Binary mutation 1#
****************************************/
template <typename tname>
class mutation_B1 : public mutation<tname>
{
public:
	/**Default constructor**/
	mutation_B1() : mutation(1, "Binary mutation") {};
	~mutation_B1() {};

	/*
	*Mutation algorithm*
	@param code source code which will be mutated
	@param fcode function class type as boundaries source
	@param gapara GA_parameter class type
	*/
	void Mutate_Bin(std::vector<bool> & code, function & fcode, GA_parameters & gapara) const;
};

/*
*Mutation algorithm*
@param code source code which will be mutated
@param fcode function class type as boundaries source
@param gapara GA_parameter class type
*/
template <typename tname>
void mutation_B1<tname>::Mutate_Bin(std::vector<bool> & code, function & fcode, GA_parameters & gapara) const
{
	//loop for the variables/code
	for (int i = 0; i < code.size(); i++)
	{
		//Check if mutation occur
		float mut_rand = Random();
		if (mut_rand <= gapara.Mut_Prob())	//check if mutation occur
		{
			if (code[i] == 0)
				code[i] = 1;
			else
				code[i] = 0;
		}
	}//end of the loop for the variables/code
}
//end of Mutation #4
#endif // !MUTATION_H
