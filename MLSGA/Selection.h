#pragma once
//****************************************
//				SELECTION header
//	Storage of different selection types
//			and basic selection class
//****************************************

#ifndef SELECTION_H
#define SELECTION_H
#include <iostream>
#include <vector>
#include "Define.h"
#include "Const.h"

//extern short MLSt;			//MLSt type currently used

template <typename tname>
class selection
{
private:
	short scode;					//Code of the selection type
	std::string name;				//Name of the selection type
public:
	/*
	*Default normal constructor*
	@param sc - code of the selection type
	@param na - name of the selection type
	*/
	selection(short sc, std::string na) { scode = sc; name = na; };
	/**Default constructor, throw ERROR - selection class cannot be empty**/
	selection() { std::cout << "ERROR#09: SELECTION - CREATION"; system("pause"); abort(); } //Default constuctor - selection type cannot be empty
	~selection() {};

	/**Returning name of the current selection type**/
	std::string Name_Show() { return name; };	//
	
	/*
	*Selection of the individuals and returning the vector*
	@param indiv - vector of individuals for the selection
	@param psize - size of the selected population
	@param ncons - number of constrains
	@param cix - index of crossover
	*/
	virtual std::vector<tname> Select(std::vector<tname> const indiv_pop, int psize, short ncons, std::vector<short> & fit_indexes) const { std::cout << "Error"; std::vector<tname> a; return a; }; 
	//virtual std::vector<tname> Select2(std::vector<tname> const indiv_pop, int psize, short ncons, short cix) const { std::cout << "Error"; std::vector<tname> a; return a; };
};

//****************************************
//			   SELECTION #1
//			  Roulette wheel
//For separate collective and individual fitness	
//****************************************
template <typename tname>
class roulette_wheel : public selection<tname>		//Roulette Wheel for 1st objective
{
public:
	/**Default constructor**/
	roulette_wheel() : selection(1, "Roulette_wheel_MLS7") {};
	~roulette_wheel() {};
	/*
	*Selection of the individuals and returning the vector*
	@param indiv - vector of individuals for the selection
	@param psize - size of the selected population
	@param ncons - number of constrains
	*/
	std::vector<tname> Select(std::vector<tname> const indiv_pop, int psize, short ncons, std::vector<short> & fit_indexes) const;
	//std::vector<tname>  Select2(std::vector<tname> const indiv_pop, int psize, short ncons, short cix) const;
};

/*
*Selection of the individuals and returning the vector*
@param indiv vector of individuals for the selection
@param psize size of the selected population
*/
template <typename tname>
std::vector<tname>  roulette_wheel<tname>::Select(std::vector<tname> const indiv_pop, int psize, short ncons, std::vector<short> & fit_indexes) const
{
	//copy the original population
	std::vector<tname> indiv = indiv_pop;

	//remove individuals that do not meet constrains
	if (ncons > 0)
	{
		for (int i = 0; i < indiv.size(); i++)
		{
			if (indiv[i].Cons_Viol_Show())
			{
				indiv.erase(indiv.begin() + i);
				i--;
			}
		}
	}
	//Copy size of the individual vector
	int indiv_size = indiv.size();

	//If all are violated just do normal selection
	if (indiv_size == 0)
	{
		indiv = indiv_pop;
		indiv_size = indiv.size();
	}


	std::vector<tname> index;							//vector of individuals - output
	double inverse_fitn_tot = 0;						//total inverse fitness
	std::vector<double> inverse_fitn(indiv_size, 0);;			//inverse fitness for each individual

	short fit_size = fit_indexes.size();
	
		
	for (int i = 0; i < indiv_size; i++)
	{
		double temp = 0;
		for (short j = 0; j < fit_size; j++)
		{
			if (TGM == false)
				temp += indiv[i].Fitness_Show(fit_indexes[j] - 1) / (double)fit_size;
			else
				temp += indiv[i].TGM_fitness[0][fit_indexes[j] - 1] / (double)fit_size;
		}

		inverse_fitn[i] = 1 / temp;
		inverse_fitn_tot += inverse_fitn[i];
	}
		
	

	inverse_fitn[0] /= inverse_fitn_tot;

	//calculation of the chance for each individual (size of the roulette wheel area)
	for (int i = 1; i < indiv_size; i++)
	{
		inverse_fitn[i] /= inverse_fitn_tot;
		inverse_fitn[i] += inverse_fitn[i - 1];				//chance calculation (area)
	}

	//chance of the last indivual have to be equal 1
	inverse_fitn[indiv_size - 1] = 1.0;

	//roulette wheel - selection process
	for (int i = 0; i < psize; i++)
	{
		int index_h;				//index of the selected individual
		float wheel_r = Random_F();	//rulette wheel random seed

									//loop for the rulette wheel
		for (int j = 0; j < indiv_size; j++)
		{
			//if the right one is found break the loop
			if (inverse_fitn[j] >= wheel_r)
			{
				index_h = j;
				break;
			}
		}

		//adding the selected individual to the output vector
		index.push_back(indiv[index_h]);
	}

	//checking if the output vector have the correct size
	if (index.size() != psize)
	{
		std::cout << "ERROR#11: SELECTION - VECTOR SIZE";				 //ERROR#11: SELECTION - VECTOR SIZE
		system("pause");
		abort();
	}
	return index;
}

/*End of Roulette Wheel MSL2*/

#endif // !SELECTION_H
