/**
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
along with this program.If not, see < https://www.gnu.org/licenses/>. 


				CLASS header
		Storage of evolutionary classes

		INDIVIDUAL - one individual with code and fitness
		POPULATION - set of individuals, store also GA_Paramteres
		COLLECTIVE - like population but contains also funtions for colelctive
		PARETO FRONT - like population but contain PF generation

*/

#pragma once


#ifndef CLASS_H
#define CLASS_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

#include "Fit_Functions.h"
#include "Mutation.h"
#include "Crossover.h"
#include "Selection.h"



/**Class for storing and operating on a single individual*/
class individual
{
private:
	std::vector<double> code;				///<Variables vector
	std::vector<bool> code_bin;				///<Variables vectorin binary encoding
	std::vector<double> fitness;			///<Fitness (objectives) storage
	std::vector<double> fitness_norm;		///<Normalised fitness (objectives) storage
	bool cons_violation;					///<Constrains violation indicator. True if constrains are violated
	std::vector<double> cons_val;			///<Current values of constrains.

	/*For NSGA II / NSGA-III*/
	double crowd_dist;						///<Crowding distance for NSGA-II	(degree for BCE);
	int rank;								///<Rank of the individual for NSGA-II (Mark in BCE)
	

public:
	/*For MOAE/D*/
	std::vector<double> namda;				///<Current weight vector
	std::vector<int>  table;				///<Current neighbourhood table - in case of MTS: Enable[i] (in case of BCE - crowding neighbourhood); clone number for HEIA
	std::vector<double> saved_fitness;			///<Saved fitness - obj_norm for BCE
	std::vector<std::vector<double>> TGM_fitness; ///<Fitness used for transgenerational memory [0] - individual fitness [1-n+1] fitness of parents up to n generation
	double utility;				///<Utility value - in case of MTS: lcount[i]; for IBEA: indicator based fitness

private:
	/**Generate code (variables) for the individual
	@param fcode - function according to which the code will be generated
	*/
	void Code_Set(function & fcode);
	/**Constrains Calculation
	@param fcode - function according to which the constraints will be calculated
	*/
	void Cons_Viol_Calc(function & fcode);
public:
	/**Default empty constructor. Throw ERROR, not used*/
	individual() { std::cout << "ERROR" << std::endl; };
	/**Individual constructor - normal
	@param fcode - function for which individual will be created
	*/
	individual(function & fcode) { extern std::string MODE; Code_Set(fcode); if (MODE == "MOEAD" ||MODE == "BCE" ||  MODE == "MOEADMSF" || MODE == "MOEADPSF" || MODE == "MOEADM2M") utility = 1.; else if (MODE == "DMOEADD") utility =0.; cons_violation = false; };
	/**Copying constructor, when the code (variables) of another individual is given
	@param Code - code which will be used to create new individual
	@param fcode - function will be used to create new individual
	*/
	individual(std::vector<double> & Code, function & fcode) { extern std::string MODE; Code_Set(fcode); if (MODE == "MOEAD" || MODE == "BCE" || MODE == "MOEADMSF" || MODE == "MOEADPSF" || MODE == "MOEADM2M") utility = 1.; else if (MODE == "DMOEADD") utility = 0.; cons_violation = false; code = Code;};
	/**Copying constructor, when the binary and real code (variables) of another individual is given
	@param Code - binary code which will be used to create new individual
	@param Code2 - real coded code which will be used to create new individual
	@param fcode - function will be used to create new individual
	*/
	individual(std::vector<bool> & Code, std::vector<double> & Code2, function & fcode) { code_bin = Code; code = Code2; };
	/**Default copying constructor. Resulting individual is a 1:1 copy.*/
	individual(const individual & indi) { code = indi.code; code_bin = indi.code_bin; fitness = indi.fitness; fitness_norm = indi.fitness_norm; TGM_fitness = indi.TGM_fitness; cons_violation = indi.cons_violation; cons_val = indi.cons_val; rank = indi.rank; crowd_dist = indi.crowd_dist; namda = indi.namda; table = indi.table; saved_fitness = indi.saved_fitness; utility = indi.utility; };
	~individual() {};

	/**Returning the fitness of individual*/
	std::vector<double> Fitness_Show() const { if (MLSGA_norm_obj == false) return fitness; else return fitness_norm; };
	/**Returning fitness for a single objective
	@param ix - index of the objective
	*/
	double Fitness_Show(int ix) const { if (MLSGA_norm_obj == false) return fitness[ix]; else return fitness_norm[ix];};

	/**Returning non normalised fitness (when normalisation is used). Used only for PF calculation*/
	std::vector<double> Fitness_Show(bool t) const { if (t) return fitness; else { abort();  return fitness_norm; } };

	/**Manually set the fitness. Used for penalty functions
	* @param obj - vector containing a new fitness.
	
	*/
	void Fitness_Set(std::vector<double>& obj) { if (MLSGA_norm_obj == false) this->fitness = obj; else this->fitness_norm = obj; }
	/**Returning binary code (variables)*/
	std::vector<bool> Code_Bin_Show()const { return code_bin; };
	/**Returning real coded code (variables)*/
	std::vector<double> Code_Show()const { return code; };
	/**Returning code (variables) for a single variable
	@param ix - index of the variable
	*/
	double Code_Show(int ix)const { return code[ix]; };
	/**Manually set the code of the individual. Used only for MTS algorithm
	*/
	std::vector<double> & Code_Set() { return this->code; };
	/**Set the code of the individual for a single variable. Used only for MTS algorithm
	@param ix - index of the variable
	*/
	double & Code_Set(int ix) { return this->code[ix]; };
	
	/** Decode the binary encoding to real values*/
	void Decode(function & fcode);

	/**Returning constrain violation parameter*/
	bool Cons_Viol_Show() const { return cons_violation; };
	/**Returning constrains values*/
	std::vector<double> Cons_Show() const { return cons_val; };
	/**Returning crowding distance*/
	double Crowd_Dist_Show() const { return crowd_dist; };
	/**Set crowd distance to given value.
	@param dist - new crowding distance
	*/
	void Crowd_Dist_Set(double dist) { this->crowd_dist = dist; };
	/**Returning the rank*/
	int Rank_Show() const { return rank; };
	/**Set rank to a given value
	@param nrank - new rank
	*/
	void Rank_Set(int nrank) { this->rank = nrank; };
	/**Mutation of the individual
	@param mcode - mutation class according to which the process will be performed
	@param fcode - function class according to which the process will be performed
	@param gapara - GA parameters according to which the process will be performed
	*/
	void Mutation(mutation<short> & mcode, function & fcode, GA_parameters & gapara)
	{
		if (ENCODING == "Real")
			mcode.Mutate(code, fcode, gapara);
		if (ENCODING == "Binary" || ENCODING == "Gray")
		{
			mcode.Mutate_Bin(code_bin, fcode, gapara);
			Decode(fcode);
		}
			
	}
	/**Rounding up the fitness values. For PF calculation*/
	inline void Round_Up()
	{
		for (short i = 0; i < fitness.size(); i++) fitness[i] = ceil(fitness[i] * pow(10, PF_prec)) / pow(10, PF_prec);
	}	
	/**Shows which individual has lower weight. Used for sorting of individuals due to weight vectors*/
	static bool Sort_NAMDA(const individual& c1, const individual& c2)
	{
		return c1.namda[0] < c2.namda[0];
	}
	/**Shows which collective has higher lcount. Used for sorting in MTS*/
	static bool Sort_MTS(const individual& c1, const individual& c2)
	{
		return c1.utility > c2.utility;
	}
	
	/**Fitness (objectives) calculation.
	@param fcode - function class according to which the process will be performed
	*/
	void Fitness_Calc(function & fcode);	
	/**Save results to the file*/
	void save();

	/**remove normalised values*/
	void Norm_Remove() { fitness_norm = fitness; };
};

/**Class for storing and operating on a whole population. It stores a group of individuals as a vector, and all parameters specific to the population*/
class population
{
private:

	
	GA_parameters * ga_para;							///<Addresses of GA parameters used
	mutation<short> * m_code;							///<Addresses of mutation type used
	selection<individual> * s_code;						///<Addresses of selection type used
	crossover<individual> * c_code;						///<Addresses of crossover type used
protected:
	std::vector<individual> indiv;						///<Individuals	
	function * f_code;									///<Addresses of functions used
	int size;											///<size of the population
public:

	std::vector<std::vector<short>> fit_index;	///<indexes of the fitnesses on each level of selection. [0] - individuals, [1] - col, [2] -col of col etc.
protected:
	/**Selection of individuals for the crossover. It is performed according to s_code.*/
	std::vector<individual> Selection() const;
public:
	/**Normal constructor.
	@param fcode - function classused for individuals creation
	@param gapara - GA parameters class cantaining all GA parameters
	@param mcode - mutation class used for mutation process
	@param scode - selection class used for selection process
	@param ccode - crossover class used for crossover procedure
	@param amt - amount of individuals which will be created by the constructor. Defines the size of population.
	*/
	population(function & fcode, GA_parameters & gapara, mutation<short> & mcode, selection<individual> & scode, crossover<individual> & ccode, int amt);
	/**Copying constructor, with new set of individuals.
	@param indi - vector of individuals which will be implemented in the new population
	@param pop - population from which all parameters will be copied
	*/
	population(std::vector<individual> &indi, population & pop);
	/**Partly copying constructor.
	@param pop - given population
	@param amt - amount of individuals that will be copied (from first one).
	*/
	population(const population & pop, int amt);		//Partly copying constructor
	/**Default copying constructor. Result is 1:1 copy.*/
	population(const population & pop);
	/**Default constructor. Throws ERROR, not used**/
	population() { std::cout << "ERROR#04: POPULATION - CREATION"; system("pause"); abort(); }; //ERROR#04
	~population() {};

	/**Return the vector of individuals*/
	std::vector<individual> Indiv_Show() const { return indiv; }
	/**Return a single individual
	@param ix - index of the individual
	*/
	individual Indiv_Show(int ix) const { return indiv[ix]; }
	/**Return the vector of individuals, but allows to change of public parameters*/
	std::vector<individual> & Indiv_Set() { return indiv; }
	/**Return a single individual, but allows to change of public parameters
	@param ix - index of the individual
	*/
	individual & Indiv_Set(int ix) { return indiv[ix]; }
	/**Return the size of the population*/
	inline int  Size_Show() const { return size; }
	/**Individuals sorting from smallest to largest fitness - according to indiv fitness*/
	void Sort_Individuals();
	/**Individuals sorting from smallest to largest fitness - according to col fitness*/
	void Sort_Individuals2();
	/*Individuals sorting from smallest to largest fitness - according to the first objective*/
	void Sort_Individuals_f1() { std::sort(indiv.begin(), indiv.end(), Sort_f1); }
	/**Individuals sorting from smallest to largest fitness - according to weight vector*/
	void Sort_Individuals_NAMDA() { std::sort(indiv.begin(), indiv.end(), individual::Sort_NAMDA); }
	/**Individuals sorting from largest to smallest lcount value - for MTS*/
	void Sort_Individuals_MTS() { std::sort(indiv.begin(), indiv.end(), individual::Sort_MTS); }

	/**Fitness calculation for whole population*/
	virtual void Fitness_Calc();
	/**Fitness calculation for a single external individual. It is used to calculate fitness of a individual according to parameters stored in population.Therefore it doesn't make any changes in this->
	@param indi - address of the individual for which fitness will be calculated
	*/
	void Fitness_Calc(individual & indi);
	/**Fitness calculation for a single individual.
	@param indi - index of the individual for which fitness will be calculated
	*/
	void Fitness_Calc(int indi);
	/**Crossover of the whole population. The process includes the parents selection procedure. Returns a vector of offspring individuals.
	@param index - index of the collective
	*/
	std::vector<individual> Crossover(short index = 0) const;
	/**Crossover of external pre-selected individuals, using population parameters. Used mostly for NSGAII
	* @param selected - vector of parents individuals
	@param index - index of the collective
	*/
	std::vector<individual> Crossover_NSGAII(std::vector<individual> selected, short index = 0) const;
	/**Mutation of the whole population*/
	void Mutation();
	/**Mutation of a single external individual.
	@param indi - address of the individual for which mutation will occur
	*/
	void Mutation(individual & indi);
	/**Mutation of a single individualwithin the population
	@param indi - index of the individual for which mutation will occur
	*/
	void Mutation(int indi);
	/**Mutation of the external population (as vector of individual). Used mostly for NSGAII. Includes routine for fitness calculation.
	@param pop - given population
	*/
	void Mutation_NSGAII(std::vector<individual> & pop);
	/**Adding of a random individual to the population*/
	void Add();											
   /**Adding a predefined individual to the population
   @param indi - individual which will be added
   */
	void Add(individual & indi);
	/**Adding a group of individuals to the population
	@param indi - vector of individuals which will be added
	*/
	void Add(std::vector<individual> & indi);
	/**Removing the last individual from population*/
	void Remove();
	/**
	*Removing a single individual from population*
	@param index - index of the individual to be removed
	*/
	void Remove(int index);
	/**Returning the address of current function - for MOEAD*/
	function * FCode_Show() const { return f_code; }
	/**Returning the address of GA_para*/
	GA_parameters * GAPara_Show() const { return ga_para; }
	/**Returning the address of the current mutation*/
	mutation<short> * MCode_Show() const { return m_code; }
	/** 
	Change the current mutation. Used for coovelution.
	@param mcode - address of new mutation class.	
	*/
	void MCode_Show(mutation<short> &mcode) { m_code = &mcode; }
	/*Returning the address of the current crossover*/
	crossover<individual> * CCode_Show() const { return c_code; }
	/**
	Change the current crossover. Used for coovelution.
	@param ccode - address of new crossover class.
	*/
	void CCode_Set(crossover<individual> & ccode) { c_code = &ccode; }
	/**Returning the address of the current selection*/
	selection<individual> * SCode_Show() const { return s_code; }

	/**Save the whole population to the file*/
	void save() { for (int i = 0; i < this->size; i++) indiv[i].save(); }
	/**Erase  the whole individual set from the population. It removes the individual vector while maintaining all other parameters. Used for collective reproduction*/
	virtual void Erase();

	/**Sort the population according to lowest values of fitness 1 (1st objective)*/
	static bool Sort_f1(const individual& c1, const individual& c2)
	{
		return c1.Fitness_Show(0) < c2.Fitness_Show(0);
	}


};
/**Class for storing and operating on a single collective (instead of population). It stores a group of individuals as a vector, and all parameters specific to the population and the collective*/
class collective : public population
{
private:
	std::vector<individual> elite; ///<Elite individuals
	short index;			///<Collective index
	double fitness;			///<Collective fitness
	double min_fitness;		///<min fitness in collective
	bool was_erased;		///<If function was previousely erased
	std::string mode;		///<GA mode used in the collective
public:
	short index_vid;
	/**Constructor for collectives creation.
	@param pop - original population, from which the collective will be created
	@param label - labels of individuals. Defines which individuals will be copied to collective from population.
	@param ix - index of the generatred collective
	@param m - GA mode. Defines which genetic algorithm will be used within collective.
	*/
	collective(const population & pop, std::vector<short> label, int ix, std::string m);
	/**Copying constructor for a new set of individuals, while maintaining collective parameters and population parameters from anothe collective.
	@param indi - set of individuals which will be implemented in new collective
	@param col - collective from which the parameters will be copied
	*/
	collective(std::vector<individual> &indi, collective & col) : population(indi, col) { index = col.index; elite = col.elite; was_erased = col.was_erased; mode = col.mode; index_vid = col.index_vid; };
	/**Default copying constructor. Result is 1:1 copy.**/
	collective(const collective & col);
	/**Default empty constructor. Throws ERROR, not used**/
	collective() { std::cout << "ERROR#05: COLLECTIVE - CREATION"; system("pause"); abort(); }
	~collective() {};

	/**Return the fitness of the collective. Fitness of the collective, not individuals inside.*/
	double Fitness_Show() { return fitness; };

	/**Return the index of the collective*/
	short Index_Show() { return index; };

	/**Change the index of the collective*/
	void Index_Set(short ix) { index = ix; };
	/**Shows which collective has lower fitness. Used for sorting of collectives*/
	static bool Sort(const collective& c1, const collective& c2) { return c1.fitness < c2.fitness; };

	/**Putting the best individuals into the elite storage*/
	void Elite_Create();
	/*Replacement of the worst individuals with  individuals taken from elite storage. It also calculate the fitness of the collective*/
	void Elite_Replace();
	/**Fitness calculation of the collective. It doesn't recalculate the fitness of individuals inside. Instead population::Fitness_Calc() should be used for that purpose.*/
	void Fitness_Calc();
	/**Erase the whole individual set and whole elite set from the collective, but other parameters are maintained.*/
	void Erase();
	/**Compare, replace (if necessary) and return the min fitness*/
	double Min_Fitness_Show(double current_min_fit);
	/**Crossover, according to the index of collective*/
	std::vector<individual> Crossover() const;
	/**Return erased state*/
	bool Was_Erased() const { return was_erased; };
	/**Set erased state to the false*/
	void Clear_Erased() { was_erased = false; };
	/**Show the GA mode used in the current collective*/
	std::string Mode_Show() { return mode; };
};
/**Class for storing and operating on a pareto optimal front. It stores a group of individuals as a vector, and all parameters specific to the population as well*/
class pareto_front : public population
{
public:
	/**Pareto front constructor. Only parameters are copied, and pareto front is empty when created.
	@param pop - population from which parameters will be copied
	*/
	pareto_front(const population & pop) : population(pop,0) {}
	/**Default empty constructor. Throws ERROR, not used**/
	pareto_front() { std::cout << "ERROR#06: PARETO FRONT- CREATION"; system("pause"); abort(); }
	~pareto_front() {};

	/**Pareto front search, by scanning the population for non-dominated individuals.
	@param pop - population for pareto front search
	*/
	void Pareto_Search(const population & pop);
	/**Pareto front search, by scanning a single population for a non-dominance.
	@param ind - ind for pareto front search
	*/
	int Pareto_Search(individual & ind);
	/**Reduce the size of the Pareto Front to the given amount of individuals
	@param size - output size of the pareto front. Default is defined in const.h)
	*/
	void Pareto_Refine(int size = PF_size);
};

#endif //!CLASS_H