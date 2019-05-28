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

#pragma once

//****************************************
//				CLASS header
//		INDIVIDUAL,POPULATION,
//		COLLECTIVE,PARETO FRONT
//****************************************

//****************************************
//		INDIVIDUAL - one individual with code and fitness
//		POPULATION - set of individuals, store also GA_Paramteres
//		COLLECTIVE - like population but contains also funtions for colelctive
//		PARETO FRONT - like population but contain PF generation
//****************************************

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




class individual
{
private:
	std::vector<double> code;				//variables
	std::vector<bool> code_bin;				//Variables in binary encoding
	std::vector<double> fitness;			//fitness storage
	std::vector<double> fitness_norm;		//normalised fitness storage
	bool cons_violation;					//constrains violation - true if constrains are violated
	std::vector<double> cons_val;			//Values of constrains

	/*For NSGA II / NSGA-III*/
	double crowd_dist;						//crowding distance	(degree for BCE) ;
	int rank;								//rank of the individual (Mark in BCE)
	

public:
	/*For MOAE/D*/
	std::vector<double> namda;				//weight vector
	std::vector<int>  table;				//neighbourhood table - in case of MTS: Enable[i] (in case of BCE - crowding neighbourhood); clone number for HEIA
	std::vector<double> saved_fitness;			//saved fitness - obj_norm for BCE
	std::vector<std::vector<double>> TGM_fitness; //Fitness used for transgenerational memory [0] - individual fitness [1-n+1] fitness of parents up to n generation
	double utility;				//utility value - in case of MTS: lcount[i]; for IBEA: indicator based fitness

private:
	/*
	*Generate code for the individual*
	@param fcode function accorting to which code will be created
	*/
	void Code_Set(function & fcode);
	/*
	*Constrains Calculation*
	@param fcode function for which constrain will be calculated
	*/
	void Cons_Viol_Calc(function & fcode);
public:
	/**Default constructor - throw ERROR, not used**/
	individual() { std::cout << "ERROR" << std::endl; };
	/*
	*Individual constructor - normal*
	@param fcode function for which individual will be created
	*/
	individual(function & fcode) { extern std::string MODE; Code_Set(fcode); if (MODE == "MOEAD" ||MODE == "BCE" ||  MODE == "MOEADMSF" || MODE == "MOEADPSF" || MODE == "MOEADM2M") utility = 1.; else if (MODE == "DMOEADD") utility =0.; cons_violation = false; };
	/*
	*Copying constructor - when code is given*
	@param Code code which have to be copied
	*/
	individual(std::vector<double> & Code, function & fcode) { code = Code;};
	/*
	*Copying constructor - when code is given*
	@param Code code which have to be copied
	*/
	individual(std::vector<bool> & Code, std::vector<double> & Code2, function & fcode) { code_bin = Code; code = Code2; };
	/**Default copying constructor**/
	individual(const individual & indi) { code = indi.code; code_bin = indi.code_bin; fitness = indi.fitness; fitness_norm = indi.fitness_norm; TGM_fitness = indi.TGM_fitness; cons_violation = indi.cons_violation; cons_val = indi.cons_val; rank = indi.rank; crowd_dist = indi.crowd_dist; namda = indi.namda; table = indi.table; saved_fitness = indi.saved_fitness; utility = indi.utility; };
	~individual() {};

	/**Returning fitness**/
	std::vector<double> Fitness_Show() const { if (MLSGA_norm_obj == false) return fitness; else return fitness_norm; };
	/*
	*Returning fitness for 1 objective *
	@param ix index of the objective
	*/
	double Fitness_Show(int ix) const { if (MLSGA_norm_obj == false) return fitness[ix]; else return fitness_norm[ix];};
	/**Returning code (genotype)**/
	std::vector<bool> Code_Bin_Show()const { return code_bin; };
	/**Returning code (genotype)**/
	std::vector<double> Code_Show()const { return code; };
	/*
	*Returning code (genotype) for 1 variable*
	@param ix index of the variable
	*/
	double Code_Show(int ix)const { return code[ix]; };
	/*
	*Set the code of the individual - for MTS*
	*/
	std::vector<double> & Code_Set() { return this->code; };
	/*
	*Set the code of the individual for 1 variable - for MTS*
	@param ix index of the variable
	*/
	double & Code_Set(int ix) { return this->code[ix]; };
	
	/* Decode the binary encoding to real values*/
	void Decode(function & fcode);

	/**Returning constrain violation parameter**/
	bool Cons_Viol_Show() const { return cons_violation; };
	/**Returning constrains values**/
	std::vector<double> Cons_Show() const { return cons_val; };
	/**Returning crowding distance**/
	double Crowd_Dist_Show() const { return crowd_dist; };
	/*
	*Set crowd distance to given value*
	@param dist - new distance
	*/
	void Crowd_Dist_Set(double dist) { this->crowd_dist = dist; };
	/**Returning rank**/
	int Rank_Show() const { return rank; };
	/*
	*Set rank to given value*
	@param nrank - new rank
	*/
	void Rank_Set(int nrank) { this->rank = nrank; };
	/*
	*Mutation of the individual*
	@param mcode mutation class type
	@param fcode function class type
	@param gapara GA parameters class type
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
	/**Rounding up fitness for PF calculation**/
	inline void Round_Up()
	{
		for (short i = 0; i < fitness.size(); i++) fitness[i] = ceil(fitness[i] * pow(10, PF_prec)) / pow(10, PF_prec);
	}	
	/**Shows which collective has lower weight - used for sorting of individuals due to weight vectors**/
	static bool Sort_NAMDA(const individual& c1, const individual& c2)
	{
		return c1.namda[0] < c2.namda[0];
	}
	/**Shows which collective has higher lcount - used for sorting in MTS**/
	static bool Sort_MTS(const individual& c1, const individual& c2)
	{
		return c1.utility > c2.utility;
	}
	/*
	*Individual creation - not used*
	@param fcode function for which individual will be created
	*/
	//void Create(function & fcode);
	/*
	*Fitness Calculation - MLS2*
	@param fcode function for which fitness will be calculated
	*/
	void Fitness_Calc(function & fcode);	
	/**Save to file**/
	void save();

	/*remove normalised values*/
	void Norm_Remove() { fitness_norm = fitness; };
};


class population
{
private:

	
	GA_parameters * ga_para;							//Addresses of GA parameters used
	mutation<short> * m_code;							//Addresses of mutation type used
	selection<individual> * s_code;						//Addresses of selection type used
	crossover<individual> * c_code;						//Addresses of crossover type used
protected:
	std::vector<individual> indiv;						//Individuals	
	function * f_code;									//Addresses of functions used
	int size;											//size of the population
public:

	std::vector<std::vector<short>> fit_index;	//indexes of fitnesses on each level used [0] - individuals, [1] - col, [2] -col of col etc.
protected:
	/**Selection of individuals for the crossover**/
	std::vector<individual> Selection() const;
	/**Selection of individuals for the crossover MLS3**/
	//std::vector<individual> Selection2(int cix) const;
public:
	/*
	*Normal constructor*
	@param fcode function class type used for individuals creation
	@param gapara GA parameters class type used for parameters
	@param mcode mutation class type used for mutation
	@param scode selection class type used for selection
	@param ccode crossover class type used for crossover
	@param amt amount of individuals which will be created by constructor
	*/
	population(function & fcode, GA_parameters & gapara, mutation<short> & mcode, selection<individual> & scode, crossover<individual> & ccode, int amt);
	/*
	*Copying constructor*
	@param indi vector of individuals which will be implemented in new population
	@param pop population from which parameters will be copyied
	*/
	population(std::vector<individual> &indi, population & pop);
	/*
	*Partly copying constructor*
	@param pop given population
	@param amt amount of individuals whom will be copied (from first one)
	*/
	population(const population & pop, int amt);		//Partly copying constructor
	/**Default copying constructor**/
	population(const population & pop);
	/**Default constructor - throw ERROR, not used**/
	population() { std::cout << "ERROR#04: POPULATION - CREATION"; system("pause"); abort(); }; //ERROR#04
	~population() {};

	/**Return the set of individuals**/
	std::vector<individual> Indiv_Show() const { return indiv; }
	/*
	*Return 1 individual*
	@param ix index of the individual
	*/
	individual Indiv_Show(int ix) const { return indiv[ix]; }
	/*Individual Set - allows to change public parameters*/
	std::vector<individual> & Indiv_Set() { return indiv; }
	individual & Indiv_Set(int ix) { return indiv[ix]; }
	/**Return size of the population**/
	inline int Size_Show() const { return size; }
	/**Individuals sorting from smallest to largest fitness - according to indiv fitness**/
	void Sort_Individuals();
	/**Individuals sorting from smallest to largest fitness - according to col fitness**/
	void Sort_Individuals2();
	/**Individuals sorting from smallest to largest fitness - according to col of col fitness**/
	//void Sort_Individuals3() { std::sort(indiv.begin(), indiv.end(), Sort3); }
	/*Individuals sorting from smalles to largest fitness - according to the first objective*/
	void Sort_Individuals_f1() { std::sort(indiv.begin(), indiv.end(), Sort_f1); }
	/**Individuals sorting from smallest to largest fitness - according to weight vector**/
	void Sort_Individuals_NAMDA() { std::sort(indiv.begin(), indiv.end(), individual::Sort_NAMDA); }
	/**Individuals sorting from largest to smallest lcount value - for MTS**/
	void Sort_Individuals_MTS() { std::sort(indiv.begin(), indiv.end(), individual::Sort_MTS); }

	/**Fitness calculation for all individuals**/
	virtual void Fitness_Calc();
	/*
	*Fitness calculation for 1 individual*
	@param indi address of the individual for which fitness will be calculated
	*/
	void Fitness_Calc(individual & indi);
	/*
	*Fitness calculation for 1 individual*
	@param indi index of the individual for which fitness will be calculated
	*/
	void Fitness_Calc(int indi);
	/*
	*Crossover of individuals*
	@param index index of the collective
	*/
	std::vector<individual> Crossover(short index = 0) const;
	/*
	*Crossover of individuals without implemented selection - for NSGAII*
	@param index index of the collective
	*/
	std::vector<individual> Crossover_NSGAII(std::vector<individual> selected, short index = 0) const;
	/**Mutation of the whole population**/
	void Mutation();
	/*
	*Mutation of 1 individual*
	@param indi address of the individual for which mutation will occur
	*/
	void Mutation(individual & indi);
	/*
	*Mutation of 1 individual*
	@param indi index of the individual for which mutation will occur
	*/
	void Mutation(int indi);
	/*
	**Mutation of the given population and evaluate fitness - for NSGAII**
	@param pop - given population
	*/
	void Mutation_NSGAII(std::vector<individual> & pop);
	/**Addition of the random individual to the population**/
	void Add();											
   /*
   *Addition of  the selected individual*
   @param indi individual which have to be added
   */
	void Add(individual & indi);
	/*
	*Addition of  the group of individuals*
	@param indi vector of individuals which have to be added
	*/
	void Add(std::vector<individual> & indi);
	/**Removal of the last individual from population**/
	void Remove();
	/*
	*Removal of the selected individual from population*
	@param index index of the individual to remove
	*/
	void Remove(int index);
	/*Showing the current function - for MOEAD*/
	function * FCode_Show() const { return f_code; }
	/*Showing the current GA_para*/
	GA_parameters * GAPara_Show() const { return ga_para; }
	/*Showing the current m_code*/
	mutation<short> * MCode_Show() const { return m_code; }
	/*Set the current m_code - for coovelution*/
	void MCode_Show(mutation<short> &mcode) { m_code = &mcode; }
	/*Showing the current c_code*/
	crossover<individual> * CCode_Show() const { return c_code; }
	/*Set the current c_code - for coovelution*/
	void CCode_Set(crossover<individual> & ccode) { c_code = &ccode; }
	/*Showing the current s_code*/
	selection<individual> * SCode_Show() const { return s_code; }

	/**Save the whole population to the file**/
	void save() { for (int i = 0; i < this->size; i++) indiv[i].save(); }
	/**Erase of the whole individual set from the population**/
	virtual void Erase();

	/**Shows which individual has lower fitness - used for sorting of individuals due to collectives fitness**/
	/*static bool Sort2(const individual& c1, const individual& c2)
	{
		if (c1.fitness.size() != 2)
			abort();
		if (MLStemp == 1 || MLStemp == 4)
			return (std::accumulate(c1.fitness.begin(), c1.fitness.end(), 0.0)) < (std::accumulate(c2.fitness.begin(), c2.fitness.end(), 0.0));
		else if (MLStemp == 2 || MLStemp == 3 || MLStemp == 5 || MLStemp == 6 || MLStemp == 7 || MLStemp == 8 || MLStemp == 9)
			return c1.fitness[fit_index_col_sel - 1] < c2.fitness[fit_index_col_sel - 1];
		else
			abort();
	};*/
	/**Shows which individual has lower fitness - used for sorting of individuals due to collectives of collectives fitness**/
	/*static bool Sort3(const individual& c1, const individual& c2)
	{
		if (MLStemp == 3)
			return c1.Code_Show(0) < c2.Code_Show(0);
		else if (MLStemp == 1 || MLStemp == 5)
			return c1.fitness[fit_index_sel - 1] < c2.fitness[fit_index_sel - 1];
		else if (MLStemp == 7 || MLStemp == 8 || MLStemp == 9)
			return (std::accumulate(c1.fitness.begin(), c1.fitness.end(), 0.0)) < (std::accumulate(c2.fitness.begin(), c2.fitness.end(), 0.0));
		else
			abort();
	}*/
	static bool Sort_f1(const individual& c1, const individual& c2)
	{
		return c1.Fitness_Show(0) < c2.Fitness_Show(0);
	}


};

class collective : public population
{
private:
	std::vector<individual> elite; //Elite individuals
	short index;			//Collective index
	double fitness;			//Collective fitness
	double min_fitness;		//min fitness in collective
	bool was_erased;		//If function was previousely erased
	std::string mode;		//mode used in collective
public:
	short index_vid;
	/*
	*Constructor for collectives creation*
	@param pop original population
	@param label labels of individuals
	@param ix index of the generatred collective
	@param 
	*/
	collective(const population & pop, std::vector<short> label, int ix, std::string m);
	/*
	*Copying constructor*
	@param indi set of individuals which will be implemented in new collective
	@param col collective from which parameters will be copied
	*/
	collective(std::vector<individual> &indi, collective & col) : population(indi, col) { index = col.index; elite = col.elite; was_erased = col.was_erased; mode = col.mode; index_vid = col.index_vid; };
	/**Default copying constructor**/
	collective(const collective & col);
	/**Default constructor - throw ERROR, not used**/
	collective() { std::cout << "ERROR#05: COLLECTIVE - CREATION"; system("pause"); abort(); }
	~collective() {};

	/**Return the fitness of the collective**/
	double Fitness_Show() { return fitness; };

	/**Return the index of the colelctive**/
	short Index_Show() { return index; };

	/*Change the index of the collective*/
	void Index_Set(short ix) { index = ix; };
	/**Shows which collective has lower fitness - used for sorting**/
	static bool Sort(const collective& c1, const collective& c2) { return c1.fitness < c2.fitness; };

	/**Putting the best individuals into elite**/
	void Elite_Create();
	/**replacement of the worst individuals with elite, also calculate the fitness of the collective**/
	void Elite_Replace();
	/**Fitness calculation of the collective**/
	void Fitness_Calc();
	/**Erase the whole individual set and whole elite set from the colelctive, but the index remains**/
	void Erase();
	/**Compare and return the min fitness**/
	double Min_Fitness_Show(double current_min_fit);
	/**Crossover**/
	std::vector<individual> Crossover() const;
	/**Return erased state**/
	bool Was_Erased() const { return was_erased; };
	/**Set erased state to the false**/
	void Clear_Erased() { was_erased = false; };
	/*Show the mode used in the current collective*/
	std::string Mode_Show() { return mode; };

	/**Collective elimination - not implemented**/
	void Elim();
	/**Collective sorting - not implemented**/
	void Sorting();
};

class pareto_front : public population
{
public:
	/*
	*Pareto front constructor*
	@param pop - population from which parameters will be copied
	*/
	pareto_front(const population & pop) : population(pop,0) {}
	/**Default constructor - throw ERROR, not used**/
	pareto_front() { std::cout << "ERROR#06: PARETO FRONT- CREATION"; system("pause"); abort(); }
	~pareto_front() {};

	/*
	*Pareto front search*
	@param pop - population for pareto front search
	*/
	void Pareto_Search(const population & pop);
	/*
	*Pareto front search*
	@param ind - ind for pareto front search
	*/
	int Pareto_Search(individual & ind);
	/*
	*Reduce the size of the Pareto Front to the given amount of individuals*
	@param size - output size of the pareto front (default PF_size - defined in const.h)
	*/
	void Pareto_Refine(int size = PF_size);
};

#endif //!CLASS_H