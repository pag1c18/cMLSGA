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

#include <sstream>			//for copying values from file
#include <fstream>			//for copying values from file
#include <cmath>

#include "Class.h"
#include "Support_Functions.h"
#include "GA_Functions.h"
#include "imported.h"
#include "Sobol.h"
#include "TGM.h"
#include "NSGAII.h"
std::vector<short> sort_vect;
time_t selec_t;				//time of selection
time_t cross_t;				//time of crossover
time_t PF_t;				//time of pareto front creation
time_t save_t;				//time of results saving
//extern short MLSt;			//MLSt type currently used
extern int nfes;			//the number of function evaluations
extern int cons_viol_count;		//How many times constrains have been violated in the current run
extern int dyn_tau;												//time state - for dynamic function
extern double dyn_t;
extern std::string MODE;	

static bool Sort(const individual& c1, const individual& c2);
/*
*Generate code for the individual*
@param fcode function accorting to which code will be created
*/
void individual::Code_Set(function & fcode)
{
	if (ENCODING == "Real")
	{
		//set variables
		for (int i = 0; i < fcode.Vars(); i++)
		{

			//assign the random value to the variable
			double x = Random();
			code.push_back(fcode.Bound(i, "lower") + x * (fcode.Bound(i, "upper") - fcode.Bound(i, "lower")));


		}
	}
	else if (ENCODING == "Binary" || ENCODING == "Gray")
	{
		//calculate the total number of bytes
		int bit_num = fcode.Vars() * Binary_string_size;

		//Create the binary string
		for (int i = 0; i < bit_num; i++)
		{
			if (Random() < 0.5)
				code_bin.push_back(0);
			else
				code_bin.push_back(1);
		}

		//create the dummy real code
		code = std::vector<double>(fcode.Vars(), 0.);

		//Decode the string - for fitness calc purposes
		Decode(fcode);
	}
	if (TGM == true)
	{
		//create empty tgm vector
		std::vector<double> temp_vect(fcode.Objs(), 0.);
		TGM_fitness = std::vector<std::vector<double>>(TGM_size + 1, temp_vect);
	}

}

void individual::Decode(function & fcode)
{
	//Calculate the maximum value of binary string
	long long int max_bit_val = pow(2, Binary_string_size) - 1;
	std::vector<bool> code_bin_temp;
	if (ENCODING == "Binary")
		code_bin_temp = code_bin;
	else if (ENCODING == "Gray")	//To decode the Gray it need to be first converted to binary in temporary sting
	{
		for (int i = 0; i < code.size(); i++)
		{
			//push the MBR
			code_bin_temp.push_back(code_bin[i*Binary_string_size]);
			//calculate and push other bits
			for (int ibit = 1; ibit < Binary_string_size; ibit++)
			{
				if (code_bin_temp[ibit - 1 + i*Binary_string_size] == code_bin[ibit + i*Binary_string_size])
					code_bin_temp.push_back(0);
				else
					code_bin_temp.push_back(1);
			}
			//check the string size
		}
		if (code_bin_temp.size() != code_bin.size())
			abort();
	}
	//Decode each variable
	for (int i = 0; i < code.size(); i++)
	{
		long long int bit_val = 0;		//The value of binary string

		//Add the value of each bit
		for (int j = 0; j < Binary_string_size; j++)
		{
			bit_val += pow(2, Binary_string_size - j - 1) * code_bin_temp[i*Binary_string_size + j];
		}
		
		//Calculate the real value
		double x = (long double)bit_val / (long double)max_bit_val;
		code[i] = fcode.Bound(i, "lower") + x * (fcode.Bound(i, "upper") - fcode.Bound(i, "lower"));
	}
}

void individual::Cons_Viol_Calc(function & fcode)
{ 
	if (fcode.Time_Dep())
		cons_val = fcode.Cons_Calc(this->code, this->fitness, dyn_t);
	else
		cons_val = fcode.Cons_Calc(this->code, this->fitness); 
	cons_violation = false;
	for (int i = 0; i < fcode.Cons(); i++)
	{
		if (cons_val[i] < cons_check_param)
		{
			cons_violation = true;
			break;
		}
	}
}

/*
*Individual creation - not used*
@param fcode function for which individual will be created
*/
/*void individual::Create(function & fcode)
{ 
	//create code for the new individual
	this->Code_Set(fcode); 
	//calcualte fitness for the new individual
	this->Fitness_Calc(fcode);
};*/

/*
*Fitness Calculation - MLS2*
@param fcode function for which fitness will be calculated
*/
void individual::Fitness_Calc(function & fcode)
{ 
	
	//Calculate and assign fitness
	if (fcode.Time_Dep())
		this->fitness = fcode.Fitness_C(code, dyn_tau, dyn_t);
	else
		this->fitness = fcode.Fitness_C(code); 
	
	


	if (TGM == true)
	{
		//Calculate the fintess according to transgenerational memory
		TGM_Calc(this->TGM_fitness, this->fitness);
	}

	if ( MLSGA_norm_obj == true)
	{
		Idealpoint_Update(this->fitness);

		if (TGM == false)
			fitness_norm = Normalize_Objective(this->fitness);
		else
			fitness_norm = Normalize_Objective(this->TGM_fitness[0]);
	}
	else if (MODE == "MOEADPSF" || MODE == "MOEADMSF" || MODE == "MOEADM2M" || MODE == "MOEAD")
		Idealpoint_Update(this->fitness);

	//Calculate constrains
	if (fcode.Cons() > 0)
	{
		Cons_Viol_Calc(fcode);
	}
	else
		cons_violation = false;
	//Check if have to save results
	if (FITNESS_ALL == true && MODE == "Normal")
		//save results to file
		save();

}

/**Save to file**/
void individual::save()
{
	if (FITNESS_ALL != true)
		return;
	//Do not save individuals that are out of constrains
	if (this->cons_violation == true)
		return;

	time_t save_t_temp = clock();							//begin time
	extern std::ofstream graph, *graph_v;					//extern file outputs - MLSGA.cpp
	extern SimpleXlsx::CWorksheet sheet1;					//extern excel output - MLSGA.cpp
	
	//std::vector<SimpleXlsx::CellDataStr> data2;
	std::vector<SimpleXlsx::CellDataDbl> data;				//excel data vector
	//SimpleXlsx::CellDataStr cellStr; //For precise data
	SimpleXlsx::CellDataDbl cellDbl;						//excel cell
	
	//save fitness to the file
	for (short h = 0; h < fitness.size(); h++)	
	{
		double temp = fitness[h];		//fitness balue
		//std::ostringstream out;			//For precise data
		//out << std::fixed;
		//out.precision(prec);
		//out << temp;
		//cellStr.value = out.str();

		//assign value to the excel cell
		cellDbl.value = temp;
		//save to main graph file
		graph << temp << " ";

		//save to video output source
		if (VIDEO == true)
			graph_v[0] << temp << " ";
		//save to excel
		if (EXCEL_EXCEPTION != true)
			data.push_back(cellDbl);
	}
	//add row to the excel
	if (EXCEL_EXCEPTION != true)
		sheet1.AddRow(data);										
	//end the line of the main graph file
	graph << std::endl;

	//end the line of the video output source
	if (VIDEO == true)
		graph_v[0] << std::endl;

	//saving the save time
	save_t += clock() - save_t_temp;
}

/**Selection of individuals for the crossover**/
std::vector<individual> population::Selection() const
{
	//select individuals for the crossover
	std::vector<short> temp_fit_index = fit_index[0];
	return s_code->Select(indiv, indiv.size(), f_code[0].Cons(), temp_fit_index);
}

/**Selection of individuals for the crossover - MLS3**/
/*std::vector<individual> population::Selection2(int cix) const
{
	if (MLSt != 3 && MLSt != 6 && MLSt != 7)
		abort();
	//select individuals for the crossover
	return s_code->Select2(indiv, indiv.size(), f_code[0].Cons(), cix);
}*/

/*
*Normal constructor*
@param fcode function class type used for individuals creation
@param gapara GA parameters class type used for parameters
@param mcode mutation class type used for mutation
@param scode selection class type used for selection
@param ccode crossover class type used for crossover
@param amt amount of individuals which will be created by constructor
*/
population::population(function & fcode, GA_parameters & gapara, mutation<short> & mcode, selection<individual> & scode, crossover<individual> & ccode, int amt)
{
	if (FILE_INPUT == false && SOBOL == false)			//normal random generator
	{
		//create the individuals and add them to the storage vector
		for (int i = 0; i < amt; i++)
		{
			indiv.push_back(individual::individual(fcode));
		}
	}
	else if (SOBOL == false && FILE_INPUT == true)		//read values from the file
	{
		//open input file
		std::ifstream file;			//input file
		file.open(input_name);

		//loop for copying value for each individual
		for (int i = 0; i < amt; i++)
		{
			std::string element_temp;			//temporary string for one value
			std::string line;					//string containing whole line from the file
			
			//copy the line from the file
			std::getline(file, line);
			std::stringstream ele(line);		//source stream of line
			std::vector<double> matrix_temp;		//temporary matrix for storaging one line of values

			for (int j = 0; j < fcode.Vars(); j++)
			{
				//get one value from the line
				std::getline(ele, element_temp, ',');
				//copy the value from string to float
				double float_temp = std::stof(element_temp);
				//assign the real value (according to boundaries)
				float_temp = fcode.Bound(j, "lower") + float_temp*(fcode.Bound(j, "upper") - fcode.Bound(j, "lower"));
				//push to the temporary matrix
				matrix_temp.push_back(float_temp);
			}
			//push temporary matrix to the final matrix
			indiv.push_back(individual(matrix_temp, fcode));
		}
	}
	else if (SOBOL == true && FILE_INPUT == false)	//generate variables by SOBOL sequence
	{
		std::vector<std::vector<double>> data = Sobol_Sequence(gapara.Pop_Size(), fcode.Vars());
		//loop for copying value for each individual
		for (int i = 0; i < amt; i++)
		{
			//copy the line from the file
			std::vector<double> matrix_temp;		//temporary matrix for storaging one line of values

			for (int j = 0; j < fcode.Vars(); j++)
			{
				//get value from data set
				double float_temp = data[j][i];
				//assign the real value (according to boundaries)
				float_temp = fcode.Bound(j, "lower") + float_temp*(fcode.Bound(j, "upper") - fcode.Bound(j, "lower"));
				//push to the temporary matrix
				matrix_temp.push_back(float_temp);
			}
			//push temporary matrix to the final matrix
			indiv.push_back(individual(matrix_temp, fcode));
		}
	}
	else											//if no right path
		abort();

	if (indiv.size() != amt)
	{
		std::cout << "ERROR#13: POPULATION - VECTOR SIZE";			//ERROR#13: POPULATION - VECTOR SIZE
		system("pause");
		abort();
	}

	//assign class parameters
	size = amt;
	f_code = &fcode;
	ga_para = &gapara;
	m_code = &mcode;
	s_code = &scode;
	c_code = &ccode;
	fit_index = std::vector<std::vector<short>>{ {1,2} };
}
/*
*Copying constructor*
@param indi vector of individuals which will be implemented in new population
@param pop population from which parameters will be copyied
*/
population::population(std::vector<individual> &indi, population & pop)
{
	//copy the indiv vector from the source
	indiv = indi;

	//copy all class parameters from the source
	size = indi.size();
	f_code = pop.f_code;
	ga_para = pop.ga_para;
	m_code = pop.m_code;
	s_code = pop.s_code;
	c_code = pop.c_code;
	fit_index = pop.fit_index;
}
/*
*Partly copying constructor*
@param pop given population
@param amt amount of individuals whom will be copied (from first one)
*/
population::population(const population & pop, int amt)
{
	//copy class parameters from the source
	f_code = pop.f_code;
	ga_para = pop.ga_para;
	m_code = pop.m_code;
	s_code = pop.s_code;
	c_code = pop.c_code;
	fit_index = pop.fit_index;

	//copy the individuals from the source - only desired amount
	for (int i = 0; i < amt; i++)
		indiv.push_back(pop.indiv[i]);

	//check the size of the indiv vector
	if (indiv.size() != amt)
	{
		std::cout << "ERROR#13: POPULATION - VECTOR SIZE #2";			//ERROR#13: POPULATION - VECTOR SIZE
		system("pause");
		abort();
	}

	//set the size of the population
	size = amt;
}
/**Default copying constructor**/
population::population(const population & pop)
{
	//copy all class parameters from the source
	f_code = pop.f_code;
	ga_para = pop.ga_para;
	m_code = pop.m_code;
	s_code = pop.s_code;
	c_code = pop.c_code;
	indiv = pop.indiv;
	size = pop.size;
	fit_index = pop.fit_index;
	//Fitness_Calc();
	//Fitness_Stat_Calc();
}

/**Fitness calculation for all individuals**/
void population::Fitness_Calc()
{

	//calculate fitness for each individual
	for (int i = 0; i < indiv.size(); i++)
		indiv[i].Fitness_Calc(f_code[0]);

}

/*
*Fitness calculation for 1 individual*
@param indi address of the individual for which fitness will be calculated
*/
void population::Fitness_Calc(individual & indi)
{

	//calculate fitness for individual
	indi.Fitness_Calc(f_code[0]);
}

/*
*Fitness calculation for 1 individual*
@param indi index of the individual for which fitness will be calculated
*/
void population::Fitness_Calc(int indi)
{

	//calculate fitness for individual
	indiv[indi].Fitness_Calc(f_code[0]);
}

/**Crossover of individuals**/
std::vector<individual> population::Crossover(short index) const
{
	//check if the population is not empty
	if (indiv.size() <= 0)
		return indiv;


	time_t temp = clock();										//starting time of the selection
	std::vector<individual> selected;							//vector containing selected individuals
	//select the individuals for the crossover
	selected = Selection();
	
	
	//calculate the time of the selection
	selec_t += clock() - temp ;

	time_t temp2 = clock();										//starting time of the crossover
	
	//do crossover and create new vectors of individuals
	std::vector<individual> new_pop = c_code->Crossover(selected,ga_para[0],f_code[0]);			//new individuals set generation
	if (new_pop.size() != indiv.size())							//check if pop sizes matches
	{
		std::cout << "ERROR#10: POPULATION/CROSSOVER - VECTOR SIZE "; //ERROR#10: CROSSOVER - VECTOR SIZE
		system("pause");
		abort();
	}

	//calculate the time of the crossover
	cross_t += clock() - temp2 ;
	
	return new_pop;
}
/**Crossover of individuals without implemented selection - for NSGAII**/
std::vector<individual> population::Crossover_NSGAII(std::vector<individual> selected, short index) const
{
	//check if the population is not empty
	if (selected.size() <= 0)
		return indiv;
	//check if NSGAII selected
	if (MODE != "NSGAII" && MODE != "BCE" && MODE !="HEIA" && MODE != "MOEADM2M" && MODE != "IBEA")
	{
		std::cout << "**********************************\n"
			<< "Wrong MODE chosen.Cannot do Crossover_NSGAII for non NSGAII. Check Define.h\nProgram will terminate"
			<< "\n**********************************\n";
		system("pause");
		abort();
	}
	time_t temp2 = clock();										//starting time of the crossover

	//do crossover and create new vectors of individuals
	std::vector<individual> new_pop = c_code->Crossover(selected, ga_para[0], f_code[0]);			//new individuals set generation
	if (new_pop.size() != selected.size())							//check if pop sizes matches
	{
		std::cout << "ERROR#10: POPULATION/CROSSOVER - VECTOR SIZE "; //ERROR#10: CROSSOVER - VECTOR SIZE
		system("pause");
		abort();
	}

	//calculate the time of the crossover
	cross_t += clock() - temp2;

	return new_pop;
}
/**Mutation of the whole population**/
void population::Mutation()
{
	//Mutate the whole population
	for (int i = 0; i < indiv.size(); i++)
	{
		//Mutate the code of the individual
		indiv[i].Mutation(m_code[0], f_code[0], ga_para[0]);

	}
	//Fitness_Calc();
}

/*
*Mutation of 1 individual*
@param indi address of the individual for which mutation will occur
*/
void population::Mutation(individual & indi)
{
	//Mutate the code of the individual
	indi.Mutation(m_code[0], f_code[0], ga_para[0]);
	//Fitness_Calc(indi);
}

/*
*Mutation of 1 individual*
@param indi index of the individual for which mutation will occur
*/
void population::Mutation(int indi)
{
	//Mutate the code of the individual
	indiv[indi].Mutation(m_code[0], f_code[0], ga_para[0]);

	//Fitness_Calc(indi);
}

/**Mutate the given population and evaluate fitness - for NSGAII**/
void population::Mutation_NSGAII(std::vector<individual> & pop)
{
	//check if NSGAII selected
	if (MODE != "NSGAII" && MODE != "DMOEADD" & MODE != "MOEADM2M")
	{
		std::cout << "**********************************\n"
			<< "Wrong MODE chosen.Cannot do Crossover_NSGAII for non NSGAII or DMOEADD. Check Define.h\nProgram will terminate"
			<< "\n**********************************\n";
		system("pause");
		abort();
	}
	//Mutate the whole population
	for (int i = 0; i < pop.size(); i++)
	{
		//Mutate the code of the individual
		pop[i].Mutation(m_code[0], f_code[0], ga_para[0]);

		//Calculate fitness of the individual
		pop[i].Fitness_Calc(f_code[0]);
		nfes++;
	}
	//Fitness_Calc();
}

/**Add random individual to the population**/
void population::Add()
{
	//Push a random individual to the indiv vector
	indiv.push_back(individual::individual(f_code[0]));

	//Increase the size of the population
	size++;
}

/*
*Add individual*
@param indi individual which have to be added
*/
void population::Add(individual & indi)
{
	//Push the individual to the indiv vector
	indiv.push_back(indi);

	//Increase the size of the population
	size++;
}

/*
*Add group of individuals*
@param indi vector of individuals which have to be added
*/
void population::Add(std::vector<individual> & indi)
{
	//Add all of individuals to the indiv vector
	for (int i = 0; i < indi.size(); i++)
	{
		//Push the individual to the indiv vector
		indiv.push_back(indi[i]);
		//Increase the size of the population
		size++;
	}
}

/**Remove the last individual from population**/
void population::Remove()
{
	//Remove the last object from the indiv vector
	indiv.pop_back();

	//Decrease the size of the population
	size--;
}

/*
*Romove individual from population*
@param index index of the individual to remove
*/
void population::Remove(int index)
{
	//Remove the indexed object from the indiv vector
	indiv.erase(indiv.begin() + index);

	//Decrease the size of the population
	size--;
}

/**Erase the whole individual set from the population**/
void  population::Erase()
{
	//Clear the indiv vector
	indiv.clear();

	//Set the size of the population to 0
	size = 0;
}

/**Individuals sorting from smallest to largest fitness - according to indiv fitness**/
void population::Sort_Individuals() { sort_vect = fit_index[0]; std::sort(indiv.begin(), indiv.end(), Sort); }
/**Individuals sorting from smallest to largest fitness - according to col fitness**/
void population::Sort_Individuals2() { sort_vect = fit_index[1]; std::sort(indiv.begin(), indiv.end(), Sort); }

/*
*Constructor for collectives creation*
@param pop original population
@param label labels of individuals
@param ix index of the generatred collective
*/
collective::collective(const population & pop, std::vector<short> label, int ix, std::string m) : population(pop,0)
{
	//check if  index is greater or equal 0
	if (ix <= 0)
	{
		//if not throw error and aborting program
		std::cout << "ERROR#05: COLLECTIVE - CREATION";
		system("pause");
		abort();
	}
	index = ix;

	//assign individuals to collective according to labels
	for (int i = 0; i < pop.Size_Show(); i++)
	{
		
		//if a label of the individual matches the label of the collective add an individual to the collective
		if (label[i] == ix)
			Add(pop.Indiv_Show(i));
		else if (label[i] == 0)
			abort();
	}

	//Create elite population (empty)
	elite = std::vector<individual>();
	was_erased = true;
	mode = m;
}

/*
*Copying constructor*
@param indi set of individuals which will be implemented in new collective
@param col collective from which parameters will be copied
*/
collective::collective(const collective & col):population(col)
{
	//Copy the collective parameters from the source
	index = col.index;
	fitness = col.fitness;
	elite = col.elite;
	size = col.size;
	min_fitness = col.min_fitness;
	was_erased = col.was_erased;
	mode = col.mode;
	index_vid = col.index_vid;
}

/**Putting the best individuals into elite**/
void collective::Elite_Create()
{
	//Calculate the elite size
	int esize = indiv.size()*MLSGA_elite_size; //size of the elite

	/*if collective is very small, set elite size to 1*/
	if (esize == 0)
		esize = 1;

	/*Sort individuals, according to indiv fitness, and assign the best ones to the elite*/
	population::Sort_Individuals();

	elite.assign(indiv.begin(), indiv.begin() + esize);
}

/**replacement of the worst individuals with elite, also calculate the fitness of the collective**/
void collective::Elite_Replace()
{
	//Copy the elite size to the temp value
	int esize = elite.size(); //size of the elite

	//Check if elite population is empty
	if (esize == 0 || size == 0)
	{
		//throw error and abort the program
		std::cout << "ERROR#15: Elitism";
		system("pause");
		abort();
	}

	/*if (f_code[0].Time_Dep())
	{
		for (int i = 0; i < esize; i++)
			elite[i].Fitness_Calc(f_code[0]);
	}*/

	//Insert the elite to the indiv vector
	indiv.insert(indiv.begin(), elite.begin(), elite.end());

	//Sort individuals, according to indiv fitness, in order to find the worst ones
	population::Sort_Individuals();

	//Delete the worst individuals
	indiv.erase(indiv.end() - esize, indiv.end());
	
	//Check if sizes are the same
	if (indiv.size() != size)
	{
		std::cout << "ERROR#15: Elitism";
		system("pause");
		abort();
	}

	//Calculate the fitness of the collective
	this->Fitness_Calc();
}

/**Fitness calculation for the collective**/
void collective::Fitness_Calc()
{
	//get the size of the population and fitness
	int pop_size = this->indiv.size();
	//copy the population fitness  vector
	std::vector<std::vector<double>>indiv_fit;
	for (int i_ind = 0; i_ind < pop_size; i_ind++)
	{
		if (TGM == false)
		{
			indiv_fit.push_back(indiv[i_ind].Fitness_Show());
		}
		else
		{
			indiv_fit.push_back(indiv[i_ind].TGM_fitness[0]);
		}
	}
	
	//check the size
	if (indiv_fit.size() != this->indiv.size())
		abort();

	//calculate the fitness of the collective
	double col_fit_temp = 0;										//temporary collective fitness

	double min_fit;			//min fitness of the collective
	//loop for population size
	for (int i = 0; i < pop_size; i++)
	{
		short fit_size = fit_index[1].size();
		double fit_temp = 0;					//temporary fitness storage
												//calculate fitness of 1 individual
		for (short j = 0; j < fit_size; j++)
			fit_temp += indiv_fit[i][fit_index[1][j]-1] / (double)fit_size;
		//add fitnees of 1 individual to the total fitness
		col_fit_temp += fit_temp;
		//check if it is the lowest one
		if (ONE_OBJ_OVERRIDE == true)
		{
			if (i == 0)
				min_fit = fit_temp;
			else if (min_fit > fit_temp)
				min_fit = fit_temp;
		}
	}
	//assign minimum fitnes
	if (ONE_OBJ_OVERRIDE == true)
		min_fitness = min_fit;
	//calculate average fitness of individuals
	col_fit_temp /= pop_size;
	//assign average fitness as a fitness of the collective
	fitness = col_fit_temp;
}

/**Erase the whole individual set and whole elite set from the colelctive, but the index remains**/
void collective::Erase()
{
	//Erase the population parent class inhibited part
	population::Erase();

	//Set fitness to 0
	this->fitness = 0;

	//Clear the elite
	elite.clear();

	was_erased = true;
}

/**Compare and return the min fitness**/
double collective::Min_Fitness_Show(double current_min_fit)
{
	//check which fitness is lower
	if (current_min_fit > this->min_fitness)
		return this->min_fitness;
	else
		return current_min_fit;
}
/**Crossover fundtion for the collective**/
std::vector<individual> collective::Crossover() const
{
	return population::Crossover(index);
}
/*std::ostream & collective::show(std::ostream & os)
{
	os << "Collective " << index << "#" << std::endl;
	population::show(os);
	return os;
}*/


/*
*Pareto front search*
@param pop population for pareto front search
*/
void pareto_front::Pareto_Search(const population & pop)
{
	//Do pareto search only for many objective optimisation
	if (ONE_OBJ_OVERRIDE == true)
		abort();

	//check if given population is not empty
	if (pop.Size_Show() == 0)
		return;
	
	//save the begin time
	time_t temp31 = clock();
	

	std::vector<individual> temp_indi = pop.Indiv_Show();


	int temp_indi_size = temp_indi.size();

	//Clear the solutions that are out of constrains
	int ncons = pop.FCode_Show()[0].Cons();
	if (ncons > 0)
	{
		for (int i = 0; i < temp_indi_size; i++)
		{

			if (temp_indi[i].Cons_Viol_Show())
			{
				temp_indi.erase(temp_indi.begin() + i);
				i--;
				temp_indi_size--;
				cons_viol_count++;
				continue;
			}
			temp_indi[i].Round_Up();
			if (MLSGA_norm_obj == true)
			{

				temp_indi[i].Norm_Remove();

			}
		}
	}
	else if (MLSGA_norm_obj == true)
	{
		for (int i = 0; i < temp_indi_size; i++)
		{
			temp_indi[i].Norm_Remove();
		}
	}

	temp_indi_size = temp_indi.size();
	//Get the number of objectives
	int n_func = pop.FCode_Show()[0].Objs();


	std::vector<std::vector<double>> Fitness_PF;		//Temporary vector for the fitness storage
	
	std::vector<std::vector<double>> Fitness_temp;		//Temporary vector for the fitness storage


	for (int i = 0; i < this->size; i++)
		Fitness_PF.push_back(indiv[i].Fitness_Show());

	for (int i = 0; i < temp_indi_size; i++)
		Fitness_temp.push_back(temp_indi[i].Fitness_Show());

	//Remove dominated solutions
	for (int j = 0; j < temp_indi_size; j++)
	{

		//Find non-dominated points
				//loop for each point
		for (int i = j+1; i < temp_indi_size; i++)
		{

			short m_i = 0;			//temporary value for checking if i variable dominates
			short m_j = 0;			//temporary value for checking if j variable dominates

			//check, for every fitness, which point is dominated
			for (short k = 0; k < n_func; k++)
			{
				//check if i dominates j
				if (Fitness_temp[i][k] <= Fitness_temp[j][k])
					m_i++;
				//j dominates i
				else
					m_j++;
			}
			if (m_i == n_func)
			{
				temp_indi.erase(temp_indi.begin() + j);
				Fitness_temp.erase(Fitness_temp.begin() + j);
				j--;
				temp_indi_size--;
				break;
			}
			else if (m_j == n_func)
			{
				temp_indi.erase(temp_indi.begin() + i);
				Fitness_temp.erase(Fitness_temp.begin() + i);
				i--;
				temp_indi_size--;
			}
		}
	}




	//Find non-dominated points
	for (int j = 0; j < temp_indi_size; j++)
	{

		//Find non-dominated points
				//loop for each point
		for (int i = 0; i < this->size; i++)
		{

			short m_i = 0;			//temporary value for checking if i variable dominates
			short m_j = 0;			//temporary value for checking if j variable dominates

			//check, for every fitness, which point is dominated
			for (short k = 0; k < n_func; k++)
			{
				//check if i dominates j
				if (Fitness_PF[i][k] <= Fitness_temp[j][k])
					m_i++;
				//j dominates i
				else
					m_j++;
			}
			if (m_i == n_func)
			{
				temp_indi.erase(temp_indi.begin()+j);
				Fitness_temp.erase(Fitness_temp.begin() + j);
				j--;
				temp_indi_size--;
				break;
			}
			else if (m_j == n_func)
			{
				this->Remove(i);
				Fitness_PF.erase(Fitness_PF.begin() + i);
				i--;
			}
		}
	}

	this->Add(temp_indi);

	//Calculate the time
	PF_t += clock() - temp31;

	//Check if PF is not too big
	if (indiv.size() > (PF_refine_size+1000))
		//refine PF to the desired size
		this->Pareto_Refine(PF_refine_size);

}
/*
*Pareto front search*
@param ind - ind for pareto front search
*/
int pareto_front::Pareto_Search(individual & ind)
{
	if (ind.Cons_Viol_Show())
		return -1;
	//Do pareto search only for many objective optimisation
	if (ONE_OBJ_OVERRIDE == true)
		abort();

	//save the begin time
	time_t temp31 = clock();
	int n_obj = ind.Fitness_Show().size();	//number of objectives

	individual temp_ind = ind;

	//remove the normalised values
	if (MLSGA_norm_obj == true)
		temp_ind.Norm_Remove();

	//copy the values
	std::vector<double> ind_fit = temp_ind.Fitness_Show();
	std::vector<double> ind_cons = temp_ind.Cons_Show();
	for (int i = 0; i < indiv.size(); i++)
	{
		int val = 0;			//temporary value for result saving
		short m_i = 0;			//temporary value for checking if i variable dominates
		short n = 0;			//temporary value for checking if variables are the same
		short m_j = 0;			//temporary value for checking if j variable dominates

		//copy the fitness
		std::vector<double> PF_ind_fit = indiv[i].Fitness_Show();


		//copy the constraints
		std::vector<double> PF_ind_cons = indiv[i].Cons_Show();

		//check, for every fitness, which point is dominated
		for (int j = 0; j < n_obj; j++)
		{
			//check if i dominates given point
			if (PF_ind_fit[j] < ind_fit[j])
				m_i++;
			//given point dominates i
			else if (PF_ind_fit[j] > ind_fit[j])
				m_j++;
			
		}
		if (m_j && !m_i)
		{
			Remove(i);
			i--;
		}
		else if ((!m_j && m_i)||(!m_j && !m_i))
		{
			//Calculate the time
			PF_t += clock() - temp31;
			return -1;
		}
	}
	Add(temp_ind);
	//Calculate the time
	PF_t += clock() - temp31;

	//Check if PF is not too big
	if (indiv.size() > (MTS_PF_refine_size))
		//refine PF to the desired size
		Pareto_Refine(MTS_PF_refine_size - 200);
	return indiv.size()-1;
}
/*
*Reduce the size of the Pareto Front to the given amount of individuals*
@param size output size of the pareto front (default PF_size - defined in const.h)
*/
void pareto_front::Pareto_Refine(int size)
{
	
	//Do pareto search only for many objective optimisation
	if (ONE_OBJ_OVERRIDE == true)
		abort();

	//check if refining is needed
	if (this->size < size)
		return;
	time_t PF_t_temp = clock();
	
	short n_obj = indiv[0].Fitness_Show().size();
	//check number of objectives
	if (n_obj == 2)
	{
		this->Sort_Individuals_f1();

		//refine the whole pareto front to the given size
		while (this->size > size)
		{
			double temp_min_storage[2] = { 0,0 };					//for storage the current min value and it's index

			//find the point in the most crowded area
			for (int i = 1; i < (this->size - 1); i++)
			{

				//copy the individual fitness
				std::vector<double> ind_fit = indiv[i].Fitness_Show();
				std::vector<double> ind_fit_prev = indiv[i - 1].Fitness_Show();
				std::vector<double> ind_fit_next = indiv[i + 1].Fitness_Show();


				//calculate the distance of the point to it's neighbour points
				double temp_min = Distance(ind_fit, ind_fit_prev);
				double temp_min2 = Distance(ind_fit, ind_fit_next);

				if (temp_min == 0 || temp_min2 == 0)
				{
					temp_min_storage[1] = i;
					break;
				}
				
				//check which distance is greater
				if (temp_min2 > temp_min)
					temp_min = temp_min2;
				//check if this distance is the min one
				if (i == 1)
				{
					//store the value of lowest distance
					temp_min_storage[0] = temp_min;

					//store the index of the point with the lowst distance
					temp_min_storage[1] = i;
				}
				else if (temp_min < temp_min_storage[0])
				{
					//store the value of lowest distance
					temp_min_storage[0] = temp_min;

					//store the index of the point with the lowst distance
					temp_min_storage[1] = i;
				}
			}
			//remove the point in the most crowded area
			Remove(temp_min_storage[1]);
		}
	}
	else
	{
		//refine based on crowding distance
		//Calculate crowding disctance
		NSGAII::Crowding_Distance_Indices_Assign(this->indiv, 0, this->size - 1);
		short count = 0;
		while (this->size > size)
		{
			//recalcualte the crowding distance every 10 removals
			if (count == 10)
			{
				NSGAII::Crowding_Distance_Indices_Assign(this->indiv, 0, this->size - 1);
				count = 0;
			}
			double temp_min = this->indiv[0].Crowd_Dist_Show();
			int index = 0;
			//find and remove the most crowded points
			for (int i = 1; i < this->size; i++)
			{
				if (this->indiv[i].Crowd_Dist_Show() < temp_min)
				{
					temp_min = this->indiv[i].Crowd_Dist_Show();
					index = i;
				}
			}
			this->Remove(index);
			count++;
		}
	}
	PF_t += clock() - PF_t_temp;
	//check if the output has the right size
	if (indiv.size() != size || this->size != size)
		abort();
}

/*for (int i = f_code->Max_Min_Fit(0, "lower"); i <= f_code->Max_Min_Fit(0, "upper"); i += i_step)
{
	double temp_min_storage[2] = { 0,0 };					//for storage the current min value and it's index

															//fine the min value for the given parameters
	for (int j = 0; j < Pareto_Set.size(); j++)
	{
		double temp_min = Pareto_Set[j].Fitness_Show(0)*((double)(i) / (double)(size - 1.0)) + Pareto_Set[j].Fitness_Show(1)*((double)(size - 1.0 - i) / (double)(size - 1.0));
		if (j == 0)
		{
			temp_min_storage[0] = temp_min;
			temp_min_storage[1] = j;
		}
		else if (temp_min < temp_min_storage[0])
		{
			temp_min_storage[0] = temp_min;
			temp_min_storage[1] = j;
		}
	}
	Add(Pareto_Set[temp_min_storage[1]]);
}*/

/*
//create real PF for a given size
std::vector<std::vector<float>> real_PF = f_code->Plot_PF(-1, size);

//refine the whole pareto front to the given size
for (int i = 0; i < size; i++)
{
	double temp_min_storage[2] = { 0,0 };					//for storage the current min value and it's index

															//find the min value for the given parameters
	for (int j = 0; j < Pareto_Set.size(); j++)
	{
		//calculate the distance beetween real PF and given PF
		double temp_min = Distance2<double, float>(Pareto_Set[j].Fitness_Show(), real_PF[i]);

		//pick up the smallest one
		if (j == 0)
		{
			temp_min_storage[0] = temp_min;
			temp_min_storage[1] = j;
		}
		else if (temp_min < temp_min_storage[0])
		{
			temp_min_storage[0] = temp_min;
			temp_min_storage[1] = j;
		}
	}
	Add(Pareto_Set[temp_min_storage[1]]);
}*/

/*int rand = Random_I(0, size - 1);

	Add(Pareto_Set[rand]);
	//refine the whole pareto front to the given size
	for (int i = 1; i < size; i++)
	{
		
		double temp_max_storage[2] = { 0,0 };					//for storage the current min value and it's index
		
		//find the min value for the given parameters
		for (int j = 0; j < Pareto_Set.size(); j++)
		{
			double temp_max = 0;
			for (int k = 0; k < this->size; k++)
			{
				//calculate the distance beetween real PF and given PF
				temp_max += Distance<double>(Pareto_Set[j].Fitness_Show(), indiv[k].Fitness_Show())/(double)this->size;

				//pick up the largest one
				if (j == 0 && k == 0)
				{
					temp_max_storage[0] = temp_max;
					temp_max_storage[1] = j;
				}
				else if (temp_max > temp_max_storage[0])
				{
					temp_max_storage[0] = temp_max;
					temp_max_storage[1] = j;
				}
			}
		}
		Add(Pareto_Set[temp_max_storage[1]]);
	}*/

/**Shows which individual has lower fitness - used for sorting of individuals**/
static bool Sort(const individual& c1, const individual& c2)
{
	double c1_val = 0, c2_val = 0;
	short fit_size = sort_vect.size();

	for (short i = 0; i < fit_size; i++)
	{
		c1_val += c1.Fitness_Show(sort_vect[i] - 1);
		c2_val += c2.Fitness_Show(sort_vect[i] - 1);
	}
	return c1_val < c2_val;
};