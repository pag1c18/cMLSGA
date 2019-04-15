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
//		 FIT FUNTIOCNS header
//	Storage of fitness funtions
//			and basic funtionc class
//****************************************

#ifndef FIT_FUNCTION_H
#define FIT_FUNCTION_H
#include <string>
#include <iostream>
#include <vector>
#include <fstream>

#include "Struct.h"
#include "Const.h"


class function
{
private:
	std::string name_func;									//Name of  function
	int num_vars;											//Number of variables in function
	short num_objs;											//Number of objectives in function
	short index;											//Function intex
	short num_cons;											//Number of constaines
	std::vector<STRUCTURES::boundaries> bound;				// Boundaries for variables
	std::vector<STRUCTURES::boundaries> max_min_fit;		// Boundaries for max and min fitness

protected:
	bool time_dep = false;									//Dynamic function - 1 Yes, 0 No
	std::vector<float> t_vector;							//dynamic t variables for changing random functions

protected:
	/*
	*Set boundaries for variables*
	@param bn vector of the variables boundaries structure
	*/
	virtual void Bound_Set(std::vector<STRUCTURES::boundaries> & bn) { bound = bn; }
	/*
	*Set boundaries*
	@param mmf vector of the fitness boundaries structure
	*/
	void Max_Min_Fit_Set(std::vector<STRUCTURES::boundaries> & mmf) { max_min_fit = mmf; }	
	/*
	*Create the t_vector*
	@param n - size of the t_vector
	*/
	void T_Vector_Create(int n) { t_vector = std::vector<float>(n, 0.f); };
public:
	/**Default constructor, throw ERROR - function class cannot be empty**/
	function() {std::cout << "ERROR#02: FUNCTION - CREATION"; system("pause"); abort();};
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	function(const char * fname, int varsn, short objsn, short ix, short consn) 
	{
		name_func = fname; num_vars = varsn; num_objs = objsn; index = ix; num_cons = consn;
	};
	~function() {};

	/*
	*Pareto Front plotting and saving to file*
	@param indexr index of the run
	*/
	virtual std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);																//Plotting real PF and saving it to file
	/*
	*Pareto Front plotting and saving to file - dynamic*
	@param indexr index of the run
	@param t current time
	*/																																																		
	virtual std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size) { std::vector<std::vector<double>> temp; return temp; };
	/*
	*Returning boundaries for the variable*
	@param i index of the variable
	@param c which boundary - "lower"/"upper"
	*/
	double Bound(int i, const char * c)														
	{ 
		if (i >= num_vars || i < 0 || (c != "lower" && c != "upper"))
		{
			std::cout << "ERROR#03: BOUNDARIES - SHOW";
			system("pause");
			abort();
		}
		if (c == "lower")
			return bound[i].lower;
		else
			return bound[i].upper;
	}
	/**Returning boundaries for the variables as a vector of structures**/
	std::vector<STRUCTURES::boundaries> Bound() { return bound; };
	/*
	*Returning boundaries for the fitness*
	@param i index of the variable
	@param c which boundary - "lower"/"upper"
	*/
	double Max_Min_Fit(int i, const char * c)												
	{
		if (i >= num_vars || i < 0 || (c != "lower" && c != "upper"))
		{
			std::cout << "ERROR#03: BOUNDARIES - SHOW";
			system("pause");
			abort();
		}
		if (c == "lower")
			return max_min_fit[i].lower;
		else
			return max_min_fit[i].upper;
	}
	/**Returning boundaries for the fitness as a vector of structures**/
	std::vector<STRUCTURES::boundaries> Max_Min_Fit() { return max_min_fit; };				
	/**Retrurning the index of the current fucntion**/
	int Index() const { return index; } ;													
	/**Returning the number of variables for the current function**/
	int Vars() const { return num_vars; };												
	/**Returning the number of objectives for the current function**/
	short Objs() const { return num_objs; };	
	/**Returning the number of constrains for the current function**/
	short Cons() const { return num_cons; };
	/**Returning time dependency**/
	bool Time_Dep() const { return time_dep; };
	/**Returning the name of the current function**/
	std::string & Name_Show() { return name_func; };	
	/*
	*Fitness calculation for the given code*
	@param code - vector with the code
	*/
	virtual std::vector<double> Fitness_C(const std::vector<double> &code) { std::cout << "ERROR#: FUNCTION - FITNESS_C"; system("pause"); abort(); };
	/*
	*Fitness calculation for the given code - dynamic*
	@param code - vector with the code
	@param t -  current time
	*/
	virtual std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t) { std::cout << "ERROR#: FUNCTION - FITNESS_C"; system("pause"); abort(); };
	
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated**
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	virtual std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0) { abort(); }

	/*
	*Min fitness calculation for one objective optimisation*
	@param pareto_front real pareto front for which fitness will be calculated
	*/
	double Min_Fitness_Get(const std::vector<std::vector<double>> & pareto_front);
	/*
	*Update the t_vector*
	@param s - if it have to be updated by one step (true) or zeroed (false)
	*/
	void T_Vector_Update(bool s = true);
};

class dynamic_function : public function
{
public:
	/**Default constructor, throw ERROR - function class cannot be empty**/
	dynamic_function() { std::cout << "ERROR#02: FUNCTION - CREATION"; system("pause"); abort(); };
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	dynamic_function(const char * fname, int varsn, short objsn, short ix, short consn) : function(fname, varsn, objsn, ix, consn)
	{
		time_dep = true;
	};
	~dynamic_function() {};

};


/********************************************
			UNCONSTRAINED FUNCTIONS
*********************************************/

//****************************************
//			Function #1
//		Zitzler–Deb–Thiele's function 1
//****************************************
class ZDT1 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	ZDT1(const char * fname = "ZDT1", int varsn = 30, short objsn = 2, short ix = 1, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~ZDT1() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//			Function #2
//		Zitzler–Deb–Thiele's function 2
//****************************************
class ZDT2 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	ZDT2(const char * fname = "ZDT2", int varsn = 30, short objsn = 2, short ix = 2, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~ZDT2() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//			Function #3 
//		Zitzler–Deb–Thiele's function 3
//****************************************
class ZDT3 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	ZDT3(const char * fname = "ZDT3", int varsn = 30, short objsn = 2, short ix = 3, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~ZDT3() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//			Function #4
//		Zitzler–Deb–Thiele's function 4
//****************************************
class ZDT4 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	ZDT4(const char * fname = "ZDT4", int varsn = 10, short objsn = 2, short ix = 4, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~ZDT4() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//			Function #5
//		Zitzler–Deb–Thiele's function 5
//****************************************
class ZDT5 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	ZDT5(const char * fname = "ZDT5", int varsn = 11, short objsn = 2, short ix = 5, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~ZDT5() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//			Function #6
//		Zitzler–Deb–Thiele's function 6
//****************************************
class ZDT6 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	ZDT6(const char * fname = "ZDT6", int varsn = 10, short objsn = 2, short ix = 6, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~ZDT6() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};



//****************************************
//			Function #7
//		Deb–Thiele-Laumanns-Zitzler–'s function 1
//****************************************
class DTLZ1 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	DTLZ1(const char * fname = "DTLZ1", int varsn = n_func_obj + 4, short objsn = n_func_obj, short ix = 7, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~DTLZ1() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//			Function #8
//		Deb–Thiele-Laumanns-Zitzler–'s function 2
//****************************************
class DTLZ2 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	DTLZ2(const char * fname = "DTLZ2", int varsn = n_func_obj + 9, short objsn = n_func_obj, short ix = 8, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~DTLZ2() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//			Function #9
//		Deb–Thiele-Laumanns-Zitzler–'s function 3
//****************************************
class DTLZ3 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	DTLZ3(const char * fname = "DTLZ3", int varsn = n_func_obj + 9, short objsn = n_func_obj, short ix = 9, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~DTLZ3() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//			Function #10
//		Deb–Thiele-Laumanns-Zitzler–'s function 4
//****************************************
class DTLZ4 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	DTLZ4(const char * fname = "DTLZ4", int varsn = n_func_obj + 9, short objsn = n_func_obj, short ix = 10, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~DTLZ4() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//			Function #11
//		Deb–Thiele-Laumanns-Zitzler–'s function 5
//****************************************
class DTLZ5 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	DTLZ5(const char * fname = "DTLZ5", int varsn = n_func_obj + 9, short objsn = n_func_obj, short ix = 11, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~DTLZ5() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//			Function #12
//		Deb–Thiele-Laumanns-Zitzler–'s function 6
//****************************************
class DTLZ6 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	DTLZ6(const char * fname = "DTLZ6", int varsn = n_func_obj + 9, short objsn = n_func_obj, short ix = 12, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~DTLZ6() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//			Function #13
//		Deb–Thiele-Laumanns-Zitzler–'s function 7
//****************************************
class DTLZ7 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	DTLZ7(const char * fname = "DTLZ7", int varsn = n_func_obj + 19, short objsn = n_func_obj, short ix = 13, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~DTLZ7() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//			Function #14
//				MOP1
//****************************************
class MOP1 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	MOP1(const char * fname = "MOP1", int varsn = 10, short objsn = 2, short ix = 14, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~MOP1() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//			Function #15
//				MOP2
//****************************************
class MOP2 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	MOP2(const char * fname = "MOP2", int varsn = 10, short objsn = 2, short ix = 15, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~MOP2() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//			Function #16
//				MOP3
//****************************************
class MOP3 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	MOP3(const char * fname = "MOP3", int varsn = 10, short objsn = 2, short ix = 16, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~MOP3() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//			Function #17
//				MOP4
//****************************************
class MOP4 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	MOP4(const char * fname = "MOP4", int varsn = 10, short objsn = 2, short ix = 17, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~MOP4() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//			Function #18
//				MOP5
//****************************************
class MOP5 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	MOP5(const char * fname = "MOP5", int varsn = 10, short objsn = 2, short ix = 18, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~MOP5() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//			Function #19
//				MOP6
//****************************************
class MOP6 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	MOP6(const char * fname = "MOP6", int varsn = 10, short objsn = 3, short ix = 19, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~MOP6() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//			Function #20
//				MOP7
//****************************************
class MOP7 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	MOP7(const char * fname = "MOP7", int varsn = 10, short objsn = 3, short ix = 20, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~MOP7() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//				Function #21
//					UF1
//****************************************
class UF1 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UF1(const char * fname = "UF1", int varsn = 30, short objsn = 2, short ix = 21, short consn = 0)
		: function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~UF1() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//				Function #22
//					UF2
//****************************************
class UF2 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UF2(const char * fname = "UF2", int varsn = 30, short objsn = 2, short ix = 22, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~UF2() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C( const std::vector<double> &code);
};

//****************************************
//				Function #23
//					UF3
//****************************************
class UF3 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UF3(const char * fname = "UF3", int varsn = 30, short objsn = 2, short ix = 23, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~UF3() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//				Function #24
//					UF4
//****************************************
class UF4 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UF4(const char * fname = "UF4", int varsn = 30, short objsn = 2, short ix = 24, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~UF4() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//				Function #25
//					UF5
//****************************************
class UF5 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UF5(const char * fname = "UF5", int varsn = 30, short objsn = 2, short ix = 25, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~UF5() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//				Function #26
//					UF6
//****************************************
class UF6 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UF6(const char * fname = "UF6", int varsn = 30, short objsn = 2, short ix = 26, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~UF6() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//				Function #27
//					UF7
//****************************************
class UF7 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UF7(const char * fname = "UF7", int varsn = 30, short objsn = 2, short ix = 27, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~UF7() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//				Function #28
//					UF8
//****************************************
class UF8 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UF8(const char * fname = "UF8", int varsn = 30, short objsn = 3, short ix = 28, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~UF8() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//				Function #29
//					UF9
//****************************************
class UF9 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UF9(const char * fname = "UF9", int varsn = 30, short objsn = 3, short ix = 29, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~UF9() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//				Function #30
//					UF10
//****************************************
class UF10 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UF10(const char * fname = "UF10", int varsn = 30, short objsn = 3, short ix = 30, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~UF10() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//				Function #31
//					WFG1
//****************************************
class WFG1 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	WFG1(const char * fname = "WFG1", int varsn = (2 * (n_func_obj - 1) + 20), short objsn = n_func_obj, short ix = 31, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~WFG1() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//				Function #32
//					WFG2
//****************************************
class WFG2 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	WFG2(const char * fname = "WFG2", int varsn = (2 * (n_func_obj - 1) + 20), short objsn = n_func_obj, short ix = 32, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~WFG2() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//				Function #33
//					WFG3
//****************************************
class WFG3 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	WFG3(const char * fname = "WFG3", int varsn = (2 * (n_func_obj - 1) + 20), short objsn = n_func_obj, short ix = 33, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~WFG3() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//				Function #34
//					WFG4
//****************************************
class WFG4 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	WFG4(const char * fname = "WFG4", int varsn = (2 * (n_func_obj - 1) + 20), short objsn = n_func_obj, short ix = 34, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~WFG4() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//				Function #35
//					WFG5
//****************************************
class WFG5 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	WFG5(const char * fname = "WFG5", int varsn = (2 * (n_func_obj - 1) + 20), short objsn = n_func_obj, short ix = 35, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~WFG5() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//				Function #36
//					WFG6
//****************************************
class WFG6 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	WFG6(const char * fname = "WFG6", int varsn = (2 * (n_func_obj - 1) + 20), short objsn = n_func_obj, short ix = 36, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~WFG6() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//				Function #37
//					WFG7
//****************************************
class WFG7 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	WFG7(const char * fname = "WFG7", int varsn = (2 * (n_func_obj - 1) + 20), short objsn = n_func_obj, short ix = 37, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~WFG7() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//				Function #38
//					WFG8
//****************************************
class WFG8 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	WFG8(const char * fname = "WFG8", int varsn = (2 * (n_func_obj - 1) + 20), short objsn = n_func_obj, short ix = 38, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~WFG8() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//				Function #39
//					WFG9
//****************************************
class WFG9 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	WFG9(const char * fname = "WFG9", int varsn = (2 * (n_func_obj - 1) + 20), short objsn = n_func_obj, short ix = 39, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~WFG9() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//				Function #40
//					IMB1
//****************************************
class IMB1 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	IMB1(const char * fname = "IMB1", int varsn = 10, short objsn = 2, short ix = 40, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~IMB1() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//				Function #41
//					IMB2
//****************************************
class IMB2 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	IMB2(const char * fname = "IMB2", int varsn = 10, short objsn = 2, short ix = 41, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~IMB2() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//				Function #42
//					IMB3
//****************************************
class IMB3 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	IMB3(const char * fname = "IMB3", int varsn = 10, short objsn = 2, short ix = 42, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~IMB3() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//				Function #43
//					IMB4
//****************************************
class IMB4 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	IMB4(const char * fname = "IMB4", int varsn = 10, short objsn = 3, short ix = 43, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~IMB4() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//				Function #44
//					IMB5
//****************************************
class IMB5 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	IMB5(const char * fname = "IMB5", int varsn = 10, short objsn = 3, short ix = 44, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~IMB5() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//				Function #45
//					IMB6
//****************************************
class IMB6 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	IMB6(const char * fname = "IMB6", int varsn = 10, short objsn = 3, short ix = 45, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~IMB6() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//				Function #46
//					IMB7
//****************************************
class IMB7 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	IMB7(const char * fname = "IMB7", int varsn = 10, short objsn = 2, short ix = 46, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~IMB7() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//				Function #47
//					IMB8
//****************************************
class IMB8 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	IMB8(const char * fname = "IMB8", int varsn = 10, short objsn = 2, short ix = 47, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~IMB8() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//				Function #48
//					IMB9
//****************************************
class IMB9 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	IMB9(const char * fname = "IMB9", int varsn = 10, short objsn = 2, short ix = 48, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~IMB9() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//				Function #49
//					IMB10
//****************************************
class IMB10 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	IMB10(const char * fname = "IMB10", int varsn = 10, short objsn = 3, short ix = 49, short consn = 0)
		: function(fname, varsn, objsn, ix, consn)
	{
		Bound_Set();
	};
	~IMB10() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};



/********************************************
			CONSTRAINED FUNCTIONS
*********************************************/

//****************************************
//				Function #50
//					CF1
//****************************************
class CF1 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	CF1(const char * fname = "CF1", int varsn = 10, short objsn = 2, short ix = 50, short consn = 1)
		: function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~CF1() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};
//****************************************
//				Function #51
//					CF2
//****************************************
class CF2 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	CF2(const char * fname = "CF2", int varsn = 10, short objsn = 2, short ix = 51, short consn = 1)
		: function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~CF2() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};
//****************************************
//				Function #52
//					CF3
//****************************************
class CF3 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	CF3(const char * fname = "CF3", int varsn = 10, short objsn = 2, short ix = 52, short consn = 1)
		: function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~CF3() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};
//****************************************
//				Function #53
//					CF4
//****************************************
class CF4 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	CF4(const char * fname = "CF4", int varsn = 10, short objsn = 2, short ix = 53, short consn = 1)
		: function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~CF4() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};
//****************************************
//				Function #54
//					CF5
//****************************************
class CF5 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	CF5(const char * fname = "CF5", int varsn = 10, short objsn = 2, short ix = 54, short consn = 1)
		: function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~CF5() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};
//****************************************
//				Function #55
//					CF6
//****************************************
class CF6 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	CF6(const char * fname = "CF6", int varsn = 10, short objsn = 2, short ix = 55, short consn = 2)
		: function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~CF6() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};
//****************************************
//				Function #56
//					CF7
//****************************************
class CF7 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	CF7(const char * fname = "CF7", int varsn = 10, short objsn = 2, short ix = 56, short consn = 2)
		: function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~CF7() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};

//****************************************
//				Function #57
//					CF8
//****************************************
class CF8 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	CF8(const char * fname = "CF8", int varsn = 10, short objsn = 3, short ix = 57, short consn = 1)
		: function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~CF8() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};

//****************************************
//				Function #58
//					CF9
//****************************************
class CF9 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	CF9(const char * fname = "CF9", int varsn = 10, short objsn = 3, short ix = 58, short consn = 1)
		: function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~CF9() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};

//****************************************
//				Function #59
//					CF10
//****************************************
class CF10 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	CF10(const char * fname = "CF10", int varsn = 10, short objsn = 3, short ix = 59, short consn = 1)
		: function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~CF10() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};

//****************************************
//				Function #60
//	Deb–Thiele-Laumanns-Zitzler–'s function 8
//****************************************
class DTLZ8 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	DTLZ8(const char * fname = "DTLZ8", int varsn = n_func_obj*10, short objsn = n_func_obj, short ix = 60, short consn = n_func_obj)
		: function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~DTLZ8() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};

//****************************************
//				Function #61
//	Deb–Thiele-Laumanns-Zitzler–'s function 9
//****************************************
class DTLZ9 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	DTLZ9(const char * fname = "DTLZ9", int varsn = n_func_obj * 10, short objsn = n_func_obj, short ix = 61, short consn = n_func_obj - 1)
		: function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~DTLZ9() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};

//****************************************
//				Function #62
//					IMB11
//****************************************
class IMB11 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	IMB11(const char * fname = "IMB11", int varsn = 10, short objsn = 2, short ix = 62, short consn = 1)
		: function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~IMB11() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};

//****************************************
//				Function #63
//					IMB12
//****************************************
class IMB12 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	IMB12(const char * fname = "IMB12", int varsn = 10, short objsn = 2, short ix = 63, short consn = 1)
		: function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~IMB12() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};

//****************************************
//				Function #64
//					IMB13
//****************************************
class IMB13 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	IMB13(const char * fname = "IMB13", int varsn = 10, short objsn = 2, short ix = 64, short consn = 1)
		: function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~IMB13() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};

//****************************************
//				Function #65
//					IMB14
//****************************************
class IMB14 : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	IMB14(const char * fname = "IMB14", int varsn = 10, short objsn = 3, short ix = 65, short consn = 1)
		: function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~IMB14() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};


/********************************************
		DYNAMIC UNCONSTRAINED FUNCTIONS
*********************************************/

//****************************************
//				Function #66
//					UDF1
//****************************************
class UDF1 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UDF1(const char * fname = "UDF1", int varsn = 10, short objsn = 2, short ix = 66, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~UDF1() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};

//****************************************
//				Function #67
//					UDF2
//****************************************
class UDF2 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UDF2(const char * fname = "UDF2", int varsn = 10, short objsn = 2, short ix = 67, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~UDF2() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};

//****************************************
//				Function #68
//					UDF3
//****************************************
class UDF3 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UDF3(const char * fname = "UDF3", int varsn = 10, short objsn = 2, short ix = 68, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~UDF3() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};

//****************************************
//				Function #69
//					UDF4
//****************************************
class UDF4 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UDF4(const char * fname = "UDF4", int varsn = 10, short objsn = 2, short ix = 69, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~UDF4() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};
//****************************************
//				Function #70
//					UDF5
//****************************************
class UDF5 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UDF5(const char * fname = "UDF5", int varsn = 10, short objsn = 2, short ix = 70, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~UDF5() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};
//****************************************
//				Function #71
//					UDF6
//****************************************
class UDF6 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UDF6(const char * fname = "UDF6", int varsn = 10, short objsn = 2, short ix = 71, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~UDF6() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};
//****************************************
//				Function #72
//					UDF8
//****************************************
class UDF8 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UDF8(const char * fname = "UDF8", int varsn = 10, short objsn = 2, short ix = 72, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
		T_Vector_Create(5);
	};
	~UDF8() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};
//****************************************
//				Function #73
//					UDF9
//****************************************
class UDF9 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UDF9(const char * fname = "UDF9", int varsn = 10, short objsn = 2, short ix = 73, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
		T_Vector_Create(5);
	};
	~UDF9() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};
//****************************************
//				Function #74
//					JY1
//****************************************
class JY1 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	JY1(const char * fname = "JY1", int varsn = 10, short objsn = 2, short ix = 74, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~JY1() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};
//****************************************
//				Function #75
//					JY2
//****************************************
class JY2 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	JY2(const char * fname = "JY2", int varsn = 10, short objsn = 2, short ix = 75, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~JY2() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};
//****************************************
//				Function #76
//					JY3
//****************************************
class JY3 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	JY3(const char * fname = "JY3", int varsn = 10, short objsn = 2, short ix = 76, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~JY3() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};
//****************************************
//				Function #77
//					JY4
//****************************************
class JY4 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	JY4(const char * fname = "JY4", int varsn = 10, short objsn = 2, short ix = 77, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~JY4() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};
//****************************************
//				Function #78
//					JY5
//****************************************
class JY5 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	JY5(const char * fname = "JY5", int varsn = 10, short objsn = 2, short ix = 78, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~JY5() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};
//****************************************
//				Function #79
//					JY6
//****************************************
class JY6 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	JY6(const char * fname = "JY6", int varsn = 10, short objsn = 2, short ix = 79, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~JY6() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};
//****************************************
//				Function #80
//					JY7
//****************************************
class JY7 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	JY7(const char * fname = "JY7", int varsn = 10, short objsn = 2, short ix = 80, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~JY7() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};
//****************************************
//				Function #81
//					JY8
//****************************************
class JY8 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	JY8(const char * fname = "JY8", int varsn = 10, short objsn = 2, short ix = 81, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~JY8() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};
//****************************************
//				Function #82
//					JY9
//****************************************
class JY9 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	JY9(const char * fname = "JY9", int varsn = 10, short objsn = 2, short ix = 82, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~JY9() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};
//****************************************
//				Function #83
//					JY10
//****************************************
class JY10 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	JY10(const char * fname = "JY10", int varsn = 10, short objsn = 2, short ix = 83, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~JY10() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};

//****************************************
//				Function #84
//					FDA1
//****************************************
class FDA1 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	FDA1(const char * fname = "FDA1", int varsn = 20, short objsn = 2, short ix = 84, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~FDA1() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};

//****************************************
//				Function #85
//					FDA2
//****************************************
class FDA2 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	FDA2(const char * fname = "FDA2", int varsn = 31, short objsn = 2, short ix = 85, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~FDA2() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};

//****************************************
//				Function #86
//					FDA3
//****************************************
class FDA3 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	FDA3(const char * fname = "FDA3", int varsn = 30, short objsn = 2, short ix = 86, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~FDA3() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};

//****************************************
//				Function #87
//					FDA4
//****************************************
class FDA4 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	FDA4(const char * fname = "FDA4", int varsn = 11, short objsn = 2, short ix = 87, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~FDA4() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};

//****************************************
//				Function #88
//					FDA5
//****************************************
class FDA5 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	FDA5(const char * fname = "FDA5", int varsn = 11, short objsn = 2, short ix = 88, short consn = 0)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~FDA5() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
};

/********************************************
		DYNAMIC CONSTRAINED FUNCTIONS
*********************************************/

//****************************************
//				Function #89
//					CDF1
//****************************************
class CDF1 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	CDF1(const char * fname = "CDF1", int varsn = 10, short objsn = 2, short ix = 89, short consn = 1)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~CDF1() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};

//****************************************
//				Function #90
//					CDF2
//****************************************
class CDF2 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	CDF2(const char * fname = "CDF2", int varsn = 10, short objsn = 2, short ix = 90, short consn = 2)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~CDF2() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};

//****************************************
//				Function #91
//					CDF3
//****************************************
class CDF3 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	CDF3(const char * fname = "CDF3", int varsn = 10, short objsn = 2, short ix = 91, short consn = 1)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~CDF3() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};

//****************************************
//				Function #92
//					CDF4
//****************************************
class CDF4 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	CDF4(const char * fname = "CDF4", int varsn = 10, short objsn = 2, short ix = 92, short consn = 1)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~CDF4() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};
//****************************************
//				Function #93
//					CDF5
//****************************************
class CDF5 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	CDF5(const char * fname = "CDF5", int varsn = 10, short objsn = 2, short ix = 93, short consn = 1)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~CDF5() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};
//****************************************
//				Function #94
//					CDF6
//****************************************
class CDF6 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	CDF6(const char * fname = "CDF6", int varsn = 10, short objsn = 2, short ix = 94, short consn = 1)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
		T_Vector_Create(5);
	};
	~CDF6() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};
//****************************************
//				Function #95
//					CDF7
//****************************************
class CDF7 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	CDF7(const char * fname = "CDF7", int varsn = 10, short objsn = 2, short ix = 95, short consn = 1)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~CDF7() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};
//****************************************
//				Function #96
//					CDF8
//****************************************
class CDF8 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	CDF8(const char * fname = "CDF8", int varsn = 10, short objsn = 2, short ix = 96, short consn = 1)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~CDF8() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};

//****************************************
//				Function #97
//					CDF9
//****************************************
class CDF9 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	CDF9(const char * fname = "CDF9", int varsn = 10, short objsn = 2, short ix = 97, short consn = 1)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~CDF9() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};
//****************************************
//				Function #98
//					CDF10
//****************************************
class CDF10 : public dynamic_function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	CDF10(const char * fname = "CDF10", int varsn = 10, short objsn = 2, short ix = 98, short consn = 1)
		: dynamic_function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~CDF10() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};

//****************************************
//				Function #99
//					GEN
//****************************************
class GEN : public function
{
protected:
	/**Set boundaries for variables**/
	void Bound_Set();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	GEN(const char * fname = "GEN", int varsn = 404668, short objsn = 2, short ix = 59, short consn = 1)
		: function(fname, varsn, objsn, ix, consn) {
		Bound_Set();
	};
	~GEN() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};
#endif // !FIT_FUNCTION_H
