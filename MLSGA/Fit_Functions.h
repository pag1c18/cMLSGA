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


		 FIT FUNCTIOCS header
	Storage of fitness funtions
			and basic funtionc class

*/


#pragma once


#ifndef FIT_FUNCTION_H
#define FIT_FUNCTION_H
#include <string>
#include <iostream>
#include <vector>
#include <fstream>

#include "Struct.h"
#include "Const.h"

/**Template for function storage class*/
class function
{
private:
	std::string name_func;									///<Name of  function
	int num_vars;											///<Number of variables in function
	short num_objs;											///<Number of objectives in function
	short num_cons;											///<Number of constaines
	std::vector<STRUCTURES::boundaries> bound;				///< Boundaries for variables
	std::vector<STRUCTURES::boundaries> max_min_fit;		///< Boundaries for max and min fitness
public: 
	std::vector<double> code_saved;							///<Currently saved code (variables) for dynamic problems
	bool const_speed;										///<If the speed is constant (for VOS)
	double service_speed;									///<Current service speed (for VOS)
	int temp_cons;											///<Current constraint
	short wind_delay;										///<The delay in weather data
	short wind_size;										///<The size of weather data
	std::vector<short> delay_time_matrix;					///<The matrix where delay weather data is stored (for dynamic problems)
protected:
	bool time_dep = false;									///<Dynamic function - 1 Yes, 0 No
	std::vector<float> t_vector;							///<dynamic t variables for changing random functions

protected:
	/**Set boundaries for variables
	
	@param bn - vector of the variables boundaries structure
	*/
	virtual void Bound_Set(std::vector<STRUCTURES::boundaries> & bn) { bound = bn; }
	/**Set boundaries for fitness
	@param mmf vector of the fitness boundaries structure
	*/
	void Max_Min_Fit_Set(std::vector<STRUCTURES::boundaries> & mmf) { max_min_fit = mmf; }	
	/**Create the t_vector for unpredictable functions
	@param n - size of the t_vector
	*/
	void T_Vector_Create(int n) { t_vector = std::vector<float>(n, 0.f); };
public:
	/**Default constructor. Throws ERROR - function class cannot be empty*/
	function() {std::cout << "ERROR#02: FUNCTION - CREATION"; system("pause"); abort();};
	/**Default normal constructor.
	@param fname - function name
	@param varsn - number of variables
	@param objsn - number of objectives
	@param consn - index of constraints
	*/
	function(const char * fname, int varsn, short objsn, short consn) 
	{
		name_func = fname; num_vars = varsn; num_objs = objsn; num_cons = consn;
	};
	~function() {};

	/**Function where the Pareto Optimal Front is defined and saved to file
	@param indexr - index of the run
	@param size - size of the pareto front. Default is defined in const.h
	*/
	virtual std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);																
	/*
	***Function where the Pareto Optimal Front is defined and saved to file for dynamic problems.
	@param indexr - index of the run
	@param t - current time step
	@param size - size of the pareto front. Default is defined in const.h
	*/																																																		
	virtual std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size) { std::vector<std::vector<double>> temp; return temp; };
	/**Returning the boundary for a single variable
	@param i - index of the variable
	@param c - defines which boundary is returned. "lower"/"upper"
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
	/**Returning the boundaries for a single variables as a vector*/
	std::vector<STRUCTURES::boundaries> Bound() { return bound; };
	/**Returning boundaries for a single objective
	@param i - index of the objective
	@param c - defines which boundary is returned. "lower"/"upper"
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
	/**Returning the boundaries for a single objective as a vector*/
	std::vector<STRUCTURES::boundaries> Max_Min_Fit() { return max_min_fit; };													
	/**Returning the number of variables for the current function*/
	int Vars() const { return num_vars; };												
	/**Returning the number of objectives for the current function*/
	short Objs() const { return num_objs; };	
	/**Returning the number of constrains for the current function*/
	short Cons() const { return num_cons; };
	/**Returning time dependency*/
	bool Time_Dep() const { return time_dep; };
	/**Returning the name of the current function*/
	std::string & Name_Show() { return name_func; };	
	/**Fitness (objectives) calculation for the given code (vector of variables).  Returning error in the template.
	@param code - vector with the code.
	*/
	virtual std::vector<double> Fitness_C(const std::vector<double> &code) { std::cout << "ERROR#: FUNCTION - FITNESS_C"; system("pause"); abort(); };
	/**Fitness (objectives) calculation for the given code (vector of variables), for dynamic problems.  Returning error in the template.
	@param code - vector with the code.
	@param tau - current time indicator
	@param t -  current time step
	*/
	virtual std::vector<double> Fitness_C(const std::vector<double> &code, int tau, double t) { std::cout << "ERROR#: FUNCTION - FITNESS_C"; system("pause"); abort(); };
	
	/**Calculate the constraints. For contrained problems: returns true if constrains are violated**
	@param code - vector with the code.
	@param fitness - vector with the fitness (objectives) values
	@param t -  current time step
	*/
	virtual std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0) { abort(); }

	/*
	*Min fitness calculation for one objective optimisation
	@param pareto_front - real pareto front for which fitness will be calculated
	*/
	double Min_Fitness_Get(const std::vector<std::vector<double>> & pareto_front);
	/**Update the t_vector for random dynamic problems
	@param s - if it have to be updated by one step (true) or zeroed (false)
	*/
	void T_Vector_Update(bool s = true);
};

class dynamic_function : public function
{
public:
	/**Default empty constructor. throw ERROR - function class cannot be empty**/
	dynamic_function() { std::cout << "ERROR#02: FUNCTION - CREATION"; system("pause"); abort(); };
	/**Default normal constructor. Calls constructor of fitness class.
	@param fname - function name
	@param varsn - number of variables
	@param objsn - number of objectives
	@param ix - index of the function
	*/
	dynamic_function(const char * fname, int varsn, short objsn, short consn) : function(fname, varsn, objsn, consn)
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
	ZDT1(const char * fname = "ZDT1", int varsn = 30, short objsn = 2, short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	ZDT2(const char * fname = "ZDT2", int varsn = 30, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	ZDT3(const char * fname = "ZDT3", int varsn = 30, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	ZDT4(const char * fname = "ZDT4", int varsn = 10, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	ZDT5(const char * fname = "ZDT5", int varsn = 11, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	ZDT6(const char * fname = "ZDT6", int varsn = 10, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	DTLZ1(const char * fname = "DTLZ1", int varsn = n_func_obj + 4, short objsn = n_func_obj,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	DTLZ2(const char * fname = "DTLZ2", int varsn = n_func_obj + 9, short objsn = n_func_obj,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	DTLZ3(const char * fname = "DTLZ3", int varsn = n_func_obj + 9, short objsn = n_func_obj,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	DTLZ4(const char * fname = "DTLZ4", int varsn = n_func_obj + 9, short objsn = n_func_obj,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	DTLZ5(const char * fname = "DTLZ5", int varsn = n_func_obj + 9, short objsn = n_func_obj,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	DTLZ6(const char * fname = "DTLZ6", int varsn = n_func_obj + 9, short objsn = n_func_obj,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	DTLZ7(const char * fname = "DTLZ7", int varsn = n_func_obj + 19, short objsn = n_func_obj,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	MOP1(const char * fname = "MOP1", int varsn = 10, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	MOP2(const char * fname = "MOP2", int varsn = 10, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	MOP3(const char * fname = "MOP3", int varsn = 10, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	MOP4(const char * fname = "MOP4", int varsn = 10, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	MOP5(const char * fname = "MOP5", int varsn = 10, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	MOP6(const char * fname = "MOP6", int varsn = 10, short objsn = 3,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	MOP7(const char * fname = "MOP7", int varsn = 10, short objsn = 3,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	UF1(const char * fname = "UF1", int varsn = 30, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn) {
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
	UF2(const char * fname = "UF2", int varsn = 30, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	UF3(const char * fname = "UF3", int varsn = 30, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	UF4(const char * fname = "UF4", int varsn = 30, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	UF5(const char * fname = "UF5", int varsn = 30, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	UF6(const char * fname = "UF6", int varsn = 30, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	UF7(const char * fname = "UF7", int varsn = 30, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	UF8(const char * fname = "UF8", int varsn = 30, short objsn = 3,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	UF9(const char * fname = "UF9", int varsn = 30, short objsn = 3,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	UF10(const char * fname = "UF10", int varsn = 30, short objsn = 3,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	WFG1(const char * fname = "WFG1", int varsn = (2 * (n_func_obj - 1) + 20), short objsn = n_func_obj,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	WFG2(const char * fname = "WFG2", int varsn = (2 * (n_func_obj - 1) + 20), short objsn = n_func_obj,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	WFG3(const char * fname = "WFG3", int varsn = (2 * (n_func_obj - 1) + 20), short objsn = n_func_obj,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	WFG4(const char * fname = "WFG4", int varsn = (2 * (n_func_obj - 1) + 20), short objsn = n_func_obj,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	WFG5(const char * fname = "WFG5", int varsn = (2 * (n_func_obj - 1) + 20), short objsn = n_func_obj,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	WFG6(const char * fname = "WFG6", int varsn = (2 * (n_func_obj - 1) + 20), short objsn = n_func_obj,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	WFG7(const char * fname = "WFG7", int varsn = (2 * (n_func_obj - 1) + 20), short objsn = n_func_obj,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	WFG8(const char * fname = "WFG8", int varsn = (2 * (n_func_obj - 1) + 20), short objsn = n_func_obj,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	WFG9(const char * fname = "WFG9", int varsn = (2 * (n_func_obj - 1) + 20), short objsn = n_func_obj,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	IMB1(const char * fname = "IMB1", int varsn = 10, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	IMB2(const char * fname = "IMB2", int varsn = 10, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	IMB3(const char * fname = "IMB3", int varsn = 10, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	IMB4(const char * fname = "IMB4", int varsn = 10, short objsn = 3,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	IMB5(const char * fname = "IMB5", int varsn = 10, short objsn = 3,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	IMB6(const char * fname = "IMB6", int varsn = 10, short objsn = 3,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	IMB7(const char * fname = "IMB7", int varsn = 10, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	IMB8(const char * fname = "IMB8", int varsn = 10, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	IMB9(const char * fname = "IMB9", int varsn = 10, short objsn = 2,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	IMB10(const char * fname = "IMB10", int varsn = 10, short objsn = 3,  short consn = 0)
		: function(fname, varsn, objsn, consn)
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
	CF1(const char * fname = "CF1", int varsn = 10, short objsn = 2,  short consn = 1)
		: function(fname, varsn, objsn, consn) {
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
	CF2(const char * fname = "CF2", int varsn = 10, short objsn = 2,  short consn = 1)
		: function(fname, varsn, objsn, consn) {
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
	CF3(const char * fname = "CF3", int varsn = 10, short objsn = 2,  short consn = 1)
		: function(fname, varsn, objsn, consn) {
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
	CF4(const char * fname = "CF4", int varsn = 10, short objsn = 2,  short consn = 1)
		: function(fname, varsn, objsn, consn) {
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
	CF5(const char * fname = "CF5", int varsn = 10, short objsn = 2,  short consn = 1)
		: function(fname, varsn, objsn, consn) {
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
	CF6(const char * fname = "CF6", int varsn = 10, short objsn = 2,  short consn = 2)
		: function(fname, varsn, objsn, consn) {
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
	CF7(const char * fname = "CF7", int varsn = 10, short objsn = 2,  short consn = 2)
		: function(fname, varsn, objsn, consn) {
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
	CF8(const char * fname = "CF8", int varsn = 10, short objsn = 3,  short consn = 1)
		: function(fname, varsn, objsn, consn) {
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
	CF9(const char * fname = "CF9", int varsn = 10, short objsn = 3,  short consn = 1)
		: function(fname, varsn, objsn, consn) {
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
	CF10(const char * fname = "CF10", int varsn = 10, short objsn = 3,  short consn = 1)
		: function(fname, varsn, objsn, consn) {
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
	DTLZ8(const char * fname = "DTLZ8", int varsn = n_func_obj*10, short objsn = n_func_obj,  short consn = n_func_obj)
		: function(fname, varsn, objsn, consn) {
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
	DTLZ9(const char * fname = "DTLZ9", int varsn = n_func_obj * 10, short objsn = n_func_obj,  short consn = n_func_obj - 1)
		: function(fname, varsn, objsn, consn) {
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
	IMB11(const char * fname = "IMB11", int varsn = 10, short objsn = 2,  short consn = 1)
		: function(fname, varsn, objsn, consn) {
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
	IMB12(const char * fname = "IMB12", int varsn = 10, short objsn = 2,  short consn = 1)
		: function(fname, varsn, objsn, consn) {
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
	IMB13(const char * fname = "IMB13", int varsn = 10, short objsn = 2,  short consn = 1)
		: function(fname, varsn, objsn, consn) {
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
	IMB14(const char * fname = "IMB14", int varsn = 10, short objsn = 3,  short consn = 1)
		: function(fname, varsn, objsn, consn) {
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

//****************************************
//				Function #66
//					DAS_CMOP1
//****************************************
class DAS_CMOP1 : public function
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
	DAS_CMOP1(const char* fname = "DAS_CMOP1", int varsn = 30, short objsn = 2,  short consn = 11)
		: function(fname, varsn, objsn, consn) {
		Bound_Set();
	};
	~DAS_CMOP1() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double>& code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};
//****************************************
//				Function #67
//					DAS_CMOP2
//****************************************
class DAS_CMOP2 : public function
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
	DAS_CMOP2(const char* fname = "DAS_CMOP2", int varsn = 30, short objsn = 2,  short consn = 11)
		: function(fname, varsn, objsn, consn) {
		Bound_Set();
	};
	~DAS_CMOP2() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double>& code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};
//****************************************
//				Function #68
//					DAS_CMOP3
//****************************************
class DAS_CMOP3 : public function
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
	DAS_CMOP3(const char* fname = "DAS_CMOP3", int varsn = 30, short objsn = 2,  short consn = 11)
		: function(fname, varsn, objsn, consn) {
		Bound_Set();
	};
	~DAS_CMOP3() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double>& code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};
//****************************************
//				Function #69
//					DAS_CMOP4
//****************************************
class DAS_CMOP4 : public function
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
	DAS_CMOP4(const char* fname = "DAS_CMOP4", int varsn = 30, short objsn = 2,  short consn = 11)
		: function(fname, varsn, objsn, consn) {
		Bound_Set();
	};
	~DAS_CMOP4() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double>& code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};
//****************************************
//				Function #70
//					DAS_CMOP5
//****************************************
class DAS_CMOP5 : public function
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
	DAS_CMOP5(const char* fname = "DAS_CMOP5", int varsn = 30, short objsn = 2,  short consn = 11)
		: function(fname, varsn, objsn, consn) {
		Bound_Set();
	};
	~DAS_CMOP5() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double>& code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};
//****************************************
//				Function #71
//					DAS_CMOP6
//****************************************
class DAS_CMOP6 : public function
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
	DAS_CMOP6(const char* fname = "DAS_CMOP6", int varsn = 30, short objsn = 2,  short consn = 11)
		: function(fname, varsn, objsn, consn) {
		Bound_Set();
	};
	~DAS_CMOP6() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double>& code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};
//****************************************
//				Function #72
//					DAS_CMOP7
//****************************************
class DAS_CMOP7 : public function
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
	DAS_CMOP7(const char* fname = "DAS_CMOP7", int varsn = 30, short objsn = 3,  short consn = 7)
		: function(fname, varsn, objsn, consn) {
		Bound_Set();
	};
	~DAS_CMOP7() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double>& code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};

//****************************************
//				Function #73
//					DAS_CMOP8
//****************************************
class DAS_CMOP8 : public function
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
	DAS_CMOP8(const char* fname = "DAS_CMOP8", int varsn = 30, short objsn = 3,  short consn = 7)
		: function(fname, varsn, objsn, consn) {
		Bound_Set();
	};
	~DAS_CMOP8() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double>& code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};

//****************************************
//				Function #74
//					DAS_CMOP9
//****************************************
class DAS_CMOP9 : public function
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
	DAS_CMOP9(const char* fname = "DAS_CMOP9", int varsn = 30, short objsn = 3,  short consn = 7)
		: function(fname, varsn, objsn, consn) {
		Bound_Set();
	};
	~DAS_CMOP9() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double>& code);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};



/********************************************
		DYNAMIC UNCONSTRAINED FUNCTIONS
*********************************************/

//****************************************
//				Function #75
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
	UDF1(const char * fname = "UDF1", int varsn = 10, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #76
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
	UDF2(const char * fname = "UDF2", int varsn = 10, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #77
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
	UDF3(const char * fname = "UDF3", int varsn = 10, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #78
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
	UDF4(const char * fname = "UDF4", int varsn = 10, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #79
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
	UDF5(const char * fname = "UDF5", int varsn = 10, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #80
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
	UDF6(const char * fname = "UDF6", int varsn = 10, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #81
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
	UDF8(const char * fname = "UDF8", int varsn = 10, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #82
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
	UDF9(const char * fname = "UDF9", int varsn = 10, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #83
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
	JY1(const char * fname = "JY1", int varsn = 10, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #84
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
	JY2(const char * fname = "JY2", int varsn = 10, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #85
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
	JY3(const char * fname = "JY3", int varsn = 10, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #86
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
	JY4(const char * fname = "JY4", int varsn = 10, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #87
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
	JY5(const char * fname = "JY5", int varsn = 10, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #88
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
	JY6(const char * fname = "JY6", int varsn = 10, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #89
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
	JY7(const char * fname = "JY7", int varsn = 10, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #90
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
	JY8(const char * fname = "JY8", int varsn = 10, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #90
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
	JY9(const char * fname = "JY9", int varsn = 10, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #92
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
	JY10(const char * fname = "JY10", int varsn = 10, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #93
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
	FDA1(const char * fname = "FDA1", int varsn = 20, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #94
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
	FDA2(const char * fname = "FDA2", int varsn = 31, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #95
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
	FDA3(const char * fname = "FDA3", int varsn = 30, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #96
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
	FDA4(const char * fname = "FDA4", int varsn = 11, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #97
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
	FDA5(const char * fname = "FDA5", int varsn = 11, short objsn = 2,  short consn = 0)
		: dynamic_function(fname, varsn, objsn, consn) {
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
//				Function #98
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
	CDF1(const char* fname = "CDF1", int varsn = 10, short objsn = 2,  short consn = 2)
		: dynamic_function(fname, varsn, objsn, consn) {
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
	std::vector<double> Fitness_C(const std::vector<double>& code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};
//****************************************
//				Function #99
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
	CDF2(const char* fname = "CDF2", int varsn = 10, short objsn = 2,  short consn = 1)
		: dynamic_function(fname, varsn, objsn, consn) {
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
	std::vector<double> Fitness_C(const std::vector<double>& code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};
//****************************************
//				Function #100
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
	CDF3(const char* fname = "CDF3", int varsn = 10, short objsn = 2,  short consn = 1)
		: dynamic_function(fname, varsn, objsn, consn) {
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
	std::vector<double> Fitness_C(const std::vector<double>& code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};
//****************************************
//				Function #101
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
	CDF4(const char* fname = "CDF4", int varsn = 10, short objsn = 2,  short consn = 1)
		: dynamic_function(fname, varsn, objsn, consn) {
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
	std::vector<double> Fitness_C(const std::vector<double>& code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};
//****************************************
//				Function #102
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
	CDF5(const char* fname = "CDF5", int varsn = 10, short objsn = 2,  short consn = 1)
		: dynamic_function(fname, varsn, objsn, consn) {
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
	std::vector<double> Fitness_C(const std::vector<double>& code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};
//****************************************
//				Function #103
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
	CDF6(const char* fname = "CDF6", int varsn = 10, short objsn = 2,  short consn = 2)
		: dynamic_function(fname, varsn, objsn, consn) {
		Bound_Set();
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
	std::vector<double> Fitness_C(const std::vector<double>& code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};
//****************************************
//				Function #104
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
	CDF7(const char * fname = "CDF7", int varsn = 10, short objsn = 2,  short consn = 1)
		: dynamic_function(fname, varsn, objsn, consn) {
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
	*/
	std::vector<double> Cons_Calc(const std::vector<double> &code, const std::vector<double> &fit, double t = 0.0);
};

//****************************************
//				Function #105
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
	CDF8(const char* fname = "CDF8", int varsn = 10, short objsn = 2,  short consn = 1)
		: dynamic_function(fname, varsn, objsn, consn) {
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
	std::vector<double> Fitness_C(const std::vector<double>& code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};
//****************************************
//				Function #106
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
	CDF9(const char* fname = "CDF9", int varsn = 10, short objsn = 2,  short consn = 2)
		: dynamic_function(fname, varsn, objsn, consn) {
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
	std::vector<double> Fitness_C(const std::vector<double>& code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};
//****************************************
//				Function #107
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
	CDF10(const char* fname = "CDF10", int varsn = 10, short objsn = 2,  short consn = 2)
		: dynamic_function(fname, varsn, objsn, consn) {
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
	std::vector<double> Fitness_C(const std::vector<double>& code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};
//****************************************
//				Function #108
//					CDF11
//****************************************
class CDF11 : public dynamic_function
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
	CDF11(const char* fname = "CDF11", int varsn = 10, short objsn = 2,  short consn = 1)
		: dynamic_function(fname, varsn, objsn, consn) {
		Bound_Set();
	};
	~CDF11() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double>& code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};

//****************************************
//				Function #109
//					CDF12
//****************************************
class CDF12 : public dynamic_function
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
	CDF12(const char* fname = "CDF12", int varsn = 10, short objsn = 2,  short consn = 1)
		: dynamic_function(fname, varsn, objsn, consn) {
		Bound_Set();
	};
	~CDF12() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double>& code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};
//****************************************
//				Function #110
//					CDF13
//****************************************
class CDF13 : public dynamic_function
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
	CDF13(const char* fname = "CDF13", int varsn = 10, short objsn = 2,  short consn = 1)
		: dynamic_function(fname, varsn, objsn, consn) {
		Bound_Set();
		T_Vector_Create(5);
	};
	~CDF13() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double>& code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};


//****************************************
//				Function #111
//					CDF14
//****************************************
class CDF14 : public dynamic_function
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
	CDF14(const char* fname = "CDF14", int varsn = 10, short objsn = 2,  short consn = 1)
		: dynamic_function(fname, varsn, objsn, consn) {
		Bound_Set();
	};
	~CDF14() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double>& code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};
//****************************************
//				Function #112
//					CDF15
//****************************************
class CDF15 : public dynamic_function
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
	CDF15(const char* fname = "CDF15", int varsn = 10, short objsn = 2,  short consn = 1)
		: dynamic_function(fname, varsn, objsn, consn) {
		Bound_Set();
	};
	~CDF15() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<double>> Plot_PF(int indexr, double t, int size = PF_real_size);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double>& code, int tau, double t);
	/*
	*Constrain check - for contrained problems: returns true if constrains are violated*
	@param code - vector with the code
	@param fitness - vector with the fitness values
	@param t -  current time
	*/
	std::vector<double> Cons_Calc(const std::vector<double>& code, const std::vector<double>& fit, double t = 0.0);
};



//****************************************
//				Function #201
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
	GEN(const char * fname = "GEN", int varsn = 3678, short objsn = 3,  short consn = 0)
		: function(fname, varsn, objsn, consn) {
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
};

//****************************************
//				Function #202
//					GEN_pat
//****************************************
class GEN_pat : public function
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
	GEN_pat(const char * fname = "GEN_pat", int varsn = 3678, short objsn = 3,  short consn = 0)
		: function(fname, varsn, objsn, consn) {
		Bound_Set();
	};
	~GEN_pat() {};
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
#endif // !FIT_FUNCTION_H
