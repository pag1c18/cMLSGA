#pragma once

//****************************************
//				SELECTION header
//	Storage of different selection types
//			and basic selection class
//****************************************

#ifndef FUNCTION_H
#define FUNCTION_H
#include <string>
#include <ctime>
#include "Struct.h"
#include <iostream>
#include <vector>
#include "Const.h"
#include <fstream>


class function
{
private:
	std::string name_func;									//Name of  function
	short num_vars;											//Number of variables in function
	short num_objs;											//Number of objectives in function
	short index;											//Function intex
	std::vector<STRUCTURES::boundaries> bound;				// Boundaries for variables
	std::vector<STRUCTURES::boundaries> max_min_fit;		// Boundaries for max and min fitness
protected:
	/*
	*Set boundaries for variables*
	@param bn vector of the variables boundaries structure
	*/
	virtual void Set_Bound(std::vector<STRUCTURES::boundaries> & bn) { bound = bn; }
	/*
	*Set boundaries*
	@param mmf vector of the fitness boundaries structure
	*/
	void Set_Max_Min_Fit(std::vector<STRUCTURES::boundaries> & mmf) { max_min_fit = mmf; }	
public:
	/**Default constructor, throw ERROR - function class cannot be empty**/
	function() {std::cout << "ERROR#02: FUNCTION - CREATION"; abort();};
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	function(const char * fname, short varsn, short objsn, short ix) 
	{ name_func = fname; num_vars = varsn; num_objs = objsn; index = ix; };
	~function() {};

	/*
	*Pareto Front plotting and saving to file*
	@param indexr index of the run
	*/
	virtual std::vector<std::vector<float>> Plot_PF(short indexr) { std::vector<std::vector<float>> temp; return temp; };																//Plotting real PF and saving it to file
	
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
			abort();
		}
		if (c == "lower")
			return bound[i].lower;
		else
			return bound[i].upper;
	}
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
	short Vars() const { return num_vars; };												
	/**Returning the number of objectives for the current function**/
	short Objs() const { return num_objs; };												
	/**Returning the name of the current function**/
	std::string & Show_Name() { return name_func; };	
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	virtual std::vector<double> Fitness_C(const std::vector<double> &code) = 0;
};
//****************************************
//				Function #1
//					UF1
//****************************************
class UF1 : public function
{
protected:
	/**Set boundaries for variables**/
	void Set_Bound();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UF1(const char * fname = "UF1", short varsn = 30, short objsn = 2, short ix = 1)
		: function(fname, varsn, objsn, ix) {
		Set_Bound();
	};
	~UF1() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<float>> Plot_PF(short indexr);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
//****************************************
//				Function #2
//					UF2
//****************************************
class UF2 : public function
{
protected:
	/**Set boundaries for variables**/
	void Set_Bound();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UF2(const char * fname = "UF2", short varsn = 30, short objsn = 2, short ix = 2)
		: function(fname, varsn, objsn, ix)
	{
		Set_Bound();
	};
	~UF2() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<float>> Plot_PF(short indexr);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C( const std::vector<double> &code);
};

//****************************************
//				Function #3
//					UF3
//****************************************
class UF3 : public function
{
protected:
	/**Set boundaries for variables**/
	void Set_Bound();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UF3(const char * fname = "UF3", short varsn = 30, short objsn = 2, short ix = 3)
		: function(fname, varsn, objsn, ix)
	{
		Set_Bound();
	};
	~UF3() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<float>> Plot_PF(short indexr);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//				Function #4
//					UF4
//****************************************
class UF4 : public function
{
protected:
	/**Set boundaries for variables**/
	void Set_Bound();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UF4(const char * fname = "UF4", short varsn = 30, short objsn = 2, short ix = 4)
		: function(fname, varsn, objsn, ix)
	{
		Set_Bound();
	};
	~UF4() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<float>> Plot_PF(short indexr);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//				Function #5
//					UF5
//****************************************
class UF5 : public function
{
protected:
	/**Set boundaries for variables**/
	void Set_Bound();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UF5(const char * fname = "UF5", short varsn = 30, short objsn = 2, short ix = 5)
		: function(fname, varsn, objsn, ix)
	{
		Set_Bound();
	};
	~UF5() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<float>> Plot_PF(short indexr);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//				Function #6
//					UF7
//****************************************
class UF7 : public function
{
protected:
	/**Set boundaries for variables**/
	void Set_Bound();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	UF7(const char * fname = "UF6", short varsn = 30, short objsn = 2, short ix = 6)
		: function(fname, varsn, objsn, ix)
	{
		Set_Bound();
	};
	~UF7() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<float>> Plot_PF(short indexr);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};



//****************************************
//				Function #7
//			Schaffer function 1
//****************************************
class SCH : public function									
{
protected:
	/**Set boundaries for variables**/
	void Set_Bound();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	SCH(const char * fname = "SCH", short varsn = 1, short objsn = 2, short ix = 7)
		: function(fname, varsn, objsn, ix)
	{
		Set_Bound();
	};
	~SCH() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<float>> Plot_PF(short indexr);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//			Function #8
//		Fonseca and Fleming function 
//****************************************
class FON : public function									
{
protected:
	/**Set boundaries for variables**/
	void Set_Bound();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	FON(const char * fname = "FON", short varsn = 3, short objsn = 2, short ix = 8)
		: function(fname, varsn, objsn, ix)
	{
		Set_Bound();
	};
	~FON() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<float>> Plot_PF(short indexr);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//			Function #9
//		Zitzler–Deb–Thiele's function 1
//****************************************
class ZDT1 : public function									
{
protected:
	/**Set boundaries for variables**/
	void Set_Bound();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	ZDT1(const char * fname = "ZDT1", short varsn = 30, short objsn = 2, short ix = 9)
		: function(fname, varsn, objsn, ix)
	{
		Set_Bound();
	};
	~ZDT1() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<float>> Plot_PF(short indexr);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//			Function #10
//		Zitzler–Deb–Thiele's function 2
//****************************************
class ZDT2 : public function
{
protected:
	/**Set boundaries for variables**/
	void Set_Bound();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	ZDT2(const char * fname = "ZDT2", short varsn = 30, short objsn = 2, short ix = 10)
		: function(fname, varsn, objsn, ix)
	{
		Set_Bound();
	};
	~ZDT2() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<float>> Plot_PF(short indexr);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//			Function #11
//		Zitzler–Deb–Thiele's function 3
//****************************************
class ZDT3 : public function
{
protected:
	/**Set boundaries for variables**/
	void Set_Bound();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	ZDT3(const char * fname = "ZDT3", short varsn = 30, short objsn = 2, short ix = 11)
		: function(fname, varsn, objsn, ix)
	{
		Set_Bound();
	};
	~ZDT3() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<float>> Plot_PF(short indexr);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//			Function #12
//		Zitzler–Deb–Thiele's function 4
//****************************************
class ZDT4 : public function
{
protected:
	/**Set boundaries for variables**/
	void Set_Bound();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	ZDT4(const char * fname = "ZDT4", short varsn = 10, short objsn = 2, short ix = 12)
		: function(fname, varsn, objsn, ix)
	{
		Set_Bound();
	};
	~ZDT4() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<float>> Plot_PF(short indexr);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//			Function #13
//		Zitzler–Deb–Thiele's function 5
//****************************************
class ZDT5 : public function
{
protected:
	/**Set boundaries for variables**/
	void Set_Bound();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	ZDT5(const char * fname = "ZDT5", short varsn = 11, short objsn = 2, short ix = 13)
		: function(fname, varsn, objsn, ix)
	{
		Set_Bound();
	};
	~ZDT5() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<float>> Plot_PF(short indexr);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};

//****************************************
//			Function #14
//		Zitzler–Deb–Thiele's function 6
//****************************************
class ZDT6 : public function									
{
protected:
	/**Set boundaries for variables**/
	void Set_Bound();
public:
	/*
	*Default normal constructor*
	@param fname function name
	@param varsn number of variables
	@param objsn number of objectives
	@param ix index of the function
	*/
	ZDT6(const char * fname = "ZDT6", short varsn = 10, short objsn = 2, short ix = 14)
		: function(fname, varsn, objsn, ix)
	{
		Set_Bound();
	};
	~ZDT6() {};
	/*
	*Plotting Pareto Front and saving to vector*
	@param indexr index of the current run
	*/
	std::vector<std::vector<float>> Plot_PF(short indexr);
	/*
	*Fitness calculation for the given code*
	@param code vector with the code
	*/
	std::vector<double> Fitness_C(const std::vector<double> &code);
};
#endif // !FUNCTION_H
