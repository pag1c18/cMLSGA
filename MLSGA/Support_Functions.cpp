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


#include "Support_Functions.h"
#include "Class.h"
#include <sstream>
#include <iomanip>

std::vector<std::vector<double>> reference_point(int pop_size, short nobj);
std::vector<std::vector<double>> recursion_point(int pop_size, short nobj);
int Get_Weight_Num(int layer, int od, int wnum);
void Get_Weight(int layer, int od, double step, std::vector<std::vector<double>> &out);
void Points_Generate_Recursion(double step, short n_obj, std::vector<std::vector<double>>  & out, std::vector<double> short_out);



std::string String_Prec(float v, int p)
{
	std::ostringstream output;			//output stream

										//send values with precision to the stream
	output << std::fixed << std::setprecision(p) << v;

	//return the string
	return output.str();

}

/**
* Calculation of the distance between 2 points.*/
template <typename tname>
double Distance<tname>(const std::vector<tname> & vec1, const std::vector<tname> & vec2)
{
	double dist = 0;		//distance - output
	//check if the size is the same
	if (vec1.size() != vec2.size())
		abort();

	//calculate the distance between 2 points
	for (unsigned int i = 0; i < vec1.size(); i++)
	{
		dist += (double)pow(vec1[i] - (double)vec2[i], 2);
	}
	return sqrt(dist);
}

/**
* Calculation of the distance between 2 points with 2 diffent values types.*/
template <typename tname, typename tname2>
double Distance2<tname,tname2>(const std::vector<tname> & vec1, const std::vector<tname2> & vec2)
{
	//check if the size is the same
	if (vec1.size() != vec2.size())
		abort();

	double dist = 0;		//distance - output

	//calculate the distance between 2 points
	for (unsigned int i = 0; i < vec1.size(); i++)
	{
		dist += (double)pow(vec1[i] - (double)vec2[i], 2);
	}
	return sqrt(dist);
}

/**Calculation of sign**/
template <typename tname>
short sgn(tname val)
{
	if (val > 0.)
		return 1;
	else if (val == 0.)
		return 0;
	else
		return -1;
}

/*
*Check Dominance*
@param ind1 - first individual
@param ind2 - 2nd individual
Return:
1 - ind1 dominates
-1 - ind2 dominates
0 - both nondominated
*/
int Dominance_Check(individual &ind1, individual &ind2, std::vector<short> &fit_indexes)
{
	int flag1, flag2;		//flags
	flag1 = flag2 = 0;
	int nobj = ind1.Fitness_Show().size();	//number of objectives

	int cons_size = ind1.Cons_Show().size();		//nubmer of constrains
	if (cons_size > 0 && !PENALTY_BASED_CONSTRAINTS)
	{
		double constr_val_ind1 = 0., constr_val_ind2 = 0.;	//temporary values of constrains valus sum
		for (int i = 0; i < cons_size; i++)
		{
			if (ind1.Cons_Show()[i] < 0.)
				constr_val_ind1 += ind1.Cons_Show()[i];
			if (ind2.Cons_Show()[i] < 0.)
				constr_val_ind2 += ind2.Cons_Show()[i];
		}
		if (constr_val_ind1 < 0. && constr_val_ind2 < 0.)
		{
			if (constr_val_ind1 > constr_val_ind2)
				return 1;
			else if (constr_val_ind1 < constr_val_ind2)
				return -1;
			else
				return 0;
		}
		else if (constr_val_ind1 < 0. && constr_val_ind2 == 0.)
			return -1;
		else if (constr_val_ind1 == 0. && constr_val_ind2 < 0.)
			return 1;
	}


	for (int i = 0; i < fit_indexes.size(); i++)
	{
		//Copy ith fitness of each individual
		double fit_ind1;
		double fit_ind2;

		short fit_ind_ix = fit_indexes[i]-1;

		if (TGM == false)
		{
			fit_ind1 = ind1.Fitness_Show(fit_ind_ix);		//ith fitness of 1st individual
			fit_ind2 = ind2.Fitness_Show(fit_ind_ix);		//ith fitness of 2nd individual
		}
		else
		{
			fit_ind1 = ind1.TGM_fitness[0][fit_ind_ix];		//ith fitness of 1st individual
			fit_ind2 = ind2.TGM_fitness[0][fit_ind_ix];		//ith fitness of 2nd individual
		}
													//Check which fitness is greater
		if (fit_ind1 < fit_ind2)
			flag1 = 1;
		else if (fit_ind1 > fit_ind2)
			flag2 = 1;
	}
	//Check which individual dominate
	if (flag1 == 1 && flag2 == 0)
		return (1);
	else if (flag1 == 0 && flag2 == 1)
		return (-1);
	else
		return (0);
}

/*
*Check Dominance - without constrained check*
@param ind1 - first individual
@param ind2 - 2nd individual
Return:
1 - ind1 dominates
-1 - ind2 dominates
0 - both nondominated
*/
int Dominance_Check_NCons(individual &ind1, individual &ind2, std::vector<short> &fit_indexes)
{
	int flag1, flag2;		//flags
	flag1 = flag2 = 0;
	int nobj = ind1.Fitness_Show().size();	//number of objectives

	for (int i = 0; i < fit_indexes.size(); i++)
	{
		//Copy ith fitness of each individual
		double fit_ind1;
		double fit_ind2;

		short fit_ind_ix = fit_indexes[i] - 1;

		if (TGM == false)
		{
			fit_ind1 = ind1.Fitness_Show(fit_ind_ix);		//ith fitness of 1st individual
			fit_ind2 = ind2.Fitness_Show(fit_ind_ix);		//ith fitness of 2nd individual
		}
		else
		{
			fit_ind1 = ind1.TGM_fitness[0][fit_ind_ix];		//ith fitness of 1st individual
			fit_ind2 = ind2.TGM_fitness[0][fit_ind_ix];		//ith fitness of 2nd individual
		}
		//Check which fitness is greater
		if (fit_ind1 < fit_ind2)
			flag1 = 1;
		else if (fit_ind1 > fit_ind2)
			flag2 = 1;
	}
	//Check which individual dominate
	if (flag1 == 1 && flag2 == 0)
		return (1);
	else if (flag1 == 0 && flag2 == 1)
		return (-1);
	else
		return (0);
}


/*
*Fitness calculation by chosen Scalarizing Function*
@param fit - fitness vector
@param namda - weight vector
@param iGen - 
*/
double Fitness_Function(std::vector <double> &fit, std::vector <double> &namda, std::vector<short> & fit_indexes, int iGen )
{
	extern std::string MODE;											//current hybrid mod


	// Chebycheff Scalarizing Function
	double max_fun;

	if (MODE == "MOEAD" || MODE == "BCE" || MODE == "MOEADM2M")
	{
		max_fun = Fitness_Function_TCH(fit, namda, fit_indexes);
	}
	else if (MODE == "MOEADPSF")
	{
		max_fun = Fitness_Function_PSF(fit, namda, fit_indexes,  iGen);
	}
	else if (MODE == "MOEADMSF")
	{
		max_fun = Fitness_Function_MSF(fit, namda, fit_indexes,  iGen);
	}
	else
		abort();

	return max_fun;
}

/*
*Fitness calculation by Chebycheff Scalarizing Function*
@param fit - fitness vector
@param namda - weight vector
*/
double Fitness_Function_TCH(std::vector <double> &fit, std::vector <double> &namda, std::vector<short> & fit_indexes)
{
	extern std::vector<double> idealpoint;			//Storage of ideal point

	// Chebycheff Scalarizing Function
	double max_fun = -1.0e+30;

	int nobj = fit_indexes.size();

	std::vector<double> idealpoint_temp(nobj,0.);

	if (MLSGA_norm_obj == false)
		idealpoint_temp = idealpoint;


	for (int n = 0; n<nobj; n++)
	{
		short fit_ix = fit_indexes[n] - 1;
		double diff = fabs(fit[fit_ix] - idealpoint_temp[fit_ix]);
		double feval;
		if (namda[fit_ix] == 0)
			feval = 0.0001*diff;
		else
			feval = diff*namda[fit_ix];
		if (feval>max_fun)
			max_fun = feval;

	}

	return max_fun;
}


/*
*Fitness calculation by PSF Scalarizing Function*
@param fit - fitness vector
@param namda - weight vector
*/
double Fitness_Function_PSF(std::vector <double> &fit, std::vector <double> &namda, std::vector<short> & fit_indexes,  int iGen)
{
	extern std::vector<double> idealpoint;			//Storage of ideal point
	extern int nfes;						//number of function evaluations

	double max_fun = -1.0e+30, min_fun = 1.0e+30;
	double au = 1, wmax = -1.0e+30, wmin = 1.0e+30;
	double eps = 1.0e-6;


	double penalty;
	double gen_val;

	if (T_con == "nfes")
	{
		gen_val = (float)nfes / max_iterations;

	}
	else
	{
		gen_val = (float)iGen / max_generations;
	}


	double func_val = 0; //output

	//Get the number of objectives
	int nobj = fit_indexes.size();
	std::vector<double> idealpoint_temp(nobj, 0.);

	if (MLSGA_norm_obj == false)
		idealpoint_temp = idealpoint;


	//if (nobj > 2) eps = 1.0e-3;
	for (int n = 0; n<nobj; n++)
	{
		short fit_ix = fit_indexes[n] - 1;
		double diff = fabs(fit[fit_ix] - idealpoint_temp[fit_ix]);
		double feval;
		if (namda[fit_ix] == 0.)
		{
			feval = (diff) / eps;
		}
		else
			feval = (diff) / namda[fit_ix];
		if (feval>max_fun)
			max_fun = feval;
		if (feval<min_fun) 
			min_fun = feval;
		if (namda[fit_ix]>wmax) 
			wmax = namda[fit_ix];
		if (namda[fit_ix] < wmin) 
			wmin = namda[fit_ix];
	}

	// normalize the weight vector (line segment)
	double nd = norm_vector(namda);

	std::vector<double> realA(nobj);
	std::vector<double> temp_namda;
	// difference beween current point and reference point
	for (int n = 0; n < nobj; n++)
	{
		short fit_ix = fit_indexes[n] - 1;
		realA[n] = (fit[fit_ix] - idealpoint_temp[fit_ix]);
		temp_namda.push_back(namda[fit_ix]);
	}


	double A1 = norm_vector(realA);
	double AB = fabs(innerproduct(realA, temp_namda));
	double d2 = pow(A1, 2.0) - pow(AB, 2.0) / pow(nd, 2.0);
	d2 = pow(abs(d2), 0.5);

	au = nobj * wmin;
	penalty = (1.0 - gen_val);
	func_val = max_fun + 10.0*au*penalty*d2;
	//func_val = max_fun;
	return func_val;
}
/*
*Fitness calculation by MSF Scalarizing Function*
@param fit - fitness vector
@param namda - weight vector
*/
double Fitness_Function_MSF(std::vector <double> &fit, std::vector <double> &namda, std::vector<short> & fit_indexes,  int iGen)
{
	//Get the number of objectives
	int nobj = fit_indexes.size();
	extern int nfes;						//number of function evaluations


	double func_val = 0; //output

	double penalty;
	double alpha = 1.;

	double gen_val;

	if (T_con == "nfes")
	{
		gen_val = (float)nfes / max_iterations;

	}
	else
	{
		gen_val = (float)iGen / max_generations;
	}

	double max_fun = -1.0e+30, min_fun = 1.0e+30;
	double au = 1, wmax = -1.0e+30, wmin = 1.0e+30;
	double eps = 1.0e-6;
	if (nobj > 2) eps = 1.0e-3;
	for (int n = 0; n<nobj; n++)
	{
		short fit_ix = fit_indexes[n] - 1;
		double diff = fabs(fit[fit_ix]);
		double feval;
		if (namda[fit_ix] == 0)
		{
			feval = (diff) / eps;
		}
		else
			feval = (diff) / namda[fit_ix];
		if (feval>max_fun) max_fun = feval;
		if (feval<min_fun) min_fun = feval;
		if (namda[fit_ix]>wmax) wmax = namda[fit_ix];
		if (namda[fit_ix] < wmin) wmin = namda[fit_ix];
	}
	au = nobj * wmin;

	penalty = (1.0 - gen_val);
	func_val = pow(max_fun / min_fun, alpha*au*penalty)*max_fun;


	return func_val;
}



std::vector<std::vector<double>> Uniform_Weights_Generate(int pop_size, short nobj, short type)
{
	std::vector<std::vector<double>> val_w;

	if (M2M_weight_method == 1)
		val_w = reference_point(pop_size, nobj);
	else if (M2M_weight_method == 2)
		val_w = recursion_point(pop_size + 1, nobj);
	else
		abort();

	if (type == 1)
	{
		for (int i = 0; i < val_w.size(); i++)
		{
			//get the value
			double temp_val = 0;
			for (int j = 0; j < val_w[i].size(); j++)
				temp_val += pow(val_w[i][j], 2);
			temp_val = sqrt(temp_val);

			for (int j = 0; j < val_w[i].size(); j++)
				val_w[i][j] /= temp_val;
		}
	}
	else if (type == 2)
	{
		for (int i = 0; i < val_w.size(); i++)
		{
			//get the value
			double temp_val = 0;
			for (int j = 0; j < val_w[i].size(); j++)
				temp_val += sqrt(val_w[i][j]);
			temp_val = temp_val/(pow(temp_val,3));

			for (int j = 0; j < val_w[i].size(); j++)
				val_w[i][j] *= temp_val;
		}
	}

	return val_w;
}

std::vector<std::vector<double>>  reference_point(int pop_size, short nobj)
{
	//Method not initialised
	abort();
						//number of objectives
	std::vector<std::vector<double>> out(1, std::vector<double>(nobj, 0));
	return out;
}

std::vector<std::vector<double>> recursion_point(int pop_size, short nobj)
{
	//number of objectives
	int temp_nobj = nobj - 1;
	int layer = 0;
	int wnum = 0;

	while (wnum < pop_size)
	{
		layer++;
		wnum = Get_Weight_Num(layer, temp_nobj, 0);
	}
	layer--;

	wnum = Get_Weight_Num(layer, temp_nobj, 0);
	std::vector<std::vector<double>> out(wnum, std::vector<double>(nobj, 0));

	double step = 1. / (layer - 1.);

	Get_Weight(layer, temp_nobj, step, out);

	for (int i = 0; i < out.size(); i++)
	{
		double temp_sum = 0;
		for (int j = 0; j < temp_nobj; j++)
			temp_sum += out[i][j];
		if (temp_sum >= 1)
			temp_sum = 1;

		out[i][temp_nobj] = 1. - temp_sum;
	}
	return out;
}

int Get_Weight_Num(int layer, int od, int wnum)
{
	int out = wnum;
	if (od == 1)
		out = wnum + layer;
	else
	{
		for (int i = 1; i <= layer; i++)
		{
			out = Get_Weight_Num(i, od - 1, out);
		}
	}
	return out;
}

void Get_Weight(int layer, int od, double step, std::vector<std::vector<double>> &out)
{
	if (out.empty())
		return;

	for (int i = 1; i <= layer; i++)
	{
		int n1 = Get_Weight_Num(i - 1, od, 0);
		int n2 = Get_Weight_Num(i, od, 0);

		for (int j = n1; j < n2; j++)
		{
			out[j][od - 1] = step*(layer - i);
		}
		if (od > 1)
		{
			std::vector<std::vector<double>> temp(out.begin() + n1, out.begin() + n2);
			Get_Weight(i, od - 1, step, temp);

			for (int j = n1, k = 0; j < n2; j++, k++)
			{
				out[j] = temp[k];
			}
		}
	}
}


//Generation of uniformly spread points
std::vector<std::vector<double>> Points_Generate(int size, short n_obj)
{
	std::vector<std::vector<double>> out;
	int temp_size = pow(size, 1. / n_obj);
	double step = 1. / (temp_size - 1.);

	std::vector<double> short_out;
	Points_Generate_Recursion(step, n_obj, out, short_out);

	return out;
}

void Points_Generate_Recursion(double step, short n_obj, std::vector<std::vector<double>>  & out, std::vector<double> short_out)
{
	short_out.push_back(0);
	for (double i = 0; i < 1 + step; i += step)
	{
		if (i > 1.)
			i = 1;
		short_out.back() = i;
		if (short_out.size() == n_obj)
			out.push_back(short_out);
		else
			Points_Generate_Recursion(step, n_obj, out, short_out);
	}
}
extern std::vector<double> idealpoint;			//Storage of ideal point
extern std::vector<double> nadirpoint;

std::vector<double> Normalize_Objective(std::vector<double> & fitness)
{

	if (idealpoint.empty() || nadirpoint.empty())
	{
		return idealpoint;
	}

	std::vector<double> norm_obj;

	int nobj = fitness.size();

	for (int n = 0; n < nobj; n++)
		norm_obj.push_back((fitness[n] - idealpoint[n]) / (nadirpoint[n] - idealpoint[n]));

	return norm_obj;
}

void Clear_Ideal_and_Nadir(short nobj)
{
	nadirpoint.clear();
	idealpoint.clear();
	idealpoint = std::vector<double>(nobj, 1.0e+30);
}

void Create_Nadir(short nobj)
{
	nadirpoint = std::vector<double>(nobj, -1.0e+30);
}

void Update_Nadirpoint(std::vector<individual> & pop, short flag)
{
	int pops = pop.size();
	int nobj = pop[0].Fitness_Show().size();


	for (int n = 0; n < nobj; n++)
	{
		if (flag == 1)
			nadirpoint[n] = -INFINITY;
		for (int i = 0; i < pops; i++)
		{
			if (pop[i].Fitness_Show(n)>nadirpoint[n])
			{
				nadirpoint[n] = pop[i].Fitness_Show(n);
			}
		}
	}
}


/*
*Update the reference - idealpoint vector*
@param ind - individual used for update
*/
void Idealpoint_Update(std::vector<double> &ind)
{
	int nobj = ind.size();

	for (int n = 0; n<nobj; n++)
	{
		if (ind[n]<idealpoint[n])
		{
			idealpoint[n] = ind[n];
		}
	}
}


double norm_vector(std::vector<double> &x)
{
	double sum = 0;
	for (int i = 0; i<x.size(); i++)
		sum = sum + x[i] * x[i];
	return sqrt(sum);
}

double innerproduct(std::vector<double> &vec1, std::vector<double> &vec2)
{
	double sum = 0;
	for (int i = 0; i<vec1.size(); i++)
		sum += vec1[i] * vec2[i];
	return sum;
}


bool Termination_Check(bool param)
{
	extern int nfes;
	extern time_t run_t;										//Time of run
	extern float run_time;
	extern int pop_size;

	if (T_con == "ngen")
		return false;
	else if (T_con == "ntime")
	{
		run_time = ((clock() - run_t) * 1000 / CLOCKS_PER_SEC);
		run_time /= 1000.;

		if (run_time >= max_time)
			return true;
	}
	else if (T_con == "nfes")
	{
		if (nfes >= max_iterations)
			return true;
		else if (param)
			if (nfes % (pop_size *T_dyn) == 0)
				return true;
	}
	else
		abort();

	return false;
}


