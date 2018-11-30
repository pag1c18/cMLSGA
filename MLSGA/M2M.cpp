#include "M2M.h"
#include "Support_Functions.h"
#include "MOEAD.h"
#include "MTS.h"
#include "NSGAII.h"
#include <ctime>

std::vector<std::vector<int>> subclass_index;
std::vector<std::vector<std::vector<double>>> all_col_class_centres;
extern std::vector<double> idealpoint;			//Storage of ideal point
extern std::vector<double> nadirpoint;
extern time_t elit_t;					//Time of 
extern int nfes;
extern short n_col;

namespace M2M
{
	std::vector<std::vector<double>> GW(std::vector<std::vector<double>> &val_w, int class_num);
	std::vector<std::vector<int>> Group_Index(std::vector<std::vector<double>> center, int group_num, std::vector<std::vector<double>> val_w);
	std::vector<individual> CMOp(collective & col, std::vector<individual> & old_col_indi);
	void Extract_Pop(collective & new_col, std::vector<individual> & old_col);
	std::vector<individual> Select_MM(std::vector<individual> & sel_pop, std::vector<std::vector<double>> & weights);
}

std::vector<individual> M2M::M2M_Calc(collective & col, int iGen)
{
	if (col.Was_Erased())
	{		
		//create the neighbourhood
		M2M::M2M_Population_Init(col.FCode_Show()[0].Objs(), col);
		col.Clear_Erased();
	}
	time_t elite_t_temp = clock();				//starting time of whole operation
	//get the index of the collective
	short col_index = col.Index_Show() - 1;
	int col_size = col.Size_Show();
	std::vector<individual> PF_out;

	std::vector<individual> old_col_indi = col.Indiv_Show();
	col.Erase();

	for (int i = 0; i < M2M_class_num; i++)
	{
		//Create the separate collective for single subproblem
		//copy the individuals
		std::vector<individual> temp_col_indi(old_col_indi.begin()+subclass_index[col_index][i], old_col_indi.begin() + subclass_index[col_index][i+1]);
		collective temp_col(temp_col_indi,col);

		//Do evolutionary routine
		std::vector<individual> temp_PF;

		if (M2M_method == "M2M")
		{
			temp_col.Clear_Erased();
			temp_PF = CMOp(temp_col, old_col_indi);
			temp_col.save();
		}
		else if (M2M_method == "MOEAD")
		{ 
			temp_PF = MOEAD::MOEAD_Calc(temp_col, iGen);
		}
		else if (M2M_method == "NSGAII")
		{
			NSGAII::NSGAII_Calc(temp_col, iGen);
			temp_PF = temp_col.Indiv_Show();
		}
		else if (M2M_method == "MTS")
		{
			MTS::MTS_Calc(temp_col);
		}
		else
			abort();




		//Push the temp_PF to output vector
		PF_out.insert(PF_out.end(), temp_PF.begin(), temp_PF.end());

		//replace the individuals in the original collective
		col.Add(temp_col.Indiv_Show());
	}
	col.Clear_Erased();
	//if (col.Size_Show() != col_size)
	//	system("pause");


	Extract_Pop(col, old_col_indi);


	col.Fitness_Calc();
	//calculate time
	elit_t += clock() - elite_t_temp;

	return PF_out;
}

void M2M::M2M_Init(short n_obj, population & pop)
{
	time_t elite_t_temp = clock();				//starting time of whole operation

	//remove the old subclasses index
	subclass_index.clear();

	
	//create new storage
	subclass_index = std::vector<std::vector<int>>(n_col, std::vector<int>(M2M_class_num, 0));

	all_col_class_centres = std::vector<std::vector<std::vector<double>>>(n_col, std::vector<std::vector<double>>());

	//calculate time
	elit_t += clock() - elite_t_temp;

	//use MOEA/D routine
	MOEAD::MOEAD_Init(n_obj, pop);
}

void M2M::M2M_Population_Init(short n_obj, collective & col)
{
	//use MOEA/D routine
	MOEAD::Population_Init(n_obj, col);

	time_t elite_t_temp = clock();				//starting time of whole operation

	int col_size = col.Size_Show();

	std::vector<std::vector<std::vector<double>>> col_class_weights;
	std::vector<std::vector<double>> col_class_centers;

	//get parameters for the subclasses
	M2M_Get_Params(col_class_weights, col_class_centers, col_size, n_obj);
	//assign the class center to global
	all_col_class_centres[col.Index_Show()-1] = col_class_centers;


	std::vector<std::vector<double>> val_w;		//storage of normalized objectives

	//Get the objective values
	for (int i = 0; i < col_size; i++)
	{
	
		val_w.push_back(col.Indiv_Show(i).Fitness_Show());
		for (int j = 0; j < n_obj; j++)
		{
			val_w[i][j] -= idealpoint[j];
		}
	}
	//Copy the individuals for ease of operations
	std::vector<individual> Col_Indi = col.Indiv_Show();


	//Calculate new indexes
	std::vector<std::vector<int>> indexes = Group_Index(col_class_centers, M2M_class_num, val_w);

	//Assign individuals to each subclass
	//Create the storage of individuals for each subclass
	std::vector<std::vector<individual>> Subclass_Individuals(M2M_class_num,std::vector<individual>());

	for (int i = 0; i < M2M_class_num; i++)
	{
		int num_prob = indexes[i].size();
		int class_size = col_class_weights[i].size();

		if (num_prob <= class_size)
		{
			while (num_prob < class_size)
			{
				indexes[i].push_back(Random_I(0, col_size - 1));
				num_prob++;
			}
			for (int j = 0; j < num_prob; j++)
			{
				Subclass_Individuals[i].push_back(Col_Indi[indexes[i][j]]);
				Subclass_Individuals[i][j].namda = col_class_weights[i][j];
			}
		}
		else
		{
			for (int j = 0; j < num_prob; j++)
			{
				Subclass_Individuals[i].push_back(Col_Indi[indexes[i][j]]);
			}
			Subclass_Individuals[i] = Select_MM(Subclass_Individuals[i], col_class_weights[i]);
		}
		if (Subclass_Individuals[i].size() != col_class_weights[i].size())
			abort();
		
	}


	//Copy the individuals to the colelctive and save indexes of each subpop
	int indiv_index = 0;
	int col_index = col.Index_Show() - 1;
	subclass_index[col_index].clear();
	subclass_index[col_index].push_back(0);
	for (int i = 0; i < M2M_class_num; i++)
	{
		for (int j = 0; j < Subclass_Individuals[i].size(); j++)
		{
			col.Indiv_Set(indiv_index) = Subclass_Individuals[i][j];
			indiv_index++;
		}
		subclass_index[col_index].push_back(indiv_index);
	}

	if (col.Size_Show() != col_size)
		abort();

	col.Clear_Erased();
	//calculate time
	elit_t += clock() - elite_t_temp;
}

void M2M::M2M_Get_Params(std::vector<std::vector<std::vector<double>>> &weights, std::vector<std::vector<double>> &center, int pop_size, short n_obj)
{
	weights.clear();
	center.clear();
	//Get weights for the whole population
	std::vector<std::vector<double>> val_w = Uniform_Weights_Generate(pop_size, n_obj, 1);


	//update the pop_size
	pop_size = val_w.size();

	//Create subproblems
	center = GW(val_w, M2M_class_num);

	std::vector<std::vector<int>> indexes = Group_Index(center, M2M_class_num, val_w);

	for (int i = 0; i < M2M_class_num; i++)
	{
		std::vector<std::vector<double>> temp_pop_weights;
		for (int j = 0; j < indexes[i].size(); j++)
		{
			temp_pop_weights.push_back(val_w[indexes[i][j]]);
		}
		weights.push_back(temp_pop_weights);
	}

}




//Get the values for centres
std::vector<std::vector<double>> M2M::GW(std::vector<std::vector<double>> &val_w, int class_num)
{
	std::vector<std::vector<double>> out;			//output

	int size = val_w.size();

	//create the aggregate vector
	std::vector<std::vector<double>> dis(size, std::vector<double>(size, 0));



	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			double value = 0;
			for (int k = 0; k < val_w[0].size(); k++)
			{
				value += val_w[i][k] * val_w[j][k];
			}
			dis[i][j] = value;
		}
	}
	double loc = Random_I(0, size - 1);
	out.push_back(val_w[loc]);
	std::vector<double> pdis = dis[loc];
	for (int i = 1; i < class_num; i++)
	{
		//get the index of min value of pdis
		double temp_min_val = 1000000;
		for (int j = 0; j < size; j++)
		{
			if (pdis[j] < temp_min_val)
			{
				temp_min_val = pdis[j];
				loc = j;
			}
		}

		//push the weights
		out.push_back(val_w[loc]);

		//update the pdis
		for (int j = 0; j < size; j++)
		{
			if (pdis[j] < dis[loc][j])
				pdis[j] = dis[loc][j];
		}

	}

	return out;
}

std::vector<std::vector<int>> M2M::Group_Index(std::vector<std::vector<double>> center, int group_num, std::vector<std::vector<double>> val_w)
{
	std::vector<std::vector<int>> out(group_num, std::vector<int>()); //output vector

	int size = val_w.size();

	std::vector<std::vector<double>> dis;

	for (int i = 0; i < size; i++)
	{
		std::vector<double> temp_vect;
		for (int j = 0; j < group_num; j++)
		{
			double value = 0;
			for (int k = 0; k < val_w[0].size(); k++)
			{
				value += val_w[i][k] * center[j][k];
			}
			temp_vect.push_back(value);
		}
		dis.push_back(temp_vect);
	}
	std::vector<int> minindex;
	//get the min indexes
	for (int i = 0; i < size; i++)
	{
		int temp_min_index = 0;
		for (int j = 1; j < group_num; j++)
		{
			if (dis[i][j] > dis[i][temp_min_index])
				temp_min_index = j;
		}
		minindex.push_back(temp_min_index);
	}

	for (int i = 0; i < size; i++)
	{
		out[minindex[i]].push_back(i);
	}

	return out;
}

std::vector<individual> M2M::CMOp(collective & col, std::vector<individual> & old_col_indi)
{
	//get the size
	int col_size = col.Size_Show();
	//Copy the individuals
	std::vector<individual> col_indi = col.Indiv_Show();

	std::vector<individual> out;

	for (int i = 0; i < col_size; i++)
	{
		std::vector<individual> Parents;
		Parents.push_back(col_indi[i]);

		if (Random_F() < M2M_select_pro)
		{
			int loc = Random_I(0, col_size - 1);
			Parents.push_back(col_indi[loc]);
		}
		else
		{
			int loc = Random_I(0, old_col_indi.size() - 1);
			Parents.push_back(old_col_indi[loc]);
		}
		/*std::vector<individual> offsprings = col.Crossover_NSGAII(Parents);

		for (int i = 0; i < 2; i++)
		{
			col.Mutation(offsprings[i]);
			col.population::Fitness_Calc(offsprings[i]);
			nfes++;
		}

		double f1, f2, f_parent;
		if (TGM == false)
		{
			f1 = Fitness_Function(offsprings[0].Fitness_Show(), Parents[0].namda, 0);
			f2 = Fitness_Function(offsprings[1].Fitness_Show(), Parents[0].namda, 0);
			f_parent = Fitness_Function(Parents[0].Fitness_Show(), Parents[0].namda, 0);
		}
		else
		{
			f1 = Fitness_Function(offsprings[0].TGM_fitness[0], Parents[0].namda, 0);
			f2 = Fitness_Function(offsprings[1].TGM_fitness[0], Parents[0].namda, 0);
			f_parent = Fitness_Function(Parents[0].Fitness_Show(), Parents[0].namda, 0);
		}
		short ix;	//index of better individual
		if (f1 < f2)
		{
			if (f1 < f_parent)
				ix = 0;
			else
				ix = 2;
		}
		else
		{
			if (f2 < f_parent)
				ix = 1;
			else
				ix = 2;
		}
		individual offspring = Parents[0];

		if (ix != 2)
			offspring = offsprings[ix];*/

		individual offspring = col.Crossover_NSGAII(Parents)[0];
		col.Mutation(offspring);

		col.population::Fitness_Calc(offspring);
		nfes++;

		offspring.namda = Parents[0].namda;
		offspring.table = Parents[0].table;
		offspring.saved_fitness = Parents[0].saved_fitness;

		if (TGM)
			offspring.TGM_fitness = Parents[0].TGM_fitness;

		

		out.push_back(offspring);
	}

	col = collective(out, col);

	return out;
}

void M2M::Extract_Pop(collective & new_col, std::vector<individual> & old_col)
{
	//Copy the individuals for ease of operations
	std::vector<individual> Col_Indi = new_col.Indiv_Show();
	Col_Indi.insert(Col_Indi.begin(), old_col.begin(), old_col.end());

	int num_nod = Col_Indi.size();
	int nobj = new_col.FCode_Show()[0].Objs();


	//get the index of the collective
	short col_index = new_col.Index_Show() - 1;

		

	std::vector < std::vector < double>> noval(num_nod, std::vector < double>(nobj, 0.));
	
	for (int i = 0; i < num_nod; i++)
	{
		//get the objectives of current individual
		std::vector<double> ind_fit = Col_Indi[i].Fitness_Show();
		for (int j = 0; j < nobj; j++)
		{
			noval[i][j] = ind_fit[j] - idealpoint[j];
		}
	}
	
	//Get the class centres for current collective
	std::vector<std::vector<double>> col_class_centers = all_col_class_centres[col_index];
	//Get the weight for current collective for each subclass
	std::vector<std::vector<std::vector<double>>> col_class_weights(M2M_class_num, std::vector<std::vector<double>>());
	for (int i = 0; i < M2M_class_num; i++)
	{
		int pop_size = subclass_index[col_index][i + 1] - subclass_index[col_index][i];
		std::vector<std::vector<double>> subclass_weights(pop_size,std::vector<double>());
		for (int j = 0; j < pop_size; j++)
		{
			subclass_weights[j] = Col_Indi[subclass_index[col_index][i] + j].namda;
		}
		col_class_weights[i] = subclass_weights;
	}


	//Calculate new indexes
	std::vector<std::vector<int>> indexes = Group_Index(col_class_centers, M2M_class_num, noval);

	//Assign individuals to each subclass
	//Create the storage of individuals for each subclass
	std::vector<std::vector<individual>> Subclass_Individuals(M2M_class_num, std::vector<individual>());

	for (int i = 0; i < M2M_class_num; i++)
	{
		int num_prob = indexes[i].size();
		int pop_size = subclass_index[col_index][i+1] - subclass_index[col_index][i];

		if (Random() < M2M_global_pro)
		{
			if (num_prob <= pop_size)
			{
				while (num_prob < pop_size)
				{
					indexes[i].push_back(Random_I(0, num_nod - 1));
					num_prob++;
				}
				for (int j = 0; j < num_prob; j++)
				{
					Subclass_Individuals[i].push_back(Col_Indi[indexes[i][j]]);
					Subclass_Individuals[i][j].namda = col_class_weights[i][j];
				}
			}
			else
			{
				for (int j = 0; j < num_prob; j++)
				{
					Subclass_Individuals[i].push_back(Col_Indi[indexes[i][j]]);
				}
				Subclass_Individuals[i] = Select_MM(Subclass_Individuals[i], col_class_weights[i]);
			}
		}
		else
		{
			Subclass_Individuals[i] = Select_MM(Col_Indi, col_class_weights[i]);
		}
		if (Subclass_Individuals[i].size() != pop_size)
			abort();
	}


	//Copy the individuals to the colelctive and save indexes of each subpop
	int indiv_index = 0;
	for (int i = 0; i < M2M_class_num; i++)
	{
		for (int j = 0; j < Subclass_Individuals[i].size(); j++)
		{
			new_col.Indiv_Set(indiv_index) = Subclass_Individuals[i][j];
			indiv_index++;
		}
	}

	if (new_col.Size_Show() != old_col.size() && DEBUG == true)
		abort();
}

std::vector<individual> M2M::Select_MM(std::vector<individual> & sel_pop, std::vector<std::vector<double>> & weights)
{


	std::vector<individual> out;		//output

	int num_nod = sel_pop.size();
	int nobj = weights[0].size();
	int out_size = weights.size();

	std::vector < std::vector < double>> noval(num_nod, std::vector < double>(nobj, 0.));
	
	
	for (int i = 0; i < num_nod; i++)
	{
		//get the objectives of current individual
		std::vector<double> ind_fit = sel_pop[i].Fitness_Show();
		for (int j = 0; j < nobj; j++)
		{
			noval[i][j] = ind_fit[j] - idealpoint[j];
		}
	}
	

	for (int i = 0; i < out_size; i++)
	{
		int index_min;
		double temp_min = 1e30;
		for (int j = 0; j < num_nod; j++)
		{
			double temp_max = -1e30;
			for (int k = 0; k < nobj; k++)
			{
				if (noval[j][k] * weights[i][k] > temp_max)
					temp_max = noval[j][k] * weights[i][k];
			}
			if (temp_max < temp_min)
			{
				index_min = j;
				temp_min = temp_max;
			}
		}

		out.push_back(sel_pop[index_min]);
		out[i].namda = weights[i];
		noval.erase(noval.begin() + index_min);
		num_nod--;
	}
	if (out.size() != weights.size())
		abort();
	return out;
}