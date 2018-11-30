#include "SVM.h"
#include <numeric>
#include "Imported\svm_i.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

/*"Usage: svm-train [options] training_set_file [model_file]\n"
"options:\n"
"-s svm_type : set type of SVM (default 0)\n"
"	0 -- C-SVC		(multi-class classification)\n"
"	1 -- nu-SVC		(multi-class classification)\n"
"	2 -- one-class SVM\n"
"	3 -- epsilon-SVR	(regression)\n"
"	4 -- nu-SVR		(regression)\n"
"-t kernel_type : set type of kernel function (default 2)\n"
"	0 -- linear: u'*v\n"
"	1 -- polynomial: (gamma*u'*v + coef0)^degree\n"
"	2 -- radial basis function: exp(-gamma*|u-v|^2)\n"
"	3 -- sigmoid: tanh(gamma*u'*v + coef0)\n"
"	4 -- precomputed kernel (kernel values in training_set_file)\n"
"-d degree : set degree in kernel function (default 3)\n"
"-g gamma : set gamma in kernel function (default 1/num_features)\n"
"-r coef0 : set coef0 in kernel function (default 0)\n"
"-c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)\n"
"-n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)\n"
"-p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)\n"
"-m cachesize : set cache memory size in MB (default 100)\n"
"-e epsilon : set tolerance of termination criterion (default 0.001)\n"
"-h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)\n"
"-b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)\n"
"-wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)\n"
"-v n: n-fold cross validation mode\n"
"-q : quiet mode (no outputs)\n"*/

static struct svm_parameter param;			//SVM parameters structure
static struct svm_problem prob, prob2;		//SVM problems structures
static struct svm_model *model;				//SVM model structure
static struct svm_node *x_space, *x_space2;	//SVM nodes for data storage
static int error;							//error of the SVM generation
static int correct;						//number of correct solutions

extern short n_col;

/*
*SVM function for given population - generate a vector of labels for separation*
@param pop population for which SVM will calculate labels
*/
std::vector<short> SVM(const population & pop)
{
	std::vector<short> real_label;				//vector for the real labels storage - output
	

	//set up SVM parameters
	SVM_PARAM();

	//get population size
	int psize = pop.Size_Show();				//source population size

	//generate training population for SVM_Train
	population pop_init{ pop, 0 };				//training population

	//Add individuals to the training population
	for (int p = 0; p < psize; p++)
		pop_init.Add();
	
	//Train SVM and generate SVM model
	SVM_TRAIN(pop_init);

	//Predict SVM and generate real labels
	real_label = SVM_PREDICT(pop);

	

	//clear the SVM data
	{
		svm_destroy_param(&param);
		free(prob.y);
		free(prob.x);
		prob.l = 0;
		free(x_space);
		free(prob2.y);
		free(prob2.x);
		prob2.l = 0;
		free(x_space2);
		svm_free_and_destroy_model(&model);
	}

	return real_label;
}

/*
*SVM training and model generation*
@param pop population which will be trained
*/
void SVM_TRAIN(const population & pop)
{
	//generate problem and data nodes for the given population
	SVM_PROB(pop);

	//make model for the problem
	model = svm_train(&prob, &param);
	
}

/*
*SVM Real label creation*
@param pop population for which label will be calculated
*/
std::vector<short> SVM_PREDICT(const population & pop)
{

	std::vector<short> real_label;		//vector for the real labels storage - output
	
	//generate problem and data nodes for the given population
	SVM_PROB2(pop);


	double predict_label;				//label predicted by the SVM
	double target_label;				//label calculated (based on distribution of individuals variables

	//clear error and correct values
	error = correct = 0;

	//generate label for each individual
	for (int i = 0; i < prob2.l; i++)
	{
		//assign problem label to targer_label
		target_label = prob2.y[i];
	
		//generate predicted label by SVM
		predict_label = svm_predict(model, prob2.x[i]);

		//push label to the vector
		real_label.push_back(predict_label);

		//check if labels are the same
		if (predict_label == target_label)
			++correct;
		//else calculate error
		else
			error += (predict_label - target_label)*(predict_label - target_label);
	}
	
	//print the statistical data
	{
		std::cout << "\nCorrect: " << correct << "/" << prob2.l << std::endl;
		std::cout << (float)correct / (float)prob2.l * 100.f << "%\n";
		std::cout << "Error: " << error << std::endl << std::endl;
	}

	return real_label;

}

/**SVM Setting up  SVM parameters**/
void SVM_PARAM()					
{
	//parameters are described in the header file
	param.svm_type = C_SVC;
	param.kernel_type = 0;
	param.degree = 3;
	param.gamma = 0;	// 1/num_features
	param.coef0 = 0;
	param.nu = 0.5;
	param.cache_size = 100;
	param.C = 1;
	param.eps = 1e-3;
	param.p = 0.1;
	param.shrinking = 1;
	param.probability = 0;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;
}

/*
*Generation of problem for SVM_TRAIN*
@param pop population for which problem will be separated
*/
void SVM_PROB(const population & pop)		//Reading population
{
	//get population size
	int psize = pop.Size_Show();			//population size

	std::vector<int> train_label;			//vector for the training-labels storage

	//generate labels for each individual
	train_label = SVM_LABEL(pop);

	//get max index of the individual coede
	int max_index = pop.Indiv_Show(1).Code_Show().size();		//max index of the individual coede

	//calculate the number of elements in the data node
	size_t elements = psize * (max_index + 1);				//number of elements in the data node

	//assign size of the population the the prob structure
	prob.l = psize;											//amount of labels

	//allocate space for each structure
	prob.y = Malloc(double, prob.l);						//value of each label				
	prob.x = Malloc(struct svm_node *, prob.l);				//data of the individual
	x_space = Malloc(struct svm_node, elements);			//data of all individuals


	int j = 0;

	for (int i = 0; i < psize; i++)
	{
		//copy the code of individual to temporary vector
		std::vector<double> code = pop.Indiv_Show(i).Code_Show();		//temporary vector for the individual's code storage
		
		//assign pointer to the start of the memory with data of the individual
		prob.x[i] = &x_space[j];

		//assign label of the individual to the problem
		prob.y[i] = train_label[i];

		//assign data of individual to the storage node
		for (int k = 0; k < max_index; k++)
		{
			x_space[j].index = k + 1;
			x_space[j].value = code[k];
			++j;
		}
		x_space[j++].index = -1;
	}

	//calculate new gamma parameter - only if defined
	if (param.gamma == 0 && max_index > 0)
		param.gamma = 1.0 / max_index;

}

/*
*Generation of problem for SVM_PREDICT*
@param pop population for which problem will be separated
*/
void SVM_PROB2(const population & pop)		//Reading population
{
	//get population size
	int psize = pop.Size_Show();			//population size

	std::vector<int> train_label;			//vector for the training-labels storage

	//generate labels for each individual
	train_label = SVM_LABEL(pop);

	//get max index of the individual coede
	int max_index = pop.Indiv_Show(1).Code_Show().size();		//max index of the individual coede

	//calculate the number of elements in the data node
	size_t elements = psize * (max_index + 1);				//number of elements in the data node

	//assign size of the population the the prob structure
	prob2.l = psize;											//amount of labels

	//allocate space for each structure
	prob2.y = Malloc(double, prob2.l);						//value of each label				
	prob2.x = Malloc(struct svm_node *, prob2.l);				//data of the individual
	x_space2 = Malloc(struct svm_node, elements);			//data of all individuals
	
	int j = 0;

	for (int i = 0; i < psize; i++)
	{
		//copy the code of individual to temporary vector
		std::vector<double> code = pop.Indiv_Show(i).Code_Show();		//temporary vector for the individual's code storage

		//assign pointer to the start of the memory with data of the individual
		prob2.x[i] = &x_space2[j];

		//assign label of the individual to the problem
		prob2.y[i] = train_label[i];

		//assign data of individual to the storage node
		for (int k = 0; k < max_index; k++)
		{
			x_space2[j].index = k + 1;
			x_space2[j].value = code[k];
			++j;
		}
		x_space2[j++].index = -1;
	}

	//calculate new gamma parameter - only if defined
	if (param.gamma == 0 && max_index > 0)
		param.gamma = 1.0 / max_index;

}

/*
*Generation of the train_train label*
@param pop population for which the label will be generated
*/
std::vector<int> SVM_LABEL(const population & pop)
{

	std::vector<double> train_avg;				//vector for storaging averages of variables for each individual
	
	//get population size
	int psize = pop.Size_Show();				//population size

	//generate the empty label vector
	std::vector<int> label (psize,0);			//vector for storaging labels - output

	//calculate average of variables for each individual
	for (int i = 0; i < psize; i++)
	{
		//copy the code of the individual
		std::vector<double> v = pop.Indiv_Show(i).Code_Show();

		//calculate the sum of all varaibles
		double sum = std::accumulate(v.begin(), v.end(), 0.0) * 100;

		//divide the sum by pop size
		double m = sum / psize ;

		//push label to the average vector
		train_avg.push_back(m);
	}

	//calculate the average for all individuals
	double sum = std::accumulate(std::begin(train_avg), std::end(train_avg), 0.0);
	double avg = sum / train_avg.size();
	double accum = 0.0;

	//calculate the standard deviaton
	std::for_each(std::begin(train_avg), std::end(train_avg), [&](const double d) 
	{
		accum += (d - avg) * (d - avg);
	});
	double sd = sqrt(accum / (train_avg.size() - 1));		//standard deviaton

	//assign labels to each individual besed on its closeness to the average (depending on the collective size there are different functions
	//4 collectives
	if (n_col == 4)				
	{
		for (int i = 0; i < psize; i++)
		{
			if (train_avg[i] < (avg - 0.6*sd))
				label[i] = 1;
			else if (train_avg[i] < avg)
				label[i] = 2;
			else if (train_avg[i] < (avg + 0.6*sd))
				label[i] = 3;
			else
				label[i] = 4;
		}
	}
	//6 collectives
	else if (n_col == 6)				
	{
		for (int i = 0; i < psize; i++)
		{
			if (train_avg[i] < (avg - sd))
				label[i] = 1;
			else if (train_avg[i] < (avg - sd / 2))
				label[i] = 2;
			else if (train_avg[i] < avg)
				label[i] = 3;
			else if (train_avg[i] < (avg + sd / 2))
				label[i] = 4;
			else if (train_avg[i] < (avg + sd))
				label[i] = 5;
			else
				label[i] = 6;
		}
	}
	//8 collectives
	else if (n_col == 8)
	{
		for (int i = 0; i < psize; i++)
		{
			if (train_avg[i] < (avg - 1.1*sd))
				label[i] = 1;
			else if (train_avg[i] < (avg - 0.6*sd))
				label[i] = 2;
			else if (train_avg[i] < (avg - 0.2*sd))
				label[i] = 3;
			else if (train_avg[i] < (avg))
				label[i] = 4;
			else if (train_avg[i] < (avg + 0.2*sd))
				label[i] = 5;
			else if (train_avg[i] < (avg + 0.6*sd))
				label[i] = 6;
			else if (train_avg[i] < (avg + 1.1*sd))
				label[i] = 7;
			else
				label[i] = 8;
		}
	}
	else
	{
		std::cout << "ERROR#14: SVM - LABEL";
		system("pause");
		abort();
	}

	return label;
}
