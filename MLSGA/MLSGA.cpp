#define _CRT_SECURE_NO_WARNINGS
/***MSLGA***/
/***Original idea by Adam Sobey and Przemyslaw Grudniewski***/
/***Written by Przemyslaw Grudniewski***/
/***Updates by Przemyslaw Grudniewski***/
#include "MLSGA.h"



time_t mut_t, pop_t, SVM_t, run_t, col_t, elit_t;							//Time of Mutation, Population initialization, SVM, Run, Collectives, Elitism

extern time_t cross_t;										//Time of Crossover
extern time_t selec_t;										//Time of Selection
extern time_t PF_t;											//Time of PF creation
extern time_t save_t;										//Time of results saving

std::ofstream file;
std::ofstream graph, *graph_v;								//File for saving graph data
std::ofstream Debug;										//File for saving debugging data

short n_col;												//current number of collectives used
int dyn_tau;												//time state - for dynamic function
double dyn_t;												//time state - for dynamic function
int pop_size;												//current population size
std::string MODE;											//current hybrid mode
std::string REINIT_MODE;									//current reinitialisation mode
std::vector<float> TGM_vect;							//current maternal effect vector



float run_time, GA_time;
extern SimpleXlsx::CWorkbook book1;							//Book - for all data	
extern SimpleXlsx::CWorksheet sheet1;						//For saving data to excel

extern int nfes;											//number of function evaluations for MOEAD and MTS
extern int cons_viol_count;		//How many times constrains have been violated in the current run



int main()
{
	//Check for errors in parameters
	Error_Check();
	double IGD_val_final = 0;												//Final IGD value
	double HV_val_final = 0;										//Final HV value
	double IGD2_val_final = 0;										//Final IGD value - for dynamic
	double HV2_val_final = 0;										//Final HV value - for dynamic
	double min_fitness = 0.;											//minimum fitness for current run - only for MLS1
	double end_generation;										//generation at which solution was found

	//Create the vector of active functions
	std::vector<bool> func_active = Func_Deactivate();



	std::vector<std::ofstream> graph_v_vector;
	std::ofstream PF, PF2, IGD_stream, HV_stream;								//Output files
	//get the comment to files and folders names
	std::string comment;										//The comment to the folder and results name
	std::cout << "Comment (to the folder and results names): ";
	getline(std::cin,comment);
	
	//set the precision of results files
	graph << std::fixed;
	graph.precision(prec);
	PF << std::fixed;
	PF.precision(PF_prec);

	//initialize TGM vector
	TGM_vect = std::vector<float>( TGM_size, 0.f );
	
	//Get the starting date and time of the program
	std::string date = Date_Get();

	if ((comment[comment.size() - 1] == ' ' || comment.size() == 0))
	{
		if (AUTO_FOLDER_COMMENT == false)
		{
			std::cout << "**********************************\n"
				<< "Wrong paramters chosen. Folder name cannot end with a space\nProgram will terminate"
				<< "\n**********************************\n";
			system("pause");
			return 1;
		}
	}
	//add the comment to the date string - if not empty
	else
		date = date + comment;
	
	//Create directories
	Directory_Create(date);

	std::cout << date;

	file.open("Results/out.txt");						//File for saving times - Temporary						
	FILE *plotPipe = _popen("c:\\gnuplot\\bin\\gnuplot.exe", "w");								//Pipe for GNUplot
	Video * generation_video;										//Pointer to video
	
	//decide which performance evaluation method will be used - based on const.h
	bool IGD_on = false;		//If IGD will be calculated
	bool HV_on = false;			//If HV will be calculated

	for (short i = 0; i < perf_eval_type.size(); i++)
	{
		if (perf_eval_type[i] == "IGD")
			IGD_on = true;
		else if (perf_eval_type[i] == "HV")
			HV_on = true;
		else
			abort();		//wrong performance evaluation defined in const.h
	}

	//set up the minimum and max number of collectives
	short n_col_b;
	short n_col_e;
	short n_col_elim;
	if (MLSGA_Hybrid == true)
	{
		n_col_b = MLSGA_n_col_b;
		n_col_e = MLSGA_n_col_e;
		n_col_elim = MLSGA_n_col_elim;
	}
	else
	{
		n_col_e = n_col_b = 1;
		n_col_elim = 0;
	}


	//Calculate time and max index
	int index_r_max = Show_Iter();
	int index_r = 1;												//Index of run
		/************Initialisation******************/
	for (short iGA_mode = 0; iGA_mode < GA_mode.size(); iGA_mode++)								//loop for different hybrid modes
	{
		std::vector<std::string> MODE_vector = GA_mode[iGA_mode];
		short coevol_amt = MODE_vector.size();

		for (short MLS = n_MLS_b; MLS <= n_MLS_e; MLS++)												//loop for different MLS types
		{
				
			//Do only selected MLSt types
			//if (MLSt != 1 && MLSt != 2 && MLSt != 7)
				//continue;
			for (n_col = n_col_b; n_col <= n_col_e; n_col += 2)										//loop for different col size
			{
				for (pop_size = Pop_size_b; pop_size <= Pop_size_e; pop_size += Pop_size_step)		//loop for different pop sizes
				{
					
					//Do only selected population sizes
					//if (pop_size != 400 && pop_size != 600 && pop_size != 800 && pop_size != 1000 && pop_size != 1200 && pop_size != 1400 && pop_size != 1600 && pop_size != 1800 && pop_size != 2000)
						//continue;
					//Initialise max_gen - termination criterion
					int max_gen;
					if (T_con == "nfes")
					{
						max_gen = max_iterations / pop_size*1000;
					}
					else if (T_con == "ngen")
						max_gen = max_generations;
						
					for (short R_Mode = 0; R_Mode < Reinit_Mode.size(); R_Mode++)														//Loop for different mode types
					{
						REINIT_MODE = Reinit_Mode[R_Mode];

						for (TGM_vect[0] = TGM_m_min; TGM_vect[0] <= TGM_m_max; TGM_vect[0] += TGM_m_step)									//loop for different maternal effect
						{
							if (TGM != true)
								TGM_vect[0] = TGM_m_max;
							for (TGM_vect[1] = TGM_g_min; TGM_vect[1] <= TGM_g_max; TGM_vect[1] += TGM_g_step)							//loop for different grandmaternal effect
							{
								if (TGM != true)
									TGM_vect[1] = TGM_g_max;

								for (short ftype = n_func_b; ftype <= n_func_e; ftype++)								//loop for different functions
								{
									if (!func_active[ftype - n_func_b])
										continue;

									if ((ftype >= 7 && ftype <= 13) || (ftype >= 60 && ftype <= 61))
										continue;
									function * Func = Func_Type(ftype);													//Function initialisation

									//skip the cases for other reinitialisation modes for non dynamic functions
									if (Func[0].Time_Dep() == false && R_Mode != 0)
										continue;

									//If dynamic create the data for the first row of excel
									std::vector<SimpleXlsx::CellDataStr> excelb1;		//First row of PF storage excel file
									if (Func[0].Time_Dep())
									{
										for (int i = 0; i < Func[0].Objs(); i++)
											excelb1.push_back("f" + std::to_string(i + 1));
										for (int i = 0; i < Func[0].Vars(); i++)
											excelb1.push_back("x" + std::to_string(i + 1));
										excelb1.push_back("IGD");								
									}
									for (short ctype = 1; ctype <= n_cross; ctype++)									//loop for different crossover types
									{
										std::vector<crossover<individual> *> Crossover_vect;	//Crossover initialisation
										for (int i_coevol = 0; i_coevol < coevol_amt; i_coevol++)
											Crossover_vect.push_back(Cross_Type(ctype, MODE_vector[i_coevol]));

										for (short mtype = 1; mtype <= n_mut; mtype++)									//loop for different mutation types
										{
											//mutation<short> * Mutation = Mut_Type(mtype);								//Mutation initialisation
											std::vector<mutation<short> *> Mutation_vect;	//Crossover initialisation
											for (int i_coevol = 0; i_coevol < coevol_amt; i_coevol++)
												Mutation_vect.push_back(Mut_Type(mtype, MODE_vector[i_coevol]));

											for (short stype = 1; stype <= n_select; stype++)										//loop for different selection types
											{
												selection<individual> * Selection = Select_Type(stype);								//Selection initialisation
												//selection<individual> * Selection = Select_Type(1);
												for (float C_Prob = cross_prob_min; C_Prob < cross_prob_max + 0.001f; C_Prob += cross_prob_step)			//Loop for different crossover probabilities
												{
													for (float M_Prob = mut_prob_min; M_Prob < mut_prob_max + 0.001f; M_Prob += mut_prob_step)				//Loop for different mutation probabilities 
													//for (float M_Prob = 1.f; M_Prob <= 1.f; M_Prob *= 10.f)
													{

														float end_fitness;															//min fitness for MLSt 1
														GA_data<float> run_data(IGD_on, HV_on);													//storing output parameters
														int GA_success = 0;														//storing how many runs was successful															
														std::vector<std::vector<double>> real_PF;									//vector storing real PF data

														std::vector<std::vector<double>> IGD_storage;							//Storing the IGD values over the runs
														std::vector<std::vector<double>> HV_storage;							//Storing the HV values over the runs
														for (int iRun = 0; iRun < n_runs; iRun++)														//loop for number of runs- Last loop
														{
															Clear_Ideal_and_Nadir(Func[0].Objs());

															std::vector<double> IGD_storage_temp;								//Storing the IGD valus in the current run
															std::vector<double> HV_storage_temp;								//Storing the HV valus in the current run
															std::vector<double> IGD_val;										//Storage for different IGD values - for dynamic problems
															std::vector<double> HV_val;										//Storage for different HV values - for dynamic problems
															dyn_tau = 0;
															//open the new excel file
															/*SimpleXlsx::CWorkbook book2;												//Book - for PF data
															SimpleXlsx::CWorksheet sheet5;*/						//For saving data to excel

															if (Func[0].Name_Show() == "UDF8" || Func[0].Name_Show() == "UDF9")
																Func[0].T_Vector_Update(false);

															cons_viol_count = 0;

															//Plot real PF
															if (Func[0].Time_Dep())
															{
																dyn_t = floor(dyn_tau / (double)T_dyn) / (double)n_steps_dyn;
																//dyn_t = 0;
																real_PF = Func->Plot_PF(index_r, dyn_t);
															}
															else
															{
																dyn_t = 0;
																//real_PF = Func->function::Plot_PF(index_r);		//
																real_PF = Func->Plot_PF(index_r);
															}

															if (ONE_OBJ_OVERRIDE == true)
																//Get the min fitness for the given function
															{
																if (iRun % n_runs == 0)
																	end_fitness = Func->Min_Fitness_Get(real_PF);										//min fitness for 
																std::cout << end_fitness << "\n\n";
																min_fitness = 1e15;
															}
															//Get the beginning time of the run
															run_t = clock();

															//Initialise the time variables
															mut_t = pop_t = cross_t = SVM_t = PF_t = SVM_t = selec_t = elit_t = col_t = save_t = 0;

															//Initialise the GA_Parameters
															double M_Prob_val = M_Prob;
															//For binary the value of mutation chance is divided by number of variables
															if (ENCODING == "Binary" || ENCODING == "Gray")
																M_Prob_val /= (Func[0].Vars()*Binary_string_size);
															GA_parameters GA_standard{ C_Prob, M_Prob_val, pop_size, max_gen };

															//Open the file for saving PF results
															PF.open(F_Name_Get("Results/", date, "/") + Name_Get("PF.csv", std::string(), index_r, dyn_t));
															if (!PF.is_open())
															{
																short count = 0;
																while (!PF.is_open())
																{
																	PF.open(F_Name_Get("Results/", date, "/") + Name_Get("PF.csv", std::string(), index_r, dyn_t));
																	count++;
																	if (count > 500)
																		abort();
																}
															}
															//Open the file for saving IGD results - if needed
															if (PERF_VAL_GEN == true )
															{
																if (IGD_on)
																	IGD_stream.open(F_Name_Get("Results/", date, "/") + Name_Get("IGD.csv", std::string(), index_r));
																if (HV_on)
																	HV_stream.open(F_Name_Get("Results/", date, "/") + Name_Get("HV.csv", std::string(), index_r));
															}

															//Open the video file and set min and max values
															generation_video = new Video{ plotPipe,"Graph_v_" + std::to_string(index_r),F_Name_Get("Results/", date, "/") + Name_Get("G_Video.avi",std::string(), index_r), Func->Objs() };					// Video class for generations
															generation_video->Set_Max_Min(Func->Max_Min_Fit());


															//Open the file for saving temprary graph data
															graph.open("Temp/" + Name_Get("graph", index_r, dyn_t));

															//Create the First rows of the axcel
															Excel_F_Row(index_r, MODE_vector, Func[0].Time_Dep(), IGD_on, HV_on);

															//Create the Second rows of the axcel
															Excel_S_Row(*Func, GA_standard, Mutation_vect[0][0], *Selection, Crossover_vect[0][0], iRun, index_r, REINIT_MODE, MODE_vector, MLS);
															
															/***** GA Initialisation *****/
															nfes = 0;

															int coevol_subpop_size = GA_standard.Pop_Size() / coevol_amt;

															//initialise the Pareto front
															//Initialise the empty population
															population p0{ *Func,GA_standard, Mutation_vect[0][0],Selection[0],Crossover_vect[0][0], 0 };
															//Initialise the Pareto Front - empty
															pareto_front Pareto_F{ p0 };

															//initialise the vector of collectives
															std::vector<collective> col;

															//initialise the vector of temporary collectives - for coevolution
															std::vector <std::vector<collective>> coevol_col_vect;

															time_t col_t_temp = clock();

															for (int i_coevol = 0; i_coevol < coevol_amt; i_coevol++)
															{
																MODE = MODE_vector[i_coevol];

																int col_temp_size = coevol_subpop_size;
																if (i_coevol + 1 == coevol_amt)
																	col_temp_size = GA_standard.Pop_Size() - ((coevol_amt-1)*coevol_subpop_size);

																//Initialise the population
																population p1{ *Func,GA_standard, Mutation_vect[i_coevol][0],*Selection,Crossover_vect[i_coevol][0], col_temp_size };
																//Clustering(p1);

																if (MODE == "MOEAD" || MODE == "MOEADMSF" || MODE == "MOEADPSF")
																{
																	//initialise the MOEAD
																	MOEAD::MOEAD_Init(Func[0].Objs(), p1);
																}
																else if (MODE == "PAES")
																{
																	//initialise the PAES
																	PAES::PAES_Init(Func[0].Objs(), p1);
																}
																else if (MODE == "MTS")
																{
																	//Initialise MTS
																	MTS::MTS_Init(p1, Pareto_F, n_col);
																}
																else if (MODE == "BCE")
																{
																	//Initialise BCE
																	BCE::BCE_Init(Func[0].Objs(), n_col, p1);
																}
																else if (MODE == "HEIA")
																{
																	//Initialise HEIA
																	HEIA::HEIA_Init(n_col, p1);
																}
																else if (MODE == "MOEADM2M")
																{
																	//Initialise M2M
																	M2M::M2M_Init(Func[0].Objs(), p1);
																}
																else if (MODE == "IBEA")
																	IBEA::IBEA_Init(n_col, p1);



																//Initialise the vector for storaging real labels
																//std::vector<int> real_label(GA_standard.Pop_Size(), 0);										//Real labels for population separation
																std::vector<short> real_label;
																//Get the starting time of the SVM
																time_t SVM_t_temp = clock();

																//Do SVM
																//If only 1 collective SVM is pointless
																if (n_col > 1)
																{
																	if (GROUPING == "SVM")
																		real_label = SVM(p1);
																	else if (GROUPING == "Cluster")
																		real_label = Clustering(p1);
																	else if (GROUPING == "Random")
																		real_label = Random_Label(p1);
																	else
																		abort();
																	//real_label = SVM_LABEL(p1);
																	//real_label = Random_Label(p1);
																}
																else if (n_col < 0)
																{
																	std::cout << "**********************************\n"
																		<< "Wrong paramters chosen. Number of collective cannot be negative. Check const.h\nProgram will terminate"
																		<< "\n**********************************\n";
																	system("pause");
																	abort();
																}
																else if (n_col == 1)
																{
																	real_label = std::vector<short>(GA_standard.Pop_Size(), 1);
																}

																//Calculate the time of SVM
																SVM_t += clock() - SVM_t_temp;

																std::cout << "Run #" << index_r << " started\n";

																std::cout << "\nSVM - complete!\n";

																/***** Collevtive generation *****/

																//Initialise the vector of temporary collectives
																std::vector<collective> temp_col;																	//vector for storage of collectives
																//int sum_c = 0;
																
																std::vector<short> biggest_col{ 0,0,0,0 };														//Vector for storing which collectives are the biggest ones
																graph_v_vector.clear();
															
																for (short iCol = 0; iCol < n_col; iCol++)														//Loop for collective generation
																{
																	//Assign the individuals to the collective
																	temp_col.push_back(collective::collective(p1, real_label, iCol + 1, MODE));

																	if (VIDEO == true && FITNESS_ALL == true && Func[0].Objs() == 2)
																		//Open the data file for the first frame of the video
																		graph_v_vector.push_back(std::ofstream("Temp/" + Name_Get("graph_v_" + std::to_string(index_r) + "_c" + std::to_string(iCol + 1), 1)));

																	//sum_c += temp_col[iCol].Size_Show();

																	//Find the biggest collective
																	if (n_col != 1)
																	{
																		if (temp_col[iCol].Size_Show() > biggest_col[1])
																		{
																			biggest_col[2] = biggest_col[0];
																			biggest_col[3] = biggest_col[1];
																			biggest_col[0] = iCol;
																			biggest_col[1] = temp_col[iCol].Size_Show();
																		}
																	}
																}
																//Rearrange collectives to remove too small ones
																for (short iCol = 0; iCol < n_col; iCol++)
																{
																	//Get the size of the collective
																	int temp_i = temp_col[iCol].Size_Show();

																	if (n_col != 1)
																	{
																		//Collective cannot be empty or too small
																		//If is too small

																		if (temp_i <= (col_temp_size / min_col_size_multi) || temp_i <= 2)
																		{
																			//Check which collective is the biggest
																			if (biggest_col[1] >= biggest_col[3])
																			{
																				//Cut the individuals from the biggest collective to the new one
																				for (int i = temp_i; i < (int)(col_temp_size / new_col_size_multi); i++)
																				{
																					//Get the individual to be cut
																					individual temp_indi = temp_col[biggest_col[0]].Indiv_Show(0);		//Inidividual which will be copied

																					//Add the individual to the new collective
																					temp_col[iCol].Add(temp_indi);

																					//Remove it from the old one
																					temp_col[biggest_col[0]].Remove();
																				}
																				//Correct the size of the biggest collective
																				biggest_col[1] = temp_col[biggest_col[0]].Size_Show();
																			}
																			else
																			{

																				for (int i = temp_i; i < (int)(col_temp_size / new_col_size_multi); i++)
																				{
																					individual temp_indi = temp_col[biggest_col[2]].Indiv_Show(0);
																					temp_col[iCol].Add(temp_indi);
																					temp_col[biggest_col[2]].Remove();
																				}
																				biggest_col[3] = temp_col[biggest_col[2]].Size_Show();
																			}
																		}
																	}
																}

																//reduce the amount of collectives
																//find the desired size
																short n_col_target = n_col / coevol_amt;
																if (n_col % coevol_amt != 0)
																{
																	if (n_col % coevol_amt >= i_coevol + 1)
																		n_col_target++;
																}
																//reduce
																while (temp_col.size() > n_col_target)
																{
																	short ix_small = 0;
																	int ix_small_size = temp_col[0].Size_Show() + temp_col[1].Size_Show();
																	//find the smalles pair- they have to be near each other in order to maintain
																	for (short iCol = 1; iCol < temp_col.size() - 1; iCol++)
																	{
																		int temp_size = temp_col[iCol].Size_Show() + temp_col[iCol + 1].Size_Show();

																		if (temp_size < ix_small_size)
																			ix_small = iCol;
																	}
																	//merge the smallest pair
																	temp_col[ix_small].Add(temp_col[ix_small + 1].Indiv_Show());

																	//remove the added collective
																	temp_col.erase(temp_col.begin() + ix_small + 1);

																}

																//check the size of population
																int temp_size = 0;
																for (short iCol = 0; iCol < temp_col.size(); iCol++)
																{
																	temp_col[iCol].Index_Set(iCol + 1);
																	temp_size += temp_col[iCol].Size_Show();
																}
																
																if (temp_size != col_temp_size && DEBUG == true)
																	abort();
																															
																//push the collectives to the overall vector
																coevol_col_vect.push_back(temp_col);
															}

															//Initialise the vector of collectives
															if (coevol_amt == 1)
															{
																col = coevol_col_vect[0];
																
															}
															else
															{
																short ix_col = 0;
																for (int i_coevol = 0; i_coevol < coevol_amt; i_coevol++)
																{
																	for (int i_col = 0; i_col < coevol_col_vect[i_coevol].size(); i_col++)
																	{
																		col.push_back(coevol_col_vect[i_coevol][i_col]);
																		col[ix_col].index_vid = ix_col + 1;																		
																		ix_col++;
																	}
																	
																}
															}
															MLS_Col_Fit_Create(col, MLS);
															
															for (short iCol = 0; iCol < col.size(); iCol++)
															{

																MODE = col[iCol].Mode_Show();

																if (coevol_amt != 1)
																	col[iCol].Index_Set(col[iCol].index_vid);
																if (VIDEO == true && FITNESS_ALL == true && Func[0].Objs() == 2)
																{
																	//Set the precision of the video data files
																	graph_v_vector[iCol] << std::fixed;
																	graph_v_vector[iCol].precision(prec);
																	graph_v = &graph_v_vector[iCol];
																}

																if (MODE == "Normal" || MODE == "NSGAII" || MODE == "DMEADD")
																{
																	//Calculate the fitness of individuals in the collective
																	col[iCol].population::Fitness_Calc();
																	nfes += col[iCol].Size_Show();
																	if (MODE == "NSGAII" || MODE == "DMOEADD")
																		col[iCol].save();
																}
																else if (MODE == "MOEAD" || MODE == "MOEADMSF" || MODE == "MOEADPSF")
																	//initialise the population for MOEAD
																	MOEAD::Population_Init(Func[0].Objs(), col[iCol]);
																else if (MODE == "PAES")
																	//initialise the population for PAES
																	PAES::PAES_Population_Init(Func[0].Objs(), col[iCol]);
																else if (MODE == "MTS")
																	//initialise the collective parameters for MTS
																	MTS::MTS_Init_Col(col[iCol]);
																else if (MODE == "BCE")
																{
																	//initialise the collective parameters for BCE
																	BCE::BCE_Pop_Init(Func[0].Objs(), col[iCol]);
																}
																else if (MODE == "HEIA")
																{
																	//Initialise HEIA
																	HEIA::HEIA_Pop_Init(col[iCol]);
																}
																else if (MODE == "MOEADM2M")
																{
																	//initialise the collective parameters for M2M
																	M2M::M2M_Population_Init(Func[0].Objs(), col[iCol]);
																}
																else if (MODE == "IBEA")
																{
																	//initialise the collective parameters for IBEA
																	IBEA::IBEA_Pop_Init(Func[0].Objs(), col[iCol]);
																}

																if (MODE == "MOEADPSF" || MODE == "MOEADMSF" || MODE == "MOEADM2M" || MLSGA_norm_obj == true)
																{
																	if (iCol == 0)
																	{
																		Create_Nadir(Func[0].Objs());
																		Update_Nadirpoint(col[iCol].Indiv_Show(), 1);
																	}
																	else if (col[iCol - 1].Mode_Show() != "MOEADPSF" & col[iCol-1].Mode_Show() != "MOEADMSF" & col[iCol-1].Mode_Show() != "MOEADM2M" & MLSGA_norm_obj != true)
																	{
																		Create_Nadir(Func[0].Objs());
																		Update_Nadirpoint(col[iCol].Indiv_Show(), 1);
																	}
																	else
																		Update_Nadirpoint(col[iCol].Indiv_Show(), 2);
																}

																if (ONE_OBJ_OVERRIDE != true)
																{
																	//Find the PF
																	Pareto_F.Pareto_Search(col[iCol]);			//PF search
																}
																else
																{
																	//calculate the fitness of the collective and grab the min one
																	col[iCol].Fitness_Calc();
																	min_fitness = col[iCol].Min_Fitness_Show(min_fitness);
																}

															}



															//check if size is ok
															if (col.size() != n_col)
																abort();
															else if (DEBUG == true)
															{
																int col_size = 0;
																for (short iCol = 0; iCol < n_col; iCol++)
																{
																	col_size += col[iCol].Size_Show();
																}
																if (col_size != GA_standard.Pop_Size())
																	abort();
															}
															int pop_size_temp = 0;
															if (DEBUG == true)
															{
																for (short iCol = 0; iCol < n_col; iCol++)
																	pop_size_temp += col[iCol].Size_Show();
																if (pop_size_temp != pop_size)
																	abort();
															}

															//Initialize the population memory in case of tgmemory, or reinitialisation
															std::vector<std::vector<collective>> col_past;		//memory of the previous populations

															//Perform only for reinitialisation active
															if (REINIT_MODE != "None" && REINIT_MODE != "Random")
															{
																for (int i = 0; i < gen_memory_size; i++)
																	col_past.push_back(col);
															}


															//If all IGD values have to be calculated
															if (ONE_OBJ_OVERRIDE != true && (PERF_VAL_GEN == true))
															{
																if (IGD_on == true)
																{
																	//calculate IGD value
																	double IGD_val_temp = IGD_calc(Pareto_F, real_PF);

																	//save to file
																	IGD_stream << 1 << " " << IGD_val_temp << std::endl;

																	//Push the value tu the storage vector
																	IGD_storage_temp.push_back(IGD_val_temp);
																}

																if (HV_on == true)
																{
																	//calculate HV value
																	double HV_val_temp = HV::HV_calc(Pareto_F, real_PF);

																	//save to file
																	HV_stream << 1 << " " << HV_val_temp << std::endl;

																	//Push the value tu the storage vector
																	HV_storage_temp.push_back(HV_val_temp);
																}

															}
															col_t += clock() - col_t_temp;
															if (VIDEO == true && FITNESS_ALL == true && Func[0].Objs() == 2)
															{
																//close all of the video data files
																for (short i = 0; i < n_col; i++)
																	graph_v_vector[i].close();

																//Make the first frame of the video
																generation_video->Make_Frame(1, index_r, dyn_t);
															}

															if (ONE_OBJ_OVERRIDE == true)
																//copy min fitness to output file
																PF << 1 << " " << min_fitness << std::endl;


															std::cout << "Collective generation - complete!\n\n";


															std::vector<individual> PF_storage;			//storage of potential PF solutions - only used for MOEAD

															unsigned short skip_gen = MLSGA_col_elim_limit;

															//override for MTS
															if (MODE == "MTS")
																skip_gen = 1;

															/************Genetic Algorithm***************/
															for (int iGen = 1; iGen < GA_standard.Max_Gen(); iGen++)	//loop for number of generations
															{
																if (DEBUG == true)
																	std::cout << iGen + 1 << std::endl;
																double IGD_val_temp;
																double HV_val_temp;
																int dyn_change_param = 0;
																//If dynamic function is used
																if (Func[0].Time_Dep())
																{

																	//If function require additional "initial" generations
																	if (ftype >= 29 && ftype <= 38 && iGen < JY_const_window)
																	{
																		//set both values to 0
																		dyn_tau = 0;
																		dyn_t = 0;
																	}
																	else
																	{
																		//if (ftype >= 29 && ftype <= 38)
																			//set dyn_tau to current generation - predef number + one state
																			//dyn_tau = iGen - 100 + T_dyn;
																		if (T_con == "ngen")
																		{
																			//set dyn_tau to current generation
																			dyn_tau = iGen;

																			//calculate the change parameter
																			dyn_change_param = dyn_tau % T_dyn;
																		}
																		else if (T_con == "nfes")
																		{
																			//set dyn_tau to current iteration number
																			dyn_tau = nfes;
																			dyn_change_param = dyn_tau % (T_dyn*pop_size);

																			//std::cout << iGen << "; " << nfes << "; " << dyn_change_param << std::endl;
																			//chack according to error
																			if (dyn_change_param <= 0.1*pop_size || dyn_change_param >= (T_dyn*pop_size) - 0.1*pop_size)
																				dyn_change_param = 0;
																		}

																		//for moead the number of iterations is lower per gneration so need more time before time change
																		//if (MODE == "MOEAD" && )
																			//dyn_tau /= MOEAD_new_sol_multi;

																		//if timestep change occur - second part of if check is for the MOEA/D only
																		if (dyn_change_param == 0)
																		{

																			if (PERF_VAL_GEN == false)
																			{

																				Pareto_F.Pareto_Search(population{ PF_storage, col[0] });
																				PF_storage.clear();

																				if (IGD_on)
																				{
																					//Calculate the IGD value
																					IGD_val_temp = IGD_calc(Pareto_F, real_PF);
																					if (IGD_val_temp != IGD_val_temp)
																						system("pause");
																					//save the IGD to the file
																					IGD_stream << iGen + 1 << " " << IGD_val_temp << std::endl;

																					//Push the value tu the storage vector
																					IGD_storage_temp.push_back(IGD_val_temp);
																				}
																				if (HV_on)
																				{
																					//Calculate the HV value
																					HV_val_temp = HV::HV_calc(Pareto_F, real_PF);
																					if (HV_val_temp != HV_val_temp)
																						system("pause");
																					//save the HV to the file
																					HV_stream << iGen + 1 << " " << HV_val_temp << std::endl;

																					//Push the value tu the storage vector
																					HV_storage_temp.push_back(HV_val_temp);
																				}
																			}
																			if (IGD_on)
																				IGD_val.push_back(IGD_val_temp);

																			if (SKIP_GRAPHS != true || iRun <= 4)
																			{

																				//Save the current PF
																				if (ONE_OBJ_OVERRIDE != true)
																				{
																					time_t save_t_temp = clock();
																					//Save the PF
																					for (int i = 0; i < Pareto_F.Size_Show(); i++)
																					{
																						//Save objectives values
																						for (short j = 0; j < Func->Objs(); j++)
																							PF << Pareto_F.Indiv_Show(i).Fitness_Show(j) << " ";
																						//Save variables values
																						for (int k = 0; k < Func->Vars(); k++)
																							PF << "," << Pareto_F.Indiv_Show(i).Code_Show(k);
																						PF << std::endl;
																					}
																					//save time
																					save_t += clock() - save_t_temp;
																				}

																				//close the data files
																				graph.close();
																				PF.close();
																			}
																			//Print graphs
																			Print(plotPipe, index_r, index_r_max, date, dyn_t, IGD_on, HV_on);

																			//Update the t_values - for UDF8 and 9
																			if (Func[0].Name_Show() == "UDF8" || Func[0].Name_Show() == "UDF9" || Func[0].Name_Show() == "CDF6")
																				Func[0].T_Vector_Update();

																			if (T_con == "ngen")
																			{
																				//calcualte current time
																				dyn_t = floor(dyn_tau / (double)T_dyn) / (double)n_steps_dyn;
																			}
																			else if (T_con == "nfes")
																			{
																				dyn_t += 1. / (double)n_steps_dyn;
																			}




																			//plot new PF
																			real_PF = Func->Plot_PF(index_r, dyn_t);

																			if (SKIP_GRAPHS != true || iRun <= 4)
																			{
																				//Open the new file for saving PF results
																				PF.open(F_Name_Get("Results/", date, "/") + Name_Get("PF.csv", std::string(), index_r, dyn_t));

																				//Open the new file for saving temprary graph data
																				graph.open("Temp/" + Name_Get("graph", index_r, dyn_t));
																			}

																		}

																	}
																}

																//Increase the index of the last generation
																end_generation = iGen + 1;

																//open graph vectors
																if (VIDEO == true && FITNESS_ALL == true && Func[0].Objs() == 2)
																{
																	for (short iCol = 0; iCol < n_col; iCol++)											//Loop for collectives
																	{

																		//Get the index of the current collective
																		short col_index = col[iCol].index_vid - 1;
																		//Open the video data file
																		graph_v_vector[col_index].open("Temp/" + Name_Get("graph_v_" + std::to_string(index_r) + "_c" + std::to_string(col_index + 1), iGen + 1));

																	}
																}

																for (short iCol = 0; iCol < n_col; iCol++)											//Loop for collectives
																{
																	MODE = col[iCol].Mode_Show();


																	//Reinitialisation/fitness calculation for dynamic change - only for dynamic problems
																	if (dyn_change_param == 0 && Func[0].Time_Dep() && !((ftype >= 29 && ftype <= 38) && iGen < JY_const_window))
																	{

																		//Reinitialisation routine
																		Reinit(col, Pareto_F, REINIT_MODE, MODE, graph_v_vector, col_past);
																		iCol = n_col;
																		continue;

																	}

																	//Get the index of the current collective
																	short col_index = col[iCol].index_vid - 1;

																	if (VIDEO == true && FITNESS_ALL == true && Func[0].Objs() == 2)
																	{
																		//Set the pointer to the current video data file
																		graph_v = &graph_v_vector[col_index];
																	}

																	if (MODE == "Normal")
																	{

																		time_t elit_t_temp = clock();
																		//Create the elite population
																		col[iCol].Elite_Create();
																		elit_t += clock() - elit_t_temp;
																		time_t pop_t_temp = clock();

																		//Create the offspring collective via selection and crossover
																		collective col_temp{ col[iCol].Crossover(),col[iCol] };
																		pop_t += clock() - pop_t_temp;
																		time_t mut_t_temp = clock();

																		//Mutate the offspring collective
																		col_temp.Mutation();
																		mut_t += clock() - mut_t_temp;
																		time_t col_t_temp = clock();

																		//Calculate the fitness of the offspring population
																		col_temp.population::Fitness_Calc();

																		nfes += col[iCol].Indiv_Show().size();

																		//Replace the old collective
																		col[iCol] = col_temp;									//Replacement
																		col_t += clock() - col_t_temp;



																		//Check if multi objective optimisation
																		if (ONE_OBJ_OVERRIDE != true)
																		{
																			//Search for pareto front
																			Pareto_F.Pareto_Search(col[iCol]);		//PF search
																		}
																		elit_t_temp = clock();
																		//Insert the elite population
																		col[iCol].Elite_Replace();

																		elit_t += clock() - elit_t_temp;

																		if (ONE_OBJ_OVERRIDE == true)
																			//Find min fitness for the collective
																			min_fitness = col[iCol].Min_Fitness_Show(min_fitness);
																	}
																	else if (MODE == "NSGAII")
																	{
																		if (ONE_OBJ_OVERRIDE == true)
																		{
																			std::cout << "**********************************\n"
																				<< "Wrong MODE chosen. Cannot use one objective to NSGAII. Check Define.h\nProgram will terminate"
																				<< "\n**********************************\n";
																			system("pause");
																			abort();
																		}

																		//Calculate the offspring collective
																		NSGAII::NSGAII_Calc(col[iCol], iGen);

																		//Get PF
																		//Check if multi objective optimisation
																		if (ONE_OBJ_OVERRIDE != true && (Func[0].Time_Dep() || PERF_VAL_GEN == true || MLS_OVERRIDE == false || (MLS_OVERRIDE == true && n_col != 1)))
																		{
																			//Search for pareto front
																			Pareto_F.Pareto_Search(col[iCol]);		//PF search
																		}
																	}
																	else if (MODE == "MOEAD" || MODE == "MOEADMSF" || MODE == "MOEADPSF" || MODE == "MOEADM2M")
																	{

																		std::vector <individual> temp_PF_storage;
																		//calculate new population and save potential PF to temp storage vector
																		if (MODE != "MOEADM2M")
																			temp_PF_storage = MOEAD::MOEAD_Calc(col[iCol], iGen);
																		else
																			temp_PF_storage = M2M::M2M_Calc(col[iCol], iGen);
																		//insert temp storage vector to strage vector
																		PF_storage.insert(PF_storage.end(), temp_PF_storage.begin(), temp_PF_storage.end());


																		if (iGen % 50 == 0)
																		{
																			MOEAD::Utility_Comp(col[iCol], iGen);
																			//iGen = 0;
																		}

																	}
																	else if (MODE == "MTS")
																	{
																		if (ONE_OBJ_OVERRIDE == true)
																		{
																			std::cout << "**********************************\n"
																				<< "Wrong MODE chosen. Cannot use one objective to NSGAII. Check Define.h\nProgram will terminate"
																				<< "\n**********************************\n";
																			system("pause");
																			abort();
																		}
																		MTS::MTS_Calc(col[iCol]);

																		//Check if max number of fit evaluations reached
																		if ((nfes >= max_iterations || (nfes % (pop_size *T_dyn) == 0 && Func[0].Time_Dep())) && T_con == "nfes")
																		{

																			//stop the process
																			iCol = n_col;
																			break;
																		}
																	}
																	else if (MODE == "DMOEADD")
																	{
																		//calculate new population and save potential PF to temp storage vector
																		std::vector <individual> temp_PF_storage = DMOEADD::DMOEADD_Calc(col[iCol]);
																		//insert temp storage vector to strage vector
																		PF_storage.insert(PF_storage.end(), temp_PF_storage.begin(), temp_PF_storage.end());

																	}
																	else if (MODE == "PAES")
																	{
																		//Get the index of the current collective
																		short col_index = col[iCol].index_vid - 1;


																		if (VIDEO == true && FITNESS_ALL == true && Func[0].Objs() == 2)
																		{
																			//Open the video data file
																			graph_v_vector[col_index].open("Temp/" + Name_Get("graph_v_" + std::to_string(index_r) + "_c" + std::to_string(col_index + 1), iGen + 1));


																			//Set the pointer to the current video data file
																			graph_v = &graph_v_vector[col_index];
																		}

																		//calculate new population and save potential PF to temp storage vector
																		std::vector <individual> temp_PF_storage = PAES_Calc(col[iCol], 0, Selection[0]);
																		//insert temp storage vector to strage vector
																		PF_storage.insert(PF_storage.end(), temp_PF_storage.begin(), temp_PF_storage.end());
																	}
																	else if (MODE == "BCE")
																	{
																		if (n_col == 1)
																			PF_storage.clear();

																		//calculate new population and save potential PF to temp storage vector
																		std::vector <individual> temp_PF_storage = BCE::BCE_Calc(col[iCol], iGen);
																		//insert temp storage vector to strage vector
																		PF_storage.insert(PF_storage.end(), temp_PF_storage.begin(), temp_PF_storage.end());

																		if (iGen % 50 == 0)
																		{
																			MOEAD::Utility_Comp(col[iCol], iGen);
																			//iGen = 0;
																		}

																	}
																	else if (MODE == "HEIA")
																	{
																		if (n_col == 1)
																			PF_storage.clear();

																		//calculate new population and save potential PF to temp storage vector
																		std::vector <individual> temp_PF_storage = HEIA::HEIA_Calc(col[iCol], iGen);
																		//insert temp storage vector to strage vector
																		PF_storage.insert(PF_storage.end(), temp_PF_storage.begin(), temp_PF_storage.end());
																	}
																	else if (MODE == "IBEA")
																	{
																		if (n_col == 1)
																			PF_storage.clear();

																		//calculate new population and save potential PF to temp storage vector
																		std::vector <individual> temp_PF_storage = IBEA::IBEA_Calc(col[iCol], iGen);
																		//insert temp storage vector to strage vector
																		PF_storage.insert(PF_storage.end(), temp_PF_storage.begin(), temp_PF_storage.end());
																	}
																	else
																	{
																		std::cout << "**********************************\n"
																			<< "Wrong MODE chosen. Check Define.h\nProgram will terminate"
																			<< "\n**********************************\n";
																		system("pause");
																		if (DEBUG == true)
																			abort();
																		else
																			return -1;
																	}
																	//Check the stopping criterion and size of PF storage for some versions
																	if (MODE == "MOEAD" || MODE == "MOEADMSF" || MODE == "MOEADPSF" || MODE == "MOEADM2M" || MODE == "DMOEADD" || MODE == "PAES" || MODE == "BCE" || MODE == "HEIA" || MODE == "IBEA")
																	{

																		//Check if multi objective optimisation
																		if (ONE_OBJ_OVERRIDE != true && (n_col != 1))   //(n_col != 1 || (n_col == 1 && Func[0].Cons() >0)))
																		{	//Search for pareto front
																			//Check potential PF only when reach certain size
																			if (PF_storage.size() >= 200)
																			{
																				Pareto_F.Pareto_Search(population{ PF_storage, col[iCol] });
																				PF_storage.clear();
																			}

																		}
																		//if original algorithm clear the PF_storage as only the final solution matters
																		else if (ONE_OBJ_OVERRIDE != true && n_col == 1 && !Func[0].Time_Dep() && !PERF_VAL_GEN && MODE != "BCE" && MODE != "HEIA" && MODE != "IBEA") // && Func[0].Cons() >0)
																		{
																			PF_storage.clear();
																		}

																		//Check if max number of fit evaluations reached
																		if ((nfes >= max_iterations || (nfes % (pop_size *T_dyn) == 0 && Func[0].Time_Dep())) && T_con == "nfes")
																		{
																			//Check potential PF
																			if (ONE_OBJ_OVERRIDE != true && (n_col != 1)) // || (n_col == 1 && Func[0].Cons() >0)))
																			{
																				Pareto_F.Pareto_Search(population{ PF_storage, col[iCol] });
																				PF_storage.clear();
																			}
																			//stop the process
																			iCol = n_col;
																			break;
																		}
																	}

																	if (MODE == "MOEADPSF" || MODE == "MOEADMSF" || MODE == "MOEADM2M" || MLSGA_norm_obj == true)
																	{
																		if (iCol == 0)
																		{
																			Create_Nadir(Func[0].Objs());
																			Update_Nadirpoint(col[iCol].Indiv_Show(), 1);
																		}
																		else if (col[iCol - 1].Mode_Show() != "MOEADPSF" & col[iCol - 1].Mode_Show() != "MOEADMSF" & col[iCol - 1].Mode_Show() != "MOEADM2M" & MLSGA_norm_obj != true)
																		{
																			Create_Nadir(Func[0].Objs());
																			Update_Nadirpoint(col[iCol].Indiv_Show(), 1);
																		}
																		else
																			Update_Nadirpoint(col[iCol].Indiv_Show(), 2);
																	}

																}	//End of loop for collectives


																time_t col_t_temp = clock();

																if (SKIP_LAST_GEN)
																{
																	if (nfes >= 200000)
																	{
																		skip_gen = max_iterations * 10;
																	}
																}



																/*****Collectives elimination and selection*****/
																if (iGen % skip_gen == 0 && n_col_elim != 0 && !(dyn_change_param == 0 && Func[0].Time_Dep() && !((ftype >= 29 && ftype <= 38) && iGen < JY_const_window)))
																{
																	//Sort the collectives according to fitness
																	std::sort(col.begin(), col.end(), collective::Sort);

																	/*Individuals should be sorted according to index of elimianted colelctive in order to copy the best individuals according to elimianted collective*/
																	for (short iCol = 0; iCol < n_col; iCol++)
																	{
																		/*//Sort the individuals according to fitness
																		if (MLSt == 3 || ((MLSt == 7 || MLSt == 9) && col[iCol].Index_Show() % 3 == 1))
																			//for MLS3 use MLSt sorting
																			col[iCol].population::Sort_Individuals3();
																		else if ((MLSt == 6 && col[iCol].Index_Show() % 2 == 0) || ((MLSt == 7 || MLSt == 9) && col[iCol].Index_Show() % 3 == 0))
																			//if MLS6 and index of col is even change the sorting type
																			col[iCol].population::Sort_Individuals();
																		else*/
																		col[iCol].population::Sort_Individuals2();
																	}
																	int iIndivSource = -1;	//index of individual which will be taken to new collective
																	int iIndivSourceMax = (pop_size / (n_col * 3)* Func[0].Cons());
																	for (short iElim = 0; iElim < n_col_elim; iElim++)
																	{

																		//Check if collectives are not the same
																		if (abs(col[n_col - 1].Fitness_Show() - col[n_col - n_col_elim - 1].Fitness_Show()) < pow(10, -14))
																			break;
																		short col_elim_index = n_col - 1 - iElim;
																		int col_size_temp = col[col_elim_index].Size_Show();			//size of the eliminated collective

																		MODE = col[col_elim_index].Mode_Show();

																		std::vector<std::vector<double>> temp_namda;		//vector for storing old namda values - for MOEAD
																		if (MODE == "MOEAD" || MODE == "MOEADMSF" || MODE == "MOEADPSF" || MODE == "MOEADM2M" || MODE == "PAES" || MODE == "BCE")
																		{
																			//copy namda
																			for (int i = 0; i < col_size_temp; i++)
																				temp_namda.push_back(col[col_elim_index].Indiv_Show(i).namda);
																		}
																		//Erase the worst collective
																		col[col_elim_index].Erase();
																		int iIndivNew = 0;	//amount of individuals in new collective
																		/*while (true)
																		{
																			for (short j = 0; j < (n_col - n_col_elim); j++)
																			{
																				double temp_rand = Random_I(0,col[j].Size_Show()-1);
																				col[col_elim_index].Add(col[j].Indiv_Show(temp_rand));
																				iIndivNew ++;
																				if (iIndivNew == col_size_temp)
																					break;
																			}
																			if (iIndivNew == col_size_temp)
																				break;
																		}*/
																		int k = 0;
																		//Create the new collective
																		int temp_const_limit = 0;		//limit how many individuals can be skip in collective generation - with constrains
																		while (iIndivNew < col_size_temp)
																		{

																			if (k == 0)
																				iIndivSource++;
																			if (iIndivSource >= iIndivSourceMax)
																				iIndivSource = 0;
																			if (iIndivSource < col[k].Size_Show())
																			{
																				if (!col[k].Indiv_Show(iIndivSource).Cons_Viol_Show() || temp_const_limit > pop_size / 2.)
																				{
																					col[col_elim_index].Add(col[k].Indiv_Show(iIndivSource));
																					if (MODE == "MOEAD" || MODE == "MOEADMSF" || MODE == "MOEADPSF" || MODE == "MOEADM2M" || MODE == "PAES" || MODE == "BCE")
																					{
																						col[col_elim_index].Indiv_Set(iIndivNew).namda = temp_namda[iIndivNew];
																					}
																					iIndivNew++;
																				}
																				else
																					temp_const_limit++;
																			}
																			else if (k != (n_col - n_col_elim - 1))
																			{
																				k++;
																				continue;
																			}
																			else
																			{
																				k = 0;
																				continue;
																			}
																			k = iIndivNew % (n_col - n_col_elim);

																		}
																		if (col[col_elim_index].Size_Show() <= 0 || col[col_elim_index].Size_Show() != col_size_temp)
																			abort();
																		//Calculate the fitness of the new collective
																		//col[col_elim_index].Fitness_Calc();
																	}

																}
																col_t += clock() - col_t_temp;
																if (VIDEO == true && FITNESS_ALL == true && Func[0].Objs() == 2)
																{
																	for (short i = 0; i < n_col; i++)
																		graph_v_vector[i].close();
																}
																if (Func[0].Objs() == 2)
																	generation_video->Make_Frame(iGen + 1, index_r, dyn_t);

																if (ONE_OBJ_OVERRIDE == true)
																{
																	//copy min fitness to the output file
																	PF << (iGen + 1) << " " << min_fitness << std::endl;

																	//check if stopping criterion met
																	if ((long)((min_fitness)*3000000.0) <= (long)(end_fitness*3000000.0))
																	{
																		//fill the min fitness for the output file
																		for (int i = iGen; i < GA_standard.Max_Gen(); i++)
																			PF << (i + 1) << " " << min_fitness << std::endl;

																		//set current generation to max generation - stop GA
																		iGen = GA_standard.Max_Gen();

																		//Solution found - Increase success number
																		GA_success++;

																	}

																}
																if (GENERATION_GRAPHS_DATA)
																{
																	//if (nfes == (max_iterations / 20) || nfes == (max_iterations / 10) || nfes == (max_iterations / 5) || nfes == (max_iterations / 2))
																	if (iGen == 49 || iGen == 99 || iGen == 249)
																	{
																		if (MLS_OVERRIDE && MODE_vector[0] != "Normal" && n_col == 1)
																			Pareto_F.Pareto_Search(col[0]);
																		std::string name;

																		name = "PF_" + std::to_string(iGen) + ".csv";
																		char * cstr = new char[name.length() + 1];
																		strcpy(cstr, name.c_str());
																		PF2.open(F_Name_Get("Results/", date, "/") + Name_Get(cstr, std::string(), index_r));
																		for (int g1 = 0; g1 < Pareto_F.Size_Show(); g1++)
																		{
																			for (short gg = 0; gg < Func->Objs(); gg++)
																			{
																				PF2 << Pareto_F.Indiv_Show(g1).Fitness_Show(gg) << ",";

																			}
																			PF2 << std::endl;
																		}
																		PF2.close();
																		delete[] cstr;
																	}
																}
																//Calculate the IGD every generation and save PF
																//If required
																if (PERF_VAL_GEN == true)
																{

																	Pareto_F.Pareto_Search(population{ PF_storage, col[0] });
																	PF_storage.clear();


																	if (IGD_on)
																	{
																		//Calculate the IGD value
																		IGD_val_temp = IGD_calc(Pareto_F, real_PF);
																		if (IGD_val_temp != IGD_val_temp)
																			system("pause");
																		//save the IGD to the file
																		IGD_stream << iGen + 1 << " " << IGD_val_temp << std::endl;

																		//Push the value tu the storage vector
																		IGD_storage_temp.push_back(IGD_val_temp);
																	}

																	if (HV_on)
																	{
																		//Calculate the HV value
																		HV_val_temp = HV::HV_calc(Pareto_F, real_PF);
																		if (HV_val_temp != HV_val_temp)
																			system("pause");
																		//save the HV to the file
																		HV_stream << iGen + 1 << " " << HV_val_temp << std::endl;

																		//Push the value tu the storage vector
																		HV_storage_temp.push_back(HV_val_temp);
																	}
																	//If Function is dynamic save the IGD values
																	//if (Func[0].Time_Dep())
																	//{


																		//if ((dyn_tau + 1) % T_dyn == 0)
																		//{
																			//Push the value to the storage vector
																			//IGD_val.push_back(IGD_val_temp);

																			/*//Add the new sheet to excel file
																			sheet5 = book2.AddSheet(std::to_string(dyn_t));
																			//Add the first row to excel sheet
																			sheet5.AddRow(excelb1);

																			//Fill the remaining cells
																			for (int i = 0; i < Pareto_F.Size_Show(); i++)
																			{
																				sheet5.BeginRow();
																				for (short j = 0; j < Func->Objs(); j++)
																					sheet5.AddCell(std::to_string(Pareto_F.Indiv_Show(i).Fitness_Show(j)));
																				for (int k = 0; k < Func->Vars(); k++)
																					sheet5.AddCell(std::to_string(Pareto_F.Indiv_Show(i).Code_Show(k)));

																				if (i == 0)
																					sheet5.AddCell(std::to_string(IGD_val_temp));
																				sheet5.EndRow();
																			}*/

																			//}
																		//}
																}
																//Additional termination criterion check - for overall loop of number of generations

																if (nfes >= max_iterations && T_con == "nfes")
																{
																	if (n_col != 1)
																	{
																		if (!(PERF_VAL_GEN == true || Func[0].Time_Dep()))
																		{
																			Pareto_F.Pareto_Search(population{ PF_storage, col[0] });
																			PF_storage.clear();
																		}
																	}
																	break;
																}

															} //end of loop for number of generations	

															//Find PF for the original algorithms
															if (n_col == 1 && MODE_vector[0] != "Normal")
																Pareto_F.Pareto_Search(col[0]);
															if (ONE_OBJ_OVERRIDE != true)
															{
																Pareto_F.Pareto_Search(population{ PF_storage, col[0] });
																PF_storage.clear();
																
															}

															if (PF_REFINE && !ONE_OBJ_OVERRIDE)
																Pareto_F.Pareto_Refine();
															time_t save_t_temp = clock();

															if (!FITNESS_ALL && !ONE_OBJ_OVERRIDE)
															{
																//Save the data
																for (int g = 0; g < Pareto_F.Size_Show(); g++)
																{
																	std::vector<SimpleXlsx::CellDataStr> data;
																	SimpleXlsx::CellDataStr cellStr;
																	for (short h = 1; h <= Func->Objs(); h++)
																	{
																		std::ostringstream out;
																		out << std::fixed;
																		out.precision(prec);
																		double fitness_temp = Pareto_F.Indiv_Show(g).Fitness_Show(h - 1);

																		out << fitness_temp;

																		cellStr.value = out.str();
																		if (VIDEO == true && Func[0].Objs() == 2)
																		{
																			graph << fitness_temp << " ";
																			graph_v[0] << fitness_temp << " ";
																		}
																		data.push_back(cellStr);
																		//pRange->Item[excel][h] = p1.Indiv_Show(g).Fitness_Show(h-1);
																	}
																	//excel++;
																	sheet1.AddRow(data);
																}
															}
															if (ONE_OBJ_OVERRIDE != true)
															{
																if (SKIP_GRAPHS != true || iRun <= 4)
																{
																	//Save the PF
																	for (int i = 0; i < Pareto_F.Size_Show(); i++)
																	{
																		//Save objectives values
																		for (short j = 0; j < Func->Objs(); j++)
																			PF << Pareto_F.Indiv_Show(i).Fitness_Show(j) << " ";
																		//Save variables values
																		for (int k = 0; k < Func->Vars(); k++)
																			PF << "," << Pareto_F.Indiv_Show(i).Code_Show(k);
																		PF << std::endl;
																	}
																}
																if (IGD_on)
																{
																	//Calculate the IGD value
																	double IGD_val_temp = IGD_calc(Pareto_F, real_PF);

																	//If dynamic calculate the average IGD during the run
																	if (Func[0].Time_Dep())
																	{
																		//Calculate the value for every function change
																		for (int i = 0; i < IGD_val.size(); i++)
																			IGD_val_final += IGD_val[i];
																		IGD_val_final /= IGD_val.size();

																		//Calculate the value for each iteration
																		for (int i = 0; i < IGD_storage_temp.size(); i++)
																			IGD2_val_final += IGD_storage_temp[i];
																		IGD2_val_final /= IGD_storage_temp.size();

																	}
																	//IF not dynamic calculate only the final IGD value
																	else
																		IGD_val_final = IGD_val_temp;

																	//clear the IGDs vectors
																	IGD_val.clear();
																	IGD_storage.push_back(IGD_storage_temp);
																	IGD_storage_temp.clear();
																}
																if (HV_on)
																{
																	//Calculate the HV value
																	double HV_val_temp = HV::HV_calc(Pareto_F, real_PF);

																	//If dynamic calculate the average HV during the run
																	if (Func[0].Time_Dep())
																	{
																		//Calculate the value for every function change
																		for (int i = 0; i < HV_val.size(); i++)
																			HV_val_final += HV_val[i];
																		HV_val_final /= HV_val.size();

																		//Calculate the value for each iteration
																		for (int i = 0; i < HV_storage_temp.size(); i++)
																			HV2_val_final += HV_storage_temp[i];
																		HV2_val_final /= HV_storage_temp.size();

																	}
																	//IF not dynamic calculate only the final HV value
																	else
																		HV_val_final = HV_val_temp;

																	//clear the HVs vectors
																	HV_val.clear();
																	HV_storage.push_back(HV_storage_temp);
																	HV_storage_temp.clear();
																}
															}
															save_t += clock() - save_t_temp;

															//close the data files
															graph.close();
															PF.close();
															IGD_stream.close();
															HV_stream.close();

															//Save the time values
															Time_Save(index_r, IGD_val_final, HV_val_final, min_fitness, end_generation, MODE_vector[0]);

															//Make video
															if (Func[0].Objs() == 2)
															{
																generation_video->Make_Video();
																// Video Saving
																//Print graphs
																Print(plotPipe, index_r, index_r_max, date, dyn_t, IGD_on, HV_on);
															}


															//Save the run data
															if (ONE_OBJ_OVERRIDE == true)
																run_data.Add(run_time, GA_time, IGD_val_final, HV_val_final, end_generation, IGD2_val_final, HV2_val_final, min_fitness);
															else
																run_data.Add(run_time, GA_time, IGD_val_final, HV_val_final, end_generation, IGD2_val_final, HV2_val_final);

															//Save the obtained PFs during run
															/*if (Func[0].Time_Dep())
																book2.Save(F_Name_Get("Results/", date, "/") + Name_Get("PF_dyn.xlsx", date, index_r));*/

															index_r++;
														} //end of loop for number of runs - Last loop

														//Calculate the standard deviation and average values for the current run
														run_data.Std_Dev_Calculation();

														//Calculate the average IGD over the runs for each generation number and save it to file and make graph
														if (PERF_VAL_GEN && n_runs > 1)
														{
															if (IGD_on == true)
																IGD_Gen_Calc(IGD_storage, index_r - 1, date, plotPipe);
															if (HV_on == true)
																HV_Gen_Calc(HV_storage, index_r - 1, date, plotPipe);
														}
														

														//Save data to the excel sheet
														Excel_GA_data(GA_success, run_data, Func[0].Time_Dep(), IGD_on, HV_on);

														static int temp_remove_i = 1;
														if (index_r / temp_data_size >= temp_remove_i)
														{
															system("rmdir /q /s Temp");
															while (!CreateDirectoryA("Temp/", NULL))
																system("rmdir /q /s Temp");
															std::cout << "remove success\n";
															temp_remove_i++;
														}

													} //end of loop for different mutation probabilities
												} //end of loop for different crossover probabilities
											} //end of loop for different selection types
										}//end of loop for different mutation types
									} //end of loop for different crossover types
								} //end of loop for different functions
							}//end loop for different grandmaternal effect
						}//end loop for different maternal effect
					} //end of loop for different mode types
				} //end of loop for different pop sizes
			} //end of loop for different col sizes
		}//end of loop for different MLSt types
		
	} //end of loop for different hybrid mode types
	//Save the excel sheet
	book1.Save(F_Name_Get("Results/", date, "/") + Name_Get("Results.xlsx", date, -1));
	file.close();

	fprintf(plotPipe, "exit");
	fflush(plotPipe);
	_pclose(plotPipe);
	system("rmdir /q /s Temp");
	std::cout << "Successful" << std::endl;
	system("pause");
	return 0;
}