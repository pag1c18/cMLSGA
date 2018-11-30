#include "functions.h"
#include "Const.h"
/*UF1 fitness and boundary calculation*/
/**********************UF1**************************/

std::vector<double> UF1::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;
	for (short i = 2; i <= Vars(); i++)
	{
		double y = code[i - 1] - sin((6 * pi * code[0]) + (i * pi) / (float)Vars());
		if (i % 2 == 1)
		{
			sizeJ1++;
			SumJ1 += pow(y, 2);
		}
		else
		{
			sizeJ2++;
			SumJ2 += pow(y, 2);
		}
	}

	std::vector<double> fitness;
	fitness.push_back(code[0] + (2.0 / sizeJ1)*SumJ1);				//fitness 1
	fitness.push_back(1 - sqrt(code[0]) + (2.0 / sizeJ2)*SumJ2);		//fitness 2
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}
	return fitness;

}

void UF1::Set_Bound()
{
	std::vector<STRUCTURES::boundaries> bound, max_min;
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);
	for (short i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}
	STRUCTURES::boundaries t1{ 6,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 6,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	function::Set_Bound(bound);
	function::Set_Max_Min_Fit(max_min);
}
std::vector<std::vector<float>> UF1::Plot_PF(short indexr)
{
	std::ofstream PF_real;
	PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + ".x1");
	
	std::vector<std::vector<float>> PF_real_val;
	for (double i = 0; i <= 1; i += 1.0 / ((double)PF_real_size - 1.0))
	{
		float temp = (1 - sqrt(i));
		std::vector<float> temp_vect;
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		PF_real_val.push_back(temp_vect);

		PF_real << i << " ";
		PF_real << temp;
		PF_real << std::endl;
	}
	PF_real.close();
	return PF_real_val;
}

/*****************UF2**********************/

std::vector<double> UF2::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;
	for (short i = 2; i <= Vars(); i++)
	{
		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += pow(code[i-1] - (0.3*pow(code[0], 2)*cos(24 * pi*code[0] + 4 * i*pi / (double)Vars()) + 0.6*code[0])*cos(6 * pi*code[0] + i*pi / (double)Vars()), 2);
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += pow(code[i-1] - (0.3*pow(code[0], 2)*cos(24 * pi*code[0] + 4 * i*pi / (double)Vars()) + 0.6*code[0])*sin(6 * pi*code[0] + i*pi / (double)Vars()), 2);
		}
	}

	std::vector<double> fitness;
	fitness.push_back(code[0] + (2 / sizeJ1)*SumJ1);				//fitness 1
	fitness.push_back(1 - sqrt(code[0]) + (2 / sizeJ2)*SumJ2);		//fitness 2
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}

	return fitness;

}

void UF2::Set_Bound()
{
	std::vector<STRUCTURES::boundaries> bound, max_min;
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);
	for (short i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}
	STRUCTURES::boundaries t1{ 4,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 4,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	function::Set_Bound(bound);
	function::Set_Max_Min_Fit(max_min);
}

std::vector<std::vector<float>> UF2::Plot_PF(short indexr)
{
	std::ofstream PF_real;
	PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + ".x1");
	
	std::vector<std::vector<float>> PF_real_val;
	for (double i = 0; i <= 1; i += 1.0 / ((double)PF_real_size - 1.0))
	{
		float temp = (1 - sqrt(i));
		std::vector<float> temp_vect;
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		PF_real_val.push_back(temp_vect);

		PF_real << i << " ";
		PF_real << temp;
		PF_real << std::endl;
	}
	PF_real.close();
	return PF_real_val;
}
/*****************UF3**********************/

std::vector<double> UF3::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2, MultipJ1, MultipJ2;
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;
	MultipJ1 = MultipJ2 = 1;
	for (short i = 2; i <= Vars(); i++)
	{
		double y = code[i-1] - pow(code[0],(0.5*(1.0+3.0*(i-2)/(double)(Vars()-2))));
		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += pow(y, 2);
			MultipJ1 *= cos((20 * y*pi) / sqrt(i));
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += pow(y, 2);
			MultipJ2 *= cos((20 * y*pi) / sqrt(i));
		}
	}

	std::vector<double> fitness;
	fitness.push_back(code[0] + (2.0 / sizeJ1)*(4.0 * SumJ1 - 2.0 * MultipJ1 + 2));				//fitness 1
	fitness.push_back(1.0 - sqrt(code[0]) + (2.0 / sizeJ2)*(4.0 * SumJ2 - 2.0 * MultipJ2 + 2));		//fitness 2
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}

	return fitness;

}

void UF3::Set_Bound()
{
	std::vector<STRUCTURES::boundaries> bound, max_min;
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);
	for (short i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}
	STRUCTURES::boundaries t1{ 6,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 6,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	function::Set_Bound(bound);
	function::Set_Max_Min_Fit(max_min);
}

std::vector<std::vector<float>> UF3::Plot_PF(short indexr)
{
	std::ofstream PF_real;
	PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + ".x1");
	
	std::vector<std::vector<float>> PF_real_val;
	for (double i = 0; i <= 1; i += 1.0 / ((double)PF_real_size - 1.0))
	{
		float temp = (1.0 - sqrt(i));
		std::vector<float> temp_vect;
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		PF_real_val.push_back(temp_vect);

		PF_real << i << " ";
		PF_real << temp;
		PF_real << std::endl;
	}
	PF_real.close();
	return PF_real_val;
}

/*****************UF4**********************/

std::vector<double> UF4::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;
	for (short i = 2; i <= Vars(); i++)
	{
		double y = code[i-1] - sin(6 * pi*code[0] + (i*pi) / (double)Vars());
		double h = abs(y) / (1.0 + exp(2 * abs(y)));
		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += h;
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += h;
		}
	}

	std::vector<double> fitness;
	fitness.push_back(code[0] + (2 / sizeJ1)*SumJ1);					//fitness 1
	fitness.push_back(1 - pow(code[0], 2) + (2 / sizeJ2)*SumJ2);		//fitness 2
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}

	return fitness;

}

void UF4::Set_Bound()
{
	std::vector<STRUCTURES::boundaries> bound, max_min;
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);
	for (short i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 2;
		temp.lower = -2;
		bound.push_back(temp);
	}
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}
	STRUCTURES::boundaries t1{ 2,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 2,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	function::Set_Bound(bound);
	function::Set_Max_Min_Fit(max_min);
}

std::vector<std::vector<float>> UF4::Plot_PF(short indexr)
{
	std::ofstream PF_real;
	PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + ".x1");
	
	std::vector<std::vector<float>> PF_real_val;
	for (double i = 0; i <= 1; i += 1.0 / ((double)PF_real_size - 1.0))
	{
		float temp = (1 - pow(i, 2));
		std::vector<float> temp_vect;
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		PF_real_val.push_back(temp_vect);

		PF_real << i << " ";
		PF_real << temp;
		PF_real << std::endl;
	}
	PF_real.close();
	return PF_real_val;
}

/*****************UF5**********************/

std::vector<double> UF5::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;
	float eps = 0.1;	// CEC '09
	short N = 10;		// CEC '09
	for (short i = 2; i <= Vars(); i++)
	{
		double y = code[i-1] - sin(6 * pi*code[0] + i*pi / (double)Vars());
		double h = 2 * pow(y, 2) - cos(4 * pi*y) + 1;
		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += h;
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += h;
		}
	}

	std::vector<double> fitness;
	double temp = (1.0 / (2.0 * N) + eps)*abs(sin(2 * N*pi*code[0]));
	fitness.push_back(code[0] + temp + (2 / sizeJ1)*SumJ1);				//fitness 1
	fitness.push_back(1 - code[0] + temp + (2 / sizeJ2)*SumJ2);			//fitness 2
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}

	return fitness;

}

void UF5::Set_Bound()
{
	std::vector<STRUCTURES::boundaries> bound, max_min;
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);
	for (short i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}
	STRUCTURES::boundaries t1{ 15,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 15,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	function::Set_Bound(bound);
	function::Set_Max_Min_Fit(max_min);
}

std::vector<std::vector<float>> UF5::Plot_PF(short indexr)
{
	short N = 10;		// CEC '09
	std::ofstream PF_real;
	PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + ".x1");
	
	std::vector<std::vector<float>> PF_real_val;
	for (double i = 0; i <= 2 * N; i++)
	{
		float temp = i / (2.0 * N);
		PF_real << temp << " ";
		PF_real << (1 - temp);
		PF_real << std::endl;

		std::vector<float> temp_vect;
		temp_vect.push_back(temp);
		temp_vect.push_back(1 - temp);
		PF_real_val.push_back(temp_vect);

	}
	PF_real.close();
	
	return PF_real_val;
}

/*****************UF7**********************/

std::vector<double> UF7::Fitness_C(const std::vector<double> & code)
{
	double sizeJ1, sizeJ2, SumJ1, SumJ2;
	sizeJ1 = sizeJ2 = SumJ1 = SumJ2 = 0;
	for (short i = 2; i <= Vars(); i++)
	{
		double y = code[i-1] - sin(6 * pi*code[0] + i*pi / (double)Vars());
		if (i % 2 == 1)
		{
			sizeJ1 += 1;
			SumJ1 += pow(y,2);
		}
		else
		{
			sizeJ2 += 1;
			SumJ2 += pow(y, 2);
		}
	}

	std::vector<double> fitness;
	fitness.push_back(pow(code[0],0.2) + (2 / sizeJ1)*SumJ1);				//fitness 1
	fitness.push_back(1 - pow(code[0],0.2) + (2 / sizeJ2)*SumJ2);		//fitness 2
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}

	return fitness;

}

void UF7::Set_Bound()
{
	std::vector<STRUCTURES::boundaries> bound, max_min;
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);
	for (short i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = -1;
		bound.push_back(temp);
	}
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}
	STRUCTURES::boundaries t1{ 6,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 6,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	function::Set_Bound(bound);
	function::Set_Max_Min_Fit(max_min);
}

std::vector<std::vector<float>> UF7::Plot_PF(short indexr)
{
	std::ofstream PF_real;
	PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + ".x1");
	
	std::vector<std::vector<float>> PF_real_val;
	for (double i = 0; i <= 1; i += 1.0 / ((double)PF_real_size - 1.0))
	{
		float temp = (1 - i);
		std::vector<float> temp_vect;
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		PF_real_val.push_back(temp_vect);

		PF_real << i << " ";
		PF_real << temp;
		PF_real << std::endl;
	}
	PF_real.close();
	return PF_real_val;
}


/**************SCH****************************/

std::vector<double> SCH::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;
	fitness.push_back(pow(code[0], 2));				//fitness 1
	fitness.push_back(pow(code[0] - 2, 2));		//fitness 2
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}

	return fitness;
}
void SCH::Set_Bound()
{
	std::vector<STRUCTURES::boundaries> bound, max_min;
	STRUCTURES::boundaries first{ 10,-10 };
	bound.push_back(first);

	STRUCTURES::boundaries t1{ pow(first.upper,2),0 };							// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ pow(first.upper + 2,2),0 };						// t2 - max and min boundary for f2
	max_min.push_back(t2);

	function::Set_Bound(bound);
	function::Set_Max_Min_Fit(max_min);
}
std::vector<std::vector<float>> SCH::Plot_PF(short indexr)
{
	std::ofstream PF_real;
	PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + ".x1");
	
	std::vector<std::vector<float>> PF_real_val;
	for (double i = 0; i <= 4; i += 4.0 / ((double)PF_real_size - 1.0))
	{
		float temp = (i - (4 * sqrt(i)) + 4);
		std::vector<float> temp_vect;
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		PF_real_val.push_back(temp_vect);

		PF_real << i << " ";
		PF_real << temp;
		PF_real << std::endl;
	}
	PF_real.close();
	return PF_real_val;
}
/**************FON****************************/
std::vector<double> FON::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;
	double sum_f1, sum_f2;
	sum_f1 = sum_f2 = 0;
	for (short i = 0; i < Vars(); i++)
	{
		sum_f1 += pow(code[i] - (1 / sqrt(Vars())),2);
		sum_f2 += pow(code[i] + (1 / sqrt(Vars())), 2);
	}
	fitness.push_back(1 - exp(-sum_f1));				//fitness 1
	fitness.push_back(1 - exp(-sum_f2));				//fitness 2
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}

	return fitness;
}
void FON::Set_Bound()
{
	std::vector<STRUCTURES::boundaries> bound, max_min;
	for (short i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 4;
		temp.lower = -4;
		bound.push_back(temp);
	}
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}

	STRUCTURES::boundaries t1{ 1,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 1,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	function::Set_Bound(bound);
	function::Set_Max_Min_Fit(max_min);
}
std::vector<std::vector<float>> FON::Plot_PF(short indexr)
{
	std::ofstream PF_real;
	PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + ".x1");
	
	std::vector<std::vector<float>> PF_real_val;
	for (double i = -1; i <= 1; i += 2.0 / ((double)PF_real_size - 1.0))
	{
		float temp2 = (1 - exp(-(pow(i - 1, 2))));
		float temp = (1 - exp(-(pow(i + 1, 2))));
		std::vector<float> temp_vect;
		temp_vect.push_back(temp2);
		temp_vect.push_back(temp);
		PF_real_val.push_back(temp_vect);

		PF_real << temp2 << " ";
		PF_real << temp;
		PF_real << std::endl;
	}
	PF_real.close();
	return PF_real_val;
}
/**************ZDT1****************************/
std::vector<double> ZDT1::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;
	double sum_f2 = 0;
	double g_x = 0;
	for (short i = 1; i < Vars(); i++)
		sum_f2 += code[i] ;
	g_x = 1 + 9 * sum_f2 / ((double)Vars() - 1);
	fitness.push_back(code[0]);								//fitness 1
	fitness.push_back(g_x*(1-sqrt((code[0]/g_x))));				//fitness 2
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}

	return fitness;
}
void ZDT1::Set_Bound()
{
	std::vector<STRUCTURES::boundaries> bound, max_min;
	for (short i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}
	STRUCTURES::boundaries t1{ 1,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	function::Set_Bound(bound);
	function::Set_Max_Min_Fit(max_min);
}

std::vector<std::vector<float>> ZDT1::Plot_PF(short indexr)
{
	std::ofstream PF_real;
	PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + ".x1");
	
	std::vector<std::vector<float>> PF_real_val;
	for (double i = 0; i <= 1; i += 1.0 / ((double)PF_real_size - 1.0))
	{
		float temp = (1 - sqrt(i));
		std::vector<float> temp_vect;
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		PF_real_val.push_back(temp_vect);

		PF_real << i << " ";
		PF_real << temp;
		PF_real << std::endl;
	}
	PF_real.close();
	return PF_real_val;
}
/**************ZDT2****************************/
std::vector<double> ZDT2::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;
	double sum_f2 = 0;
	double g_x = 0;
	for (short i = 1; i < Vars(); i++)
		sum_f2 += code[i];
	g_x = 1 + 9 * sum_f2 / ((double)Vars() - 1);
	fitness.push_back(code[0]);										//fitness 1
	fitness.push_back(g_x*(1 - pow((code[0] / g_x),2)));				//fitness 2
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}

	return fitness;
}
void ZDT2::Set_Bound()
{
	std::vector<STRUCTURES::boundaries> bound, max_min;
	for (short i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}
	STRUCTURES::boundaries t1{ 1,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	function::Set_Bound(bound);
	function::Set_Max_Min_Fit(max_min);
}

std::vector<std::vector<float>> ZDT2::Plot_PF(short indexr)
{
	std::ofstream PF_real;
	PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + ".x1");
	
	std::vector<std::vector<float>> PF_real_val;
	for (double i = 0; i <= 1; i += 1.0 / ((double)PF_real_size - 1.0))
	{
		float temp = (1 - pow(i, 2));
		std::vector<float> temp_vect;
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		PF_real_val.push_back(temp_vect);

		PF_real << i << " ";
		PF_real << temp;
		PF_real << std::endl;
	}
	PF_real.close();
	return PF_real_val;
}
/**************ZDT3****************************/
std::vector<double> ZDT3::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;
	double sum_f2 = 0;
	double g_x = 0;
	double f2,f1;
	f1 = code[0];
	for (short i = 1; i < Vars(); i++)
		sum_f2 += code[i];
	g_x = 1 + 9 * sum_f2 / ((double)Vars() - 1);
	f2 = g_x*(1 - sqrt(f1 / g_x) - f1 / g_x*sin(10 * pi*f1));
	fitness.push_back(f1);									//fitness 1
	fitness.push_back(f2);				//fitness 2
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}

	return fitness;
}
void ZDT3::Set_Bound()
{
	std::vector<STRUCTURES::boundaries> bound, max_min;
	for (short i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}
	STRUCTURES::boundaries t1{ 1,0};									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,-1 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	function::Set_Bound(bound);
	function::Set_Max_Min_Fit(max_min);
}

std::vector<std::vector<float>> ZDT3::Plot_PF(short indexr)				//non continous
{
	std::ofstream PF_real;
	PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + ".x1");

	std::vector<std::vector<float>> PF_real_val;
	double i_step = ((0.8518328654 - 0.8233317983) + (0.6525117038 - 0.6183967944) + (0.4538821041 - 0.4093136748) + (0.2577623634 - 0.1822287280) + (0.0830015349)) / ((double)PF_real_size - 1.0);
	for (double i = 0; i <= 0.0830015349; i += i_step)
	{
		float temp = (1 - sqrt(i) - i*sin(10 * pi*i));
		std::vector<float> temp_vect;
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		PF_real_val.push_back(temp_vect);

		PF_real << i << " ";
		PF_real << temp;
		PF_real << std::endl;;
	}
	for (double i = 0.1822287280; i <= 0.2577623634; i += i_step)
	{
		float temp = (1 - sqrt(i) - i*sin(10 * pi*i));
		std::vector<float> temp_vect;
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		PF_real_val.push_back(temp_vect);

		PF_real << i << " ";
		PF_real << temp;
		PF_real << std::endl;
	}
	for (double i = 0.4093136748; i <= 0.4538821041; i += i_step)
	{
		float temp = (1 - sqrt(i) - i*sin(10 * pi*i));
		std::vector<float> temp_vect;
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		PF_real_val.push_back(temp_vect);

		PF_real << i << " ";
		PF_real << temp;
		PF_real << std::endl;
	}
	for (double i = 0.6183967944; i <= 0.6525117038; i += i_step)
	{
		float temp = (1 - sqrt(i) - i*sin(10 * pi*i));
		std::vector<float> temp_vect;
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		PF_real_val.push_back(temp_vect);

		PF_real << i << " ";
		PF_real << temp;
		PF_real << std::endl;
	}
	for (double i = 0.8233317983; i <= 0.8518328654; i += i_step)
	{
		float temp = (1 - sqrt(i) - i*sin(10 * pi*i));
		std::vector<float> temp_vect;
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		PF_real_val.push_back(temp_vect);

		PF_real << i << " ";
		PF_real << temp;
		PF_real << std::endl;
	}
	
	PF_real.close();
	return PF_real_val;
}
/**************ZDT4****************************/
std::vector<double> ZDT4::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;
	double sum_f2 = 0;
	double g_x = 0;
	double f2, f1;
	f1 = code[0];
	for (short i = 1; i < Vars(); i++)
		sum_f2 += (pow(code[i], 2) - 10 * cos(4 * pi*code[i]));
	g_x = 1 + 10 * ((double)Vars() - 1) + sum_f2;
	//f2 = g_x*(1 - pow(f1 / g_x, 2));  //Website
	f2 = g_x*(1 - sqrt(f1 / g_x));  //paper
	fitness.push_back(f1);									//fitness 1
	fitness.push_back(f2);									//fitness 2
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}

	return fitness;
}
void ZDT4::Set_Bound()
{
	std::vector<STRUCTURES::boundaries> bound, max_min;
	STRUCTURES::boundaries first{ 1,0 };
	bound.push_back(first);
	for (short i = 1; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 5;
		temp.lower = -5;
		bound.push_back(temp);
	}
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}
	STRUCTURES::boundaries t1{ 1,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 150,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	function::Set_Bound(bound);
	function::Set_Max_Min_Fit(max_min);
}

std::vector<std::vector<float>> ZDT4::Plot_PF(short indexr)
{
	std::ofstream PF_real;
	PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + ".x1");

	std::vector<std::vector<float>> PF_real_val;
	for (double i = 0; i <= 1; i += 1.0 / ((double)PF_real_size - 1.0))
	{
		float temp = (1 - sqrt(i));
		std::vector<float> temp_vect;
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		PF_real_val.push_back(temp_vect);

		PF_real << i << " ";
		PF_real << temp;
		PF_real << std::endl;
	}
	PF_real.close();
	return PF_real_val;
}
/**************ZDT5****************************/
std::vector<double> ZDT5::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;
	double sum_f2 = 0;
	double g_x = 0;
	for (short i = 1; i < Vars(); i++)
		sum_f2 += code[i];
	g_x = 1 + 9 * sum_f2 / ((double)Vars() - 1);
	fitness.push_back(code[0]);								//fitness 1
	fitness.push_back(g_x*(1 - sqrt((code[0] / g_x))));				//fitness 2
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}

	return fitness;
}
void ZDT5::Set_Bound()
{
	std::vector<STRUCTURES::boundaries> bound, max_min;
	for (short i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}
	STRUCTURES::boundaries t1{ 1,0 };									// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	function::Set_Bound(bound);
	function::Set_Max_Min_Fit(max_min);
}

std::vector<std::vector<float>> ZDT5::Plot_PF(short indexr)
{
	std::ofstream PF_real;
	PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + ".x1");

	std::vector<std::vector<float>> PF_real_val;
	for (double i = 0; i <= 1; i += 1.0 / ((double)PF_real_size - 1.0))
	{
		float temp = (1 - sqrt(i));
		std::vector<float> temp_vect;
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		PF_real_val.push_back(temp_vect);

		PF_real << i << " ";
		PF_real << temp;
		PF_real << std::endl;
	}
	PF_real.close();
	return PF_real_val;
}
/**************ZDT6****************************/
std::vector<double> ZDT6::Fitness_C(const std::vector<double> & code)
{
	std::vector<double> fitness;
	double sum_f2 = 0;
	double f1, f2;
	double g_x = 0;
	for (short i = 1; i < Vars(); i++)
		sum_f2 += code[i];
	g_x = 1 + 9 * pow((sum_f2 / ((double)Vars() - 1)),0.25);
	f1 = 1 - exp(-4.0 * code[0])*pow(sin(6 * pi*code[0]), 6);
	f2 = g_x * (1 - pow((f1 / g_x), 2));			//From paper
	//f2 = 1 - pow((f1 / g_x), 2);					//From website - ETH - wrong!!

	fitness.push_back(f1);										//fitness 1
	fitness.push_back(f2);										//fitness 2
	if (fitness.size() != Objs())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE";				//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}

	return fitness;
}
void ZDT6::Set_Bound()
{
	std::vector<STRUCTURES::boundaries> bound, max_min;
	for (short i = 0; i < Vars(); i++)
	{
		STRUCTURES::boundaries temp;
		temp.upper = 1;
		temp.lower = 0;
		bound.push_back(temp);
	}
	if (bound.size() != Vars())
	{
		std::cout << "ERROR#12: FUNCTION - VECTOR SIZE; SET BOUND";		//ERROR#12: FUNCTION - VECTOR SIZE
		abort();
	}
	STRUCTURES::boundaries t1{ 1,0.28 };								// t1 - max and min boundary for f1
	max_min.push_back(t1);
	STRUCTURES::boundaries t2{ 10,0 };									// t2 - max and min boundary for f2
	max_min.push_back(t2);

	function::Set_Bound(bound);
	function::Set_Max_Min_Fit(max_min);
}
std::vector<std::vector<float>> ZDT6::Plot_PF(short indexr)
{
	std::ofstream PF_real;
	PF_real.open("Temp/" + real_PF_out + "_" + std::to_string(indexr) + ".x1");
	
	std::vector<std::vector<float>> PF_real_val;
	double i_min, i_max, i_step;
	i_min = 0.2807753191;
	i_max = 1.0;
	i_step = (i_max - i_min) / ((double)PF_real_size - 1.0);
	for (double i = i_min; i <= i_max; i += i_step)
	{
		float temp = (1 - (i*i));
		std::vector<float> temp_vect;
		temp_vect.push_back(i);
		temp_vect.push_back(temp);
		PF_real_val.push_back(temp_vect);

		PF_real << i << " ";
		PF_real << temp;
		PF_real << std::endl;
	}
	PF_real.close();
	return PF_real_val;
}
/*
void functions::collective()
{
	short pop_size = 5;
	//pop_size=size(pop_real,1);
	short col_index1 = 0;
	short col_index2 = 0;
	short col_index3 = 0;
	short col_index4 = 0;
	short col_index5 = 0;
	short col_index6 = 0;
	short col_index7 = 0;
	short col_index8 = 0;
	
	for (short ipop = 0; ipop < pop_size; i++)
	{
		if (real_label(ipop) == 1)
		{
			col_index1 = col_index1 + 1;
			collective(1).pop_real(col_index1, :) = pop_real(ipop, :);
		}
		else if (real_label(ipop) == 2)
		{
			col_index2 = col_index2 + 1;
			collective(2).pop_real(col_index2, :) = pop_real(ipop, :);
		}
		else if (real_label(ipop) == 3)
		{
			col_index3 = col_index3 + 1;
			collective(3).pop_real(col_index3, :) = pop_real(ipop, :);
		}
		else if (real_label(ipop) == 4)
		{
			col_index4 = col_index4 + 1;
			collective(4).pop_real(col_index4, :) = pop_real(ipop, :);
		}
		else if (real_label(ipop) == 5)
		{
			col_index5 = col_index5 + 1;
			collective(5).pop_real(col_index5, :) = pop_real(ipop, :);
		}
		else if (real_label(ipop) == 6)
		{
			col_index6 = col_index6 + 1;
			collective(6).pop_real(col_index6, :) = pop_real(ipop, :);
		}
		else if (real_label(ipop) == 7)
		{
			col_index7 = col_index7 + 1;
			collective(7).pop_real(col_index7, :) = pop_real(ipop, :);
		}
		else
		{
			col_index8 = col_index8 + 1;
			collective(8).pop_real(col_index8, :) = pop_real(ipop, :);
		}
	}

}
*/