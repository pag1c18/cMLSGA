#include "Contour_Plot.h"
#include <iostream>
#include "imported.h"
#include "Const.h"
#include "Define.h"
#include <string>
static double max_cp;	//max z value for the contour plot

/*
*Contour plot generation*
@param indexr index of the run for the current contour plot
*/
void Contour_Plot(int indexr)
{
	//Check if we have data to calculate contour plot
	if (FITNESS_ALL == false)
	{
		std::cout << "No data points for countour plot!\n"
			<< "Define FITNESS_ALL !\n";
		system("pause");
		return;
	}// !FITNESS_ALL

	//print a message
	std::cout << "Contour plot generation - may take a while!\n";


	//Initialise the kernel
	KDE* kde = new KDE();				//pointer to kernel
	std::string output, line;				//temp strings for data storage 
	kde->set_kernel_type(3);			//Kernel type (3 - is the best)
	short pdf = 1;						//parameter
	
	//open source file
	std::ifstream file("Temp/Graph_" + std::to_string(indexr) + ".x1");	//source file

	//read the input file
	while (std::getline(file, line))			
	{
		std::vector<double> data2;		//vector for data storage
		std::stringstream ss(line);		//whole line

		//read the values
		while (std::getline(ss, output, ' '))
			data2.push_back(atof(output.c_str()));
		kde->add_data(data2);
	}

	//close the source file
	file.close();

	//open the output data file
	std::ofstream file2("Temp/Graph2_" + std::to_string(indexr) + ".x1");
	{
		//copy the resolution
		short nd = c_plot_res;						//resolution of the contour plot - define step size
		
		//copy the degree
		short ad = c_plot_deg;						//contour plot degree
		
		//get min and max value for each axis
		double min_x = kde->get_min(0);				//min x-axis value
		double max_x = kde->get_max(0);				//max x-axis value
		double min_y = kde->get_min(1);				//min y-axis value
		double max_y = kde->get_max(1);				//max y-axis value

		//calculate the step for each axis
		double x_increment = (max_x - min_x) / nd;	//x-axis step
		double y_increment = (max_y - min_y) / nd;	//y-axis step
		
		double y;	//y-axis progress
		double x = min_x;	//x-axis progress

		std::vector<std::vector<double>> results;		//vector for results saving (all)
		std::vector<double> par_results;				//vector for results saving (one step)
		
		float temp2 = 0.0001f;							//temporary value for defining max z value
		int proc = 0;									//progress value						
		size_t count = 0;								//current step
		
		//print the current progress
		std::cout << "Progress: " << proc << "%";

		//calculate the z,x,y values of contour plot
		for (int i = 0; i < nd; i++, x += x_increment)			//along x axis
		{
			y = min_y;
			//calculate the z,y values for current x
			for (int j = 0; j < nd; j++, y += y_increment)		//along y axis
			{
				double temp_z;						//z value
				//next step
				count++;

				//calculate the z value for given x,y
				if (pdf == 1)
					temp_z = kde->pdf(x, y);
				else
					temp_z = kde->cdf(x, y);
				//add z value to the sum
				temp2 += temp_z;

				//push x,y and z to the vector
				par_results.push_back(x);
				par_results.push_back(y);
				par_results.push_back(temp_z);
				
				//push the current results vector to the vector with all results
				results.push_back(par_results);
				par_results.clear();

				//print the progress bar
				if (count % (nd*nd / 100) == 0)
				{
					if (proc < 10)
						std::cout << "\b\b";
					else
						std::cout << "\b\b\b";
					proc++;
					std::cout << proc << "%";
				}
			}
		}
		std::cout << std::endl;

		//copy the size of the rsults
		size_t results_size = results.size();			//size of the results
		
		//calculate the max z value of the contour plot
		temp2 = (temp2 * ad) / results_size;
		max_cp = temp2 * 1.01;
		
		//save the results to the file
		for (int i = 0; i < results_size; i++)
		{
			file2 << results[i][0] << " "
				<< results[i][1] << " ";
			if (results[i][2] > (temp2))
				file2 <<  temp2 << " " << std::endl;
			else
				file2 << results[i][2] << " " << std::endl;
			if (i!= (results_size - 1) && results[i + 1][0] != results[i][0])
				file2 << std::endl;
		}	
	}
	//close the output data file
	file2.close();

	//delete[] kde;
}

/**returning max z value for the contour plot**/
double Get_Max_CP()
{
	return max_cp;
}
/*
{
	file2 << x << " ";
	file2 << y << " ";
	if (pdf == 1) {
		file2 << temp;
	}
	else {
		file2 << kde->cdf(x, y);
	}
	file2 << endl;
}
			}
			file2 << endl;*/