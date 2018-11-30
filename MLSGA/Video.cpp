#include "Video.h"
#include "Support_Functions.h"
#include "Define.h"
#include <opencv2\opencv.hpp>
extern short n_col;
extern std::string MODE;

Video::Video() { std::cout << "ERROR - Video class cannot be empty"; system("pause"); abort(); }
/*
*Getting name for the video frame generation*
@param n generation index number
@param indexr index of the current run
*/
void Video::Name_Get(int n, int indexr, double t)
{
	std::string x;					//string with the name - output

	//check if it is 2 objective optimisation
	if (max_min_fit.size() == 2)
	{
		//set styles for the different line types
		x = "set style line 1 lt rgb \"red\" \n ";	//Style line for points
		x += "set style line 50 lc rgb \"green\" pt 7 ps 1 \n ";	//Style line for real PF legend
		x += "set style line 3 lc rgb \"green\" pt 6 ps 0.25 \n ";	//Style line for real PF
		x += "set style line 4 lt 1 lc rgb \"#FF1493\" \n ";	//Style line for points
		x += "set style line 5 lt 2 lc rgb \"#8B0000\" \n ";	//Style line for points
		x += "set style line 6 lt 5 lc rgb \"#FF8C00\" \n ";	//Style line for points
		x += "set style line 7 lt 13 lc rgb \"#9400D3\" \n ";	//Style line for points
		x += "set style line 8 lt 9 lc rgb \"#00FFFF\" \n ";	//Style line for points
		x += "set style line 9 lt 11 lc rgb \"#0000CD\" \n ";	//Style line for points
		x += "set style line 10 lt 20 lc rgb \"#006400\" \n ";	//Style line for points
		x += "set style line 11 lt 16 lc rgb \"#778899\" \n ";	//Style line for points
		
		//set legend position
		x += "set key outside\n";

		//set labels for x and y axis
		x += "set xlabel \"Fitness 1\"\nset ylabel \"Fitness 2\"\n";

		//set title
		x += "set title \" Generation #" + std::to_string(n) + "\"\n";
		
		//set the output to the png file
		x += "set term png size ";
		x += std::to_string(size.width);
		x += ",";
		x += std::to_string(size.height);
		
		//set location of the png file
		x += "\nset output \"Temp\\\\";
		x += input_file + "_" + std::to_string(n);
		
		//set x and y axis ranges
		x += ".jpg\"\nset xr[";
		x += std::to_string(max_min_fit[0].lower);
		x += ":";
		x += std::to_string(max_min_fit[0].upper);
		x += "] \nset yr[";
		x += std::to_string(max_min_fit[1].lower);
		x += ":";
		x += std::to_string(max_min_fit[1].upper);
		x += "] \n";
		
		//plot functions
		x += "plot ";

		//if one objective do not print PF
		if (ONE_OBJ_OVERRIDE != true)
		{
			x += "\"Temp\\\\";
			x += real_PF_out + "_" + std::to_string(indexr) + "_" + String_Prec(t, 1) + ".x1";
			x += "\" using 1:2 with points ls 3 notitle";
			x += ", 1/0 with points ls 50 title \"Real PF\", ";
		}

		//print plots for each collective
		for (short i = 0; i < n_col; i++)
		{
			x += "\"Temp\\\\";
			x += input_file + "_c" + std::to_string(i + 1) + "_" + std::to_string(n) + ".x1\"" ;
			x += " using 1:2 with points ls ";
			x += std::to_string(i + 4);
			x += " title \"Col #" + std::to_string(i + 1) + "\"";
			if (i != n_col - 1)
				x += ", ";
		}
		x += "\nset term wxt 1 \n";
	}
	//if not throw error
	else
	{
		std::cout << "ERROR#14: VIDEO - LABEL SIZE";
		system("pause");
		abort();
	}

	//assign the output string to the print name
	name_print = x;
}

/*
*Default real constructor*
@param pipename pipe for image/graph generator
@param ifile input file name template
@param ofile output file name template
@param framew width of the video frame - default value defined in Const.h
@param frameh height of the video frame - default value defined in Const.h
@param fpsn Frame per second in the video  - default value defined in Const.h
@param colour true - colour video, false - grayscale
*/
Video::Video(FILE * pipename, std::string ifile, std::string ofile, short nobj, short framew, short frameh, double fpsn, bool colour)
{
	//check if the video should be created
	if (VIDEO == false || FITNESS_ALL == false || nobj != 2)
		return;

	//set the video parameters
	size.width = framew;
	size.height = frameh;
	fps = fpsn;
	is_colour = colour;
	Pipe = pipename;
	input_file = ifile;
	output_file = ofile;
	frames = 0;

	int Fcc = video_v.fourcc('M', 'J', 'P', 'G');
	const char * temp_c = output_file.c_str();

	//open the video file for saving
	video_v.open(temp_c, Fcc, fps, size, is_colour);	
}

/*
*Making frame for video and saving it as .jpg*
@param n generation index number
@param indexr index of the current run
@param t time state of the current run
*/
void Video::Make_Frame(int n, int indexr, double t)
{
	if (n % frame_skip != 1 && (MODE == "PAES" || MODE == "MOEAD" || MODE == "MOEADMSF" || MODE == "MOEADPSF" || MODE == "MOEADM2M"))
		return;

	//check if the video should be created
	if (VIDEO == false || FITNESS_ALL == false)
		return;

	//get the string for the frame printing
	Name_Get(n,indexr, t);
	const char * temp_c = name_print.c_str();

	//print the frame
	fprintf(Pipe, temp_c);
	fflush(Pipe);

	//increase the number of frames
	frames++;
	
}

void Video::Make_Video()
{
	//check if the video should be created
	if (VIDEO == false || FITNESS_ALL == false)
		return;

	//add the frames to the movie
	for (long i = 0; i < frames; i++)
	{
		int frame_skip_temp = 1;
		if (MODE == "MOEAD" || MODE == "PAES" || MODE == "MOEADMSF" || MODE == "MOEADPSF" || MODE == "MOEADM2M")
			frame_skip_temp = frame_skip;
		//get the frame location
		std::string img_name = "Temp/";
		img_name += input_file;
		img_name += "_" + std::to_string(frame_skip_temp*i+1) + ".jpg";

		//read the image and make the proper frame from it
		cv::Mat frame = cv::imread(img_name);
		while (frame.rows != size.height || frame.cols != size.width)
		{
			int z = frame.type();
			frame = cv::Mat{ size.height, size.width, z };
			frame = cv::imread(img_name);
		}

		//add frames to the video
		for (long i = 0; i < frame_img_multi; i++)
			video_v << frame;
	}

	//release the video and save it to file
	video_v.release();
}


