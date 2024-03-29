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

  VIDEO header

Storage of the class used for video creation

*/

#pragma once
#ifndef VIDEO_H
#define VIDEO_H

#include <opencv2\videoio.hpp>
#include <string>
#include "Struct.h"
#include "Const.h"

/**class used for video creation*/
class Video
{
private:
	cv::VideoWriter video_v;								///<Video writer file from OpenCV
	cv::Size size;											///<Size of the movie frame
	bool is_colour;											///<Colour / Black and white; Default is True (colour)
	double fps;												///<FPS of the movie;
	int frames;												///<Number of frames in the movie
	FILE * Pipe;											///<Output pipe - for GNUPlot
	std::string name_print;									///<String for frame creation
	std::string input_file;									///<Input file name - data file
	std::string output_file;								///<Output file name - image file
	std::vector<STRUCTURES::boundaries> max_min_fit;		///<Boundaries for max and min fitness
protected:
	/**Get the name of input/outpit for the video frame generation.
	@param n - generation index
	@param indexr - index of the current run
	@param t - time step of the current run
	*/
	void Name_Get(int n, int indexr, double t);								
public:
	/**Default empty constructor*/
	Video();
	/**Default normal constructor
	@param pipename - pipe for image/graph generator
	@param ifile - input file name template
	@param ofile - output file name template
	@param framew - width of the video frame. Default value defined in Const.h
	@param frameh - height of the video frame. Default value defined in Const.h
	@param fpsn - Frame per second in the video. Default value defined in Const.h
	@param colour - true - colour video, false - grayscale.
	*/
	Video(FILE * pipename, std::string ifile, std::string ofile, short nobj, short framew = frame_w, short frameh = frame_h, double fpsn = FPS, bool colour = true);
	~Video() {}
	/**Setting max and min fitness for the video frame generation
	@param mmf - vector of the boundaries with max and min fitness values
	*/
	void Set_Max_Min(std::vector<STRUCTURES::boundaries> mmf) { max_min_fit = mmf; }
	/**Making frame for the video and saving it as .jpg
	@param n - generation index
	@param indexr - index of the current run
	@param t - time state of the current run
	*/
	void Make_Frame(int n, int indexr, double t);
	/**Making video from frames and saving to the file*/
	void Make_Video();							
};


#endif // !VIDEO_H
