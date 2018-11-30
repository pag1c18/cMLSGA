#include "Random_Label.h"
#include "Random.h"
/*
*Generate random labels for the population*
@param pop source population for which label will be found
*/

extern short n_col;

std::vector<short> Random_Label(const population & pop)
{
	int psize = pop.Size_Show();		//size of the population
	std::vector<short> label (psize,0);			//vector for labels storage - output

	//assign labels by random
	for (int i = 0; i < psize; i++)
	{
		label[i] = Random_I(1, n_col);
	}
	
	//return the vector
	return label;
}