#pragma once
#ifndef SVM_H
#define SVM_H

#include <vector>
#include "Class.h"


/*
*SVM function for given population - generate a vector of labels for separation*
@param pop population for which SVM will calculate labels
*/
std::vector<short> SVM(const population & pop);
/*
*SVM training and model generation*
@param pop population which will be trained
*/
static void SVM_TRAIN(const population & pop);
/*
*SVM Real label creation*
@param pop population for which label will be calculated
*/
static std::vector<short> SVM_PREDICT(const population & pop);	
/**SVM Setting up  SVM parameters**/
static void SVM_PARAM();
/*
*Generation of problem for SVM_TRAIN*
@param pop population for which problem will be separated
*/
static void SVM_PROB(const population & pop);
/*
*Generation of problem for SVM_PREDICT*
@param pop population for which problem will be separated
*/
static void SVM_PROB2(const population & pop);		
/*
*Generation of the train_train label*
@param pop population for which the label will be generated
*/
std::vector<int> SVM_LABEL(const population & pop);		
#endif // !SVM_H

