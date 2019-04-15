/*
Copyright (c) 2000-2018 Chih-Chung Chang and Chih-Jen Lin
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither name of copyright holders nor the names of its contributors
may be used to endorse or promote products derived from this software
without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/


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

