/*--------------/
ALH
base.h
CMSASCPP: Construct, Merge, Solve & Adapt (CMSA) algorithm for the SCPP
Article: Evolutionary Methods for the Score-Constrained Packing Problem, A.L. Hawa, R. Lewis, J.M. Thompson
17/05/2019
/--------------*/

#ifndef BASE_H
#define BASE_H

#include <iostream>
#include <fstream>
#include <vector>


void ProgramInfo();

void ArgumentCheck();

void CreateFromFile(std::ifstream& ifs, int& numScores, double& totalItemWidth, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
                    std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers);

void OutputCMSA(int lowerBound);


#endif


