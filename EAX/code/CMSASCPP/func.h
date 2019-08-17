/*--------------/
ALH
func.h
CMSASCPP: Construct, Merge, Solve & Adapt (CMSA) algorithm for the SCPP
Article: Evolutionary Methods for the Score-Constrained Packing Problem, A.L. Hawa, R. Lewis, J.M. Thompson
17/05/2019
/--------------*/

#ifndef FUNC_H
#define FUNC_H

#include <string>
#include "consts.h"
#include "dlx.h"


int LowerBound(double totalItemWidth);

double Fitness(std::vector<std::vector<int> >& strip, std::vector<int>& stripSum);

void Construct(int numScores, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
                std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers);

void Adapt();

void MFFPlus(int numScores, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
            std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<int>& stripSum, std::vector<std::vector<int> >& strip);

void PreAHCMFFP(int i1, int j1, int& feasible, std::vector<int>& allScores, std::vector<int>& partners, std::vector<int>& itemOrder,
                std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<int>& stripSum,
                std::vector<std::vector<int> >& strip);

void AHC(int& feasible, std::vector<int>& scores, std::vector<int>& original, std::vector<int>& final);

void InitInstance(int nScores, std::vector<std::vector<int> >& adjMat, std::vector<int>& scores, std::vector<int>& order, std::vector<int>& B);

void MCM(int nScores, int& matchSize, std::vector<std::vector<int> >& adjMat, std::vector<int>& B, std::vector<int>& matchList,
         std::vector<int>& cycleVertex);

void MPS(int nScores, int& nCycles, std::vector<int>& B, std::vector<int>& matchList, std::vector<std::vector<int> >& mpStructure);

void BCR(int nScores, int& feasible, int matchSize, int nCycles, std::vector<int>& B, std::vector<int>& matchList,
         std::vector<int>& cycleVertex, std::vector<std::vector<int> >& mpStructure, std::vector<std::vector<int> >& adjMat, std::vector<int>& altHam);

void CP(int nScores, int full, bool multiple, std::vector<int>& matchList, std::vector<int>& B, std::vector<int>& altHam, std::vector<std::vector<int> >& RStar);

void CMSA(const Timer& timer, int numScores, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
          std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers);


#endif
