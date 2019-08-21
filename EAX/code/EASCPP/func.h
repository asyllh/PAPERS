/*--------------/
ALH
func.h
EASCPP: Evolutionary Algorithm (EA) for the SCPP
Article: Evolutionary Methods for the Score-Constrained Packing Problem, A.L. Hawa, R. Lewis, J.M. Thompson
26/03/2019
/--------------*/

#ifndef FUNC_H
#define FUNC_H

#include <string>
#include "consts.h"


void Swap(int& a, int& b);

int LowerBound(double totalItemWidth);

double Fitness(std::vector<int>& stripSum, std::vector<std::vector<int> >& strip);

void PermuteStrips(std::vector<std::vector<int> >& strip, std::vector<int>& stripSum);

void CreateInitPop(int numScores, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
                   std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<std::vector<int> >& populationSum,
                   std::vector<std::vector<std::vector<int> > >& population);

void MFFPlus(int numScores, int order, std::vector<int>& allScores, std::vector<int>& partners, std::vector<int>& partialItem,
             std::vector<std::vector<int> >& adjMatrix, std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers,
             std::vector<int>& stripSum, std::vector<std::vector<int> >& strip);

void PreAHCMFFP(int i1, int j1, int& feasible, std::vector<int>& allScores, std::vector<int>& partners, std::vector<int>& itemOrder,
                std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<int>& stripSum,
                std::vector<std::vector<int> >& strip);

void Mutation(int numScores, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
              std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<int>& stripSum,
              std::vector<std::vector<int> >& strip);

void LocalSearch(int numScores, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
                 std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<int>& stripSum,
                 std::vector<std::vector<int> >& strip, std::vector<int>& stripSumX, std::vector<std::vector<int> >& stripX,
                 std::vector<int>& stripSumY, std::vector<std::vector<int> >& stripY);

void PairPair(std::vector<int>& allScores, std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers,
              std::vector<int>& stripSumX, std::vector<std::vector<int> >& stripX, std::vector<int>& stripSumY, std::vector<std::vector<int> >& stripY);

void PairSin(std::vector<int>& allScores, std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers,
             std::vector<int>& stripSumX, std::vector<std::vector<int> >& stripX, std::vector<int>& stripSumY, std::vector<std::vector<int> >& stripY);

void SinSin(std::vector<int>& allScores, std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers,
            std::vector<int>& stripSumX, std::vector<std::vector<int> >& stripX, std::vector<int>& stripSumY, std::vector<std::vector<int> >& stripY);

void MoveSin(int& feasible, std::vector<int>& allScores, std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers,
             std::vector<int>& stripSumX, std::vector<std::vector<int> >& stripX, std::vector<int>& stripSumY, std::vector<std::vector<int> >& stripY);

void PreAHCLS(int i1, int a1, int b1, int j1, int c1, int d1, int& feasible, int swapType, int moveType, std::vector<int>& allScores,
              std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<int>& stripSumX,
              std::vector<std::vector<int> >& stripX, std::vector<int>& stripSumY, std::vector<std::vector<int> >& stripY);

void AHC(int& feasible, std::vector<int>& scores, std::vector<int>& original, std::vector<int>& final);

void InitInstance(int nScores, std::vector<std::vector<int> >& adjMat, std::vector<int>& scores, std::vector<int>& order, std::vector<int>& B);

void MCM(int nScores, int& matchSize, std::vector<std::vector<int> >& adjMat, std::vector<int>& B, std::vector<int>& matchList,
         std::vector<int>& cycleVertex);

void MPS(int nScores, int& nCycles, std::vector<int>& B, std::vector<int>& matchList, std::vector<std::vector<int> >& mpStructure);

void BCR(int nScores, int& feasible, int matchSize, int nCycles, std::vector<int>& B, std::vector<int>& matchList,
         std::vector<int>& cycleVertex, std::vector<std::vector<int> >& mpStructure, std::vector<std::vector<int> >& adjMat, std::vector<int>& altHam);

void CP(int nScores, int full, bool multiple, std::vector<int>& matchList, std::vector<int>& B, std::vector<int>& altHam, std::vector<std::vector<int> >& RStar);

void EA(const Timer& timer, int numScores, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
        std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<std::vector<int> >& populationSum,
        std::vector<std::vector<std::vector<int> > >& population);

void GGA(int numScores, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
         std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<int>& offspringSum,
         std::vector<std::vector<int> >& offspring, std::vector<int>& stripSumX, std::vector<std::vector<int> >& stripX,
         std::vector<int>& stripSumY, std::vector<std::vector<int> >& stripY);

void AGX(int numScores, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
         std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<int>& offspringSum,
         std::vector<std::vector<int> >& offspring, std::vector<int>& stripSumX, std::vector<std::vector<int> >& stripX,
         std::vector<int>& stripSumY, std::vector<std::vector<int> >& stripY);

void AGXDash(int numScores, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
             std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<int>& offspringSum,
             std::vector<std::vector<int> >& offspring, std::vector<int>& stripSumX, std::vector<std::vector<int> >& stripX,
             std::vector<int>& stripSumY, std::vector<std::vector<int> >& stripY);


#endif
