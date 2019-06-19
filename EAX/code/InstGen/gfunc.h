/*--------------/
ALH
gfunc.h
InstGen: Problem Instance Generator for the SCPP
Article: Evolutionary Methods for the Score-Constrained Packing Problem, A.L. Hawa, R. Lewis, J.M. Thompson
11/01/2019
/--------------*/

#ifndef G_FUNC_H
#define G_FUNC_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <string>
#include <random>

class Timer {
private:
    std::chrono::high_resolution_clock::time_point startTime, endTime;
    std::chrono::duration<float> totalTime;

public:
    Timer() {
        startTime = std::chrono::high_resolution_clock::now();
    }

    ~Timer() {
        endTime = std::chrono::high_resolution_clock::now();
        totalTime = endTime - startTime;
        float totalTimeMs = totalTime.count() * 1000.0f;
        float totalTimeMin = totalTime.count() / 60.0f;
        std::cout << "\nCPU Time: " << totalTimeMs << "ms -- " << totalTime.count() << "s -- " << totalTimeMin << "mins" << std::endl;
    }
};

extern int instance;
extern int tau;
extern int numItem;
extern int minWidth;
extern int maxWidth;
extern int minItemWidth;
extern int maxItemWidth;
extern int seed;
extern int instType;
extern int numTypes;

void ProgramInfo();


void CreateAInstance(int numScores, double& totalItemWidth, std::vector<int>& allScores, std::vector<int>& partners,
                     std::vector<std::vector<int> >& adjMatrix, std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers);

void CreateRInstance(int numScores, double& totalItemWidth, std::vector<int>& allScores, std::vector<int>& partners,
                     std::vector<std::vector<int> >& adjMatrix, std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers,
                     std::vector<std::vector<int> >& itemTypes, std::vector<int>& typeNumber);

void OutputProbInst(int numScores, double totalItemWidth, std::vector<int>& allScores, std::vector<int>& partners,
                    std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers,
                    std::vector<std::vector<int> >& itemTypes, std::vector<int>& typeNumber);


#endif
