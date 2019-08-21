/*--------------/
ALH
consts.h
EASCPP: Evolutionary Algorithm (EA) for the SCPP
Article: Evolutionary Methods for the Score-Constrained Packing Problem, A.L. Hawa, R. Lewis, J.M. Thompson
26/03/2019
/--------------*/

#ifndef CONSTS_H
#define CONSTS_H

#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <set>


class Timer {
private:
    std::chrono::high_resolution_clock::time_point m_start;

public:
    Timer()
            : m_start(std::chrono::high_resolution_clock::now())
    {

    }

    void reset()
    {
        m_start = std::chrono::high_resolution_clock::now();
    }

    double elapsed() const
    {
        return std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - m_start).count();
    }

    /*~Timer()
    {
        std::cout << "Total Time: "
                  << std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - m_start).count()
                << "seconds." << std::endl;
    }*/
};


extern int instance;
extern int timeLimit;
extern int tau;
extern int numItem;
extern int stripWidth;
extern int numPop;
extern int xOver;
extern int seed;
extern int instType;
extern int numIterations;
extern std::ofstream outputEAStream;

extern int bestSize;
extern double bestFit;
extern long numFeas;
extern long numInfeas;
extern long numSubSCP;
extern std::ofstream timeLogStream; //main.cpp, eafunc.cpp
extern double bestTime;

#endif
