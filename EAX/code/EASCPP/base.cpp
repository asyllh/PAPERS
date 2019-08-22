/*--------------/
ALH
base.cpp
EASCPP: Evolutionary Algorithm (EA) for the SCPP
Article: Evolutionary Methods for the Score-Constrained Packing Problem, A.L. Hawa, R. Lewis, J.M. Thompson
26/03/2019
/--------------*/

#include <iterator>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <string>
#include "consts.h"
#include "base.h"

int instance; //Number of instances
int timeLimit = 600; //How long EA will run for (seconds)
int tau = 70; //Minimum scoring distance
int numItem = 0; //Number of items n in the set I in the problem instance file
int stripWidth = 2500; //Width of strips
int numPop = 25; //Number of initial solutions in population
int xOver = 1; //Crossover Type, 1 = GGA, 2 = AGX, 3 = AGX'
int seed = 1;
int instType = 0;
int numIterations = 0; //Counts number of iterations of EA within the given time limit
std::ofstream outputEAStream;

void ProgramInfo() {

    std::cout << "EA - Evolutionary Algorithm for the SCSPP:\n-------------\n"
              << "INPUT FILE:\n"
              << "       .txt file must be first argument.\n"
              << "PARAMETERS:\n"
              << "       -e <int>    [Time limit for EA. Default = 600s.]\n"
              << "       -t <int>    [Constraint value. Default = 70.]\n"
              << "       -W <int>    [Width of strips. Default = 2500.]\n"
              << "       -p <int>    [Number of solutions in population. Default = 25.]\n"
              << "       -x <int>    [Crossover operator. 1: GGA. 2: AGX. 3: AGX'. Default = GGA.]\n"
              << "       -s <int>    [Seed. Default = 1.]\n"
              << "---------------\n\n";
}

void ArgumentCheck() {

    bool error = false;

    std::cout << "EA - Evolutionary Algorithm for the SCSPP\n------------------------------\n";
    if(timeLimit > 0 && timeLimit < 60) {
        std::cout << "[WARNING]: Time limit is set to " << timeLimit << " seconds, potentially insufficient amount of time for the program.\n";
    }
    if(timeLimit == 0) {
        std::cerr << "[ERROR]: Time limit cannot be zero.\n";
        error = true;
    }
    if(tau == 0) {
        std::cout << "[WARNING]: Constraint value is zero, problem is equivalent to BPP.\n";
    }
    if(stripWidth == 0) {
        std::cerr << "[ERROR]: Strip cannot have width zero.\n";
        error = true;
    }
    if(numPop < 2) {
        std::cerr << "[ERROR]: Insufficient number of solutions in population.\n";
        error = true;
    }
    if(xOver != 1 && xOver != 2 && xOver != 3) {
        std::cerr << "[ERROR]: Invalid choice of recombination operator. Please choose either 1: GGA, 2: AGX, or 3: AGX'.\n";
        error = true;
    }
    if(instType != 1 && instType != 2) {
        std::cerr << "[ERROR]: Invalid instance type in file. Please choose either 1: Artificial, or 2: Real.\n";
        error = true;
    }

    if(error) {
        std::cerr << "[EXIT PROGRAM.]\n";
        exit(1);
    }

}

void CreateFromFile(std::ifstream& ifs, int& numScores, double& totalItemWidth, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
                    std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers) {


    ifs >> instance;
    ifs >> numItem;
    ifs >> instType;
    numScores = numItem * 2;

    partners.resize(numScores, 0);

    adjMatrix.resize(numScores, std::vector<int>(numScores, 0));

    itemNumbers.resize(numScores, std::vector<int>(numScores, 0));

    itemWidths.resize(numScores, std::vector<int>(numScores, 0));

    for(int i = 0; i < numScores; ++i) {
        int a;
        ifs >> a;
        allScores.push_back(a);
    }

    for(int i = 0; i < numScores; ++i) {
        int p;
        ifs >> p;
        partners[i] = p;
    }


    for(int i = 0; i < numScores; ++i) {
        int w;
        ifs >> w;
        itemWidths[i][partners[i]] = w;
    }

    for(int i = 0; i < numScores; ++i) {
        int n;
        ifs >> n;
        itemNumbers[i][partners[i]] = n;
    }

    ifs >> totalItemWidth;

    //end of reading in from file

    for(int i = 0; i < allScores.size() - 1; ++i) {
        for(int j = i + 1; j < allScores.size(); ++j) {
            if(allScores[i] + allScores[j] >= tau) {
                adjMatrix[i][j] = 1;
                adjMatrix[j][i] = 1;
            }
        }
    }

    for(int i = 0; i < partners.size(); ++i) {
        adjMatrix[i][partners[i]] = 2;
    }


}

void OutputEA(int lowerBound) {

    double solnQuality = static_cast<double>(bestSize) / static_cast<double>(lowerBound);
    double propFeas = static_cast<double>(numFeas) / static_cast<double>(numSubSCP);

    outputEAStream << instance << "\t" << lowerBound << "\t" << bestSize << "\t" << solnQuality << "\t" << bestFit << "\t" << numIterations << "\t" << propFeas << "\t" << bestTime << std::endl;

}
