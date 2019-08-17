/*--------------/
ALH
base.cpp
CMSASCPP: Construct, Merge, Solve & Adapt (CMSA) algorithm for the SCPP
Article: Evolutionary Methods for the Score-Constrained Packing Problem, A.L. Hawa, R. Lewis, J.M. Thompson
17/05/2019
/--------------*/

#include <iterator>
#include <iomanip>
#include <cmath>
#include <string>
#include "consts.h"
#include "base.h"

int instance; //Number of instances
int timeLimit = 3600; //How long EA will run for (seconds)
int tau = 70; //Minimum scoring distance
int numItem = 0; //Number of items n in the set I in the problem instance file
int stripWidth = 2500; //Width of strips
int numPop = 5; //Number of initial solutions in population
int maxAge = 5; //Max number of iterations of CMSA that a component can stay in BStar
int seed = 1;
int instType = 0;
int numIterations = 0; //Counts number of iterations of EA within the given time limit
std::ofstream outputStream;

void ProgramInfo() {

    std::cout << "CMSA - Construct, Merge, Solve and Adapt:\n-------------\n"
              << "INPUT FILE:\n"
              << "       .txt file must be first argument.\n"
              << "PARAMETERS:\n"
              << "       -e <int>    [Time limit for CMSA. Default = 3600s.]\n"
              << "       -t <int>    [Constraint value. Default = 70.]\n"
              << "       -W <int>    [Width of strips. Default = 2500.]\n"
              << "       -p <int>    [Number of solutions in population. Default = 5.]\n"
              << "       -a <int>    [Max age. Default = 5.]\n"
              << "       -s <int>    [Seed. Default = 1.]\n"
              << "       -xc <int>   [Time limit for Exact Cover. Default = 600s.]\n"
              << "---------------\n\n";
}

void ArgumentCheck() {

    bool error = false;

    if(tau == 0) {
        std::cerr << "[ERROR]: Constraint value cannot be zero.\n";
        error = true;
    }
    if(stripWidth == 0) {
        std::cerr << "[ERROR]: Strip cannot have width zero.\n";
        error = true;
    }
    if(numPop < 1) {
        std::cerr << "[ERROR]: Insufficient number of solutions in population.\n";
        error = true;
    }
    if(maxAge == 0){
        std::cerr << "[ERROR]: Max age cannot be zero.\n";
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

void OutputCMSA(int lowerBound) {

    int n = numItem / 100;
    int w = floor(stripWidth / 1000);

    std::string outputFileName;
    if(instType == 1) {
        outputFileName = "CA" + std::to_string(n) + std::to_string(w) + std::to_string(numPop) + std::to_string(maxAge) + ".txt";
    }
    else {
        outputFileName = "CR" + std::to_string(n) + std::to_string(w) + std::to_string(numPop) + std::to_string(maxAge) + ".txt";
    }
    outputStream.open(outputFileName.c_str(), std::ios::app);
    if(outputStream.fail()) {
        std::cerr << "[ERROR]: Cannot write to outputStream file." << std::endl;
        exit(1);
    }
    //outputStream << numItem << "\t" << instType << "\t" << stripWidth << "\t"  << numPop << "\t" << maxAge << "\t" << timeLimit << "\t" << xcTimeLimit << std::endl;

    double solnQuality = static_cast<double>(Sbsf.size()) / static_cast<double>(lowerBound);

    //outputStream << "Inst\tlb\t#S\tq\tFitness\t#it\n";
    outputStream << instance << "\t" << lowerBound << "\t" << Sbsf.size() << "\t" << solnQuality << "\t" << bsfFit << "\t" << numIterations << std::endl;

    outputStream.close();

}
