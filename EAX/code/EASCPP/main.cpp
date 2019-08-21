/*--------------/
ALH
main.cpp
EASCPP: Evolutionary Algorithm (EA) for the SCPP
Article: Evolutionary Methods for the Score-Constrained Packing Problem, A.L. Hawa, R. Lewis, J.M. Thompson
26/03/2019
/--------------*/

#include <cmath>
#include <cstring>
#include "base.h"
#include "func.h"

int main(int argc, char** argv) {
    if(argc <= 1) {
        ProgramInfo();
        exit(1);
    }

    //region User Arguments
    for(int j = 1; j < argc; ++j) {
        if(strcmp("-e", argv[j]) == 0) { timeLimit = atoi(argv[++j]); }
        else if(strcmp("-t", argv[j]) == 0) { tau = atoi(argv[++j]); }
        else if(strcmp("-W", argv[j]) == 0) { stripWidth = atoi(argv[++j]); }
        else if(strcmp("-p", argv[j]) == 0) { numPop = atoi(argv[++j]); }
        else if(strcmp("-x", argv[j]) == 0) { xOver = atoi(argv[++j]); }
        else if(strcmp("-s", argv[j]) == 0) { seed = atoi(argv[++j]); }
    }
    //endregion

    //region Variables
    int numScores; //Number of score widths
    double totalItemWidth; //Sum of all n item widths in the set I
    std::vector<int> allScores; //Vector containing all score widths (size = numScores)
    std::vector<int> partners;
    std::vector<std::vector<int> > adjMatrix;
    std::vector<std::vector<int> > itemWidths;
    std::vector<std::vector<int> > itemNumbers;
    std::vector<std::vector<std::vector<int> > > population;//3D matrix containing solutions
    std::vector<std::vector<int> > populationSum; //Contains sum of each strip for every solution in the population
    //endregion

    srand(seed);

    std::ifstream ifs(argv[1]);
    if(!ifs) {
        std::cerr << "[ERROR]: Cannot read from file." << std::endl;
        exit(1);
    }
    CreateFromFile(ifs, numScores, totalItemWidth, allScores, partners, adjMatrix, itemWidths, itemNumbers);
    ifs.close();

    ArgumentCheck();

    int lowerBound = LowerBound(totalItemWidth);

    CreateInitPop(numScores, allScores, partners, adjMatrix, itemWidths, itemNumbers, populationSum, population);

    std::string filename1;
    int n = numItem / 100;
    int w = floor(stripWidth / 1000);
    if(instType == 1) {
        filename1 = "ETA" + std::to_string(n) + std::to_string(w) + std::to_string(xOver) + "_" + std::to_string(instance) + ".txt";
    }
    else {
        filename1 = "ETR" + std::to_string(n) + std::to_string(w) + std::to_string(xOver) + "_" + std::to_string(instance) + ".txt";
    }
    timeLogStream.open(filename1.c_str());
    if(timeLogStream.fail()) {
        std::cerr << "[ERROR]: Cannot write to timeLogStream file." << std::endl;
        exit(1);
    }
    timeLogStream << instance << "\t" << numItem << "\t" << instType << "\t" << stripWidth << "\t" << xOver << "\t" << lowerBound << std::endl;
    timeLogStream << "#S\tFit\tTime" << std::endl;
    timeLogStream << bestSize << "\t" << bestFit << "\t0.0" << std::endl;

    std::string filename2;
    if(instType == 1) {
        filename2 = "EA" + std::to_string(n) + std::to_string(w) + std::to_string(xOver) + ".txt";
    }
    else {
        filename2 = "ER" + std::to_string(n) + std::to_string(w) + std::to_string(xOver) + ".txt";
    }
    outputEAStream.open(filename2.c_str(), std::ios::app);
    if(outputEAStream.fail()) {
        std::cerr << "[ERROR]: Cannot write to outputNSStream file." << std::endl;
        exit(1);
    }

    Timer timer;

    while(timer.elapsed() < timeLimit) {
        EA(timer, numScores, allScores, partners, adjMatrix, itemWidths, itemNumbers, populationSum, population);
        ++numIterations;
    }

    timeLogStream.close();

    OutputEA(lowerBound);

    outputEAStream.close();


}//END INT MAIN
