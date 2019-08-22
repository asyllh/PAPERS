/*--------------/
ALH
main.cpp
CMSASCPP: Construct, Merge, Solve & Adapt (CMSA) algorithm for the SCPP
Article: Evolutionary Methods for the Score-Constrained Packing Problem, A.L. Hawa, R. Lewis, J.M. Thompson
17/05/2019
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
        else if(strcmp("-a", argv[j]) == 0) { maxAge = atoi(argv[++j]); }
        else if(strcmp("-s", argv[j]) == 0) { seed = atoi(argv[++j]); }
        else if(strcmp("-xc", argv[j]) == 0) { xcTimeLimit = atoi(argv[++j]); }
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
    //endregion

    srand(seed);

    std::ifstream ifs(argv[1]);
    if(!ifs) {
        std::cerr << "[ERROR]: Cannot read from input file." << std::endl;
        exit(1);
    }
    CreateFromFile(ifs, numScores, totalItemWidth, allScores, partners, adjMatrix, itemWidths, itemNumbers);
    ifs.close();

    ArgumentCheck();

    int lowerBound = LowerBound(totalItemWidth);

    int n = numItem / 100;
    int w = floor(stripWidth / 1000);

    std::string timeLogFileName;
    if(instType == 1) {
        timeLogFileName = "CTA" + std::to_string(n) + std::to_string(w) + std::to_string(numPop) + std::to_string(maxAge) + "_" + std::to_string(instance) + ".txt";
    }
    else {
        timeLogFileName = "CTR" + std::to_string(n) + std::to_string(w) + std::to_string(numPop) + std::to_string(maxAge) +"_" + std::to_string(instance) + ".txt";
    }
    timeLogStream.open(timeLogFileName.c_str());
    if(timeLogStream.fail()) {
        std::cerr << "[ERROR]: Cannot write to timeLogStream file." << std::endl;
        exit(1);
    }
    timeLogStream << instance << "\t" << numItem << "\t" << instType << "\t" << stripWidth << "\t" << numPop
                  << "\t" << maxAge << "\t" << timeLimit << "\t" << xcTimeLimit << std::endl;
    timeLogStream << "It\t#S\tFitness\tTime" << std::endl;
    timeLogStream << "0\t0\t0.000000\t0.0" << std::endl;

    std::string xcLogFileName;
    if(instType == 1) {
        xcLogFileName = "XCA" + std::to_string(n) + std::to_string(w) + std::to_string(numPop) + std::to_string(maxAge) +"_" + std::to_string(instance) + ".txt";
    }
    else {
        xcLogFileName = "XCR" + std::to_string(n) + std::to_string(w) + std::to_string(numPop) + std::to_string(maxAge) +"_" + std::to_string(instance) + ".txt";
    }
    xcLogStream.open(xcLogFileName.c_str());
    if(xcLogStream.fail()) {
        std::cerr << "[ERROR]: Cannot write to xcLogStream file." << std::endl;
        exit(1);
    }
    xcLogStream << instance << "\t" << numItem << "\t" << instType << "\t" << stripWidth << "\t" << numPop
                  << "\t" << maxAge << "\t" << timeLimit << "\t" << xcTimeLimit << std::endl;
    xcLogStream << "It\t#B\tbL\tbFitness\t#soln\tTime" << std::endl;


    Timer timer;

    CMSA(timer, numScores, allScores, partners, adjMatrix, itemWidths, itemNumbers);

    OutputCMSA(lowerBound);

    timeLogStream.close();
    xcLogStream.close();


}//END INT MAIN
