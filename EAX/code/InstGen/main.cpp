/*--------------/
ALH
main.cpp
InstGen: Problem Instance Generator for the SCPP
Article: Evolutionary Methods for the Score-Constrained Packing Problem, A.L. Hawa, R. Lewis, J.M. Thompson
11/01/2019
/--------------*/

#include <iostream>
#include <cstring>
#include "gfunc.h"

int main(int argc, char** argv) {
    if(argc <= 1) {
        ProgramInfo();
        exit(1);
    }


    //Reading in user arguments
    for(int j = 1; j < argc; ++j) {
        if(strcmp("-i", argv[j]) == 0) { instance = atoi(argv[++j]); }
        else if(strcmp("-t", argv[j]) == 0) { tau = atoi(argv[++j]); }
        else if(strcmp("-n", argv[j]) == 0) { numItem = atoi(argv[++j]); }
        else if(strcmp("-a", argv[j]) == 0) { minWidth = atoi(argv[++j]); }
        else if(strcmp("-b", argv[j]) == 0) { maxWidth = atoi(argv[++j]); }
        else if(strcmp("-m", argv[j]) == 0) { minItemWidth = atoi(argv[++j]); }
        else if(strcmp("-M", argv[j]) == 0) { maxItemWidth = atoi(argv[++j]); }
        else if(strcmp("-s", argv[j]) == 0) { seed = atoi(argv[++j]); }
        else if(strcmp("-c", argv[j]) == 0) { instType = atoi(argv[++j]); }
    }


    int numScores = numItem * 2; //Number of score widths
    double totalItemWidth; //Sum of all n item widths in the set I
    std::vector<int> allScores; //Vector containing all score widths (size = numScores)
    std::vector<int> partners(numScores, 0); //Vector containing partners
    std::vector<std::vector<int> > itemWidths(numScores, std::vector<int>(numScores, 0)); //Matrix containing item widths
    std::vector<std::vector<int> > adjMatrix(numScores, std::vector<int>(numScores, 0)); //Matrix showing score widths that meet VSC or that are partners
    std::vector<std::vector<int> > itemNumbers(numScores, std::vector<int>(numScores, 0)); //Matrix numbering all items from 0 to numItem - 1 (allItems)
    std::vector<std::vector<int> > itemTypes(4);
    std::vector<int> typeNumber(numScores);

    srand(seed);

    Timer timer;

    if(instType == 1){
        CreateAInstance(numScores, totalItemWidth, allScores, partners, adjMatrix, itemWidths, itemNumbers);
    }
    else {
        CreateRInstance(numScores, totalItemWidth, allScores, partners, adjMatrix, itemWidths, itemNumbers, itemTypes, typeNumber);
    }

    OutputProbInst(numScores, totalItemWidth, allScores, partners, itemWidths, itemNumbers, itemTypes, typeNumber);


}// END INT MAIN