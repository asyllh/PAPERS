/*--------------/
ALH
func.cpp
EASCPP: Evolutionary Algorithm (EA) for the SCPP
Article: Evolutionary Methods for the Score-Constrained Packing Problem, A.L. Hawa, R. Lewis, J.M. Thompson
26/03/2019
/--------------*/

#include <algorithm>
#include <cmath>
#include <climits>
#include "func.h"

int bestSize = 0;
double bestFit = 0.0;
long numFeas = 0;
long numInfeas = 0;
long numSubSCP = 0;
std::ofstream timeLogStream; //stream to open file to log timestamps when better soln found
double bestTime = 0.0;

void Swap(int& a, int& b) {
    int temp = a;
    a = b;
    b = temp;
}

int LowerBound(double totalItemWidth) {
    int lBound = ceil(totalItemWidth / stripWidth);
    return lBound;
}

double Fitness(std::vector<int>& stripSum, std::vector<std::vector<int> >& strip) {
    double total = 0.0;

    for(int i = 0; i < strip.size(); ++i) {
        double a = static_cast<double>(stripSum[i]) / static_cast<double>(stripWidth);
        total += pow(a, 2);
    }
    double final = total / static_cast<double>(strip.size());

    return final;
}

void PermuteStrips(std::vector<std::vector<int> >& strip, std::vector<int>& stripSum){

    for(int i = strip.size()-1; i >= 0; i--){
        int r = rand() % (i+1);
        strip[i].swap(strip[r]);
        Swap(stripSum[i], stripSum[r]);
    }
}

void CreateInitPop(int numScores, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
                   std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<std::vector<int> >& populationSum,
                   std::vector<std::vector<std::vector<int> > >& population) {

    std::vector<std::vector<int> > strip(numItem);
    std::vector<int> stripSum(numItem, 0);
    std::vector<int> temp; //throwaway vector, not used

    MFFPlus(numScores, 1, allScores, partners, temp, adjMatrix, itemWidths, itemNumbers, stripSum, strip);

    Mutation(numScores, allScores, partners, adjMatrix, itemWidths, itemNumbers, stripSum, strip);

    double tempFit = Fitness(stripSum, strip);

    population.push_back(strip);
    populationSum.push_back(stripSum);

    bestSize = population[0].size();
    bestFit = tempFit;


    for(int i = 1; i < numPop; ++i) {
        strip.clear();
        strip.resize(numItem);
        stripSum.clear();
        stripSum.resize(numItem, 0);

        MFFPlus(numScores, 2, allScores, partners, temp, adjMatrix, itemWidths, itemNumbers, stripSum, strip);

        Mutation(numScores, allScores, partners, adjMatrix, itemWidths, itemNumbers, stripSum, strip);

        population.push_back(strip);
        populationSum.push_back(stripSum);

        tempFit = Fitness(stripSum, strip);
        if(population[i].size() < bestSize){
            bestSize = population[i].size();
            bestFit = tempFit;

        }
        else if(population[i].size() == bestSize){
            if(tempFit > bestFit) {
                bestSize = population[i].size();
                bestFit = tempFit;

            }
        }

    }

}

void MFFPlus(int numScores, int order, std::vector<int>& allScores, std::vector<int>& partners, std::vector<int>& partialItem,
             std::vector<std::vector<int> >& adjMatrix, std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers,
             std::vector<int>& stripSum, std::vector<std::vector<int> >& strip) {

    //Modified First-Fit Plus (shell)

    int feasible;
    std::vector<int> itemOrder;
    std::vector<int> checked(numScores, 0);

    if(order == 1) { //Sort all items in DECREASING order (non-increasing), only used once when creating initial population
        int mini;
        int min = 0;
        int max = 1000; //1000 = maxItemWidth
        while(itemOrder.size() < numItem) {
            for(int i = 0; i < numScores; ++i) {
                if(checked[i] == 1) {
                    continue;
                }
                if(itemWidths[i][partners[i]] > min && itemWidths[i][partners[i]] <= max) {
                    min = itemWidths[i][partners[i]];
                    mini = i;
                }
            }
            itemOrder.push_back(mini);
            checked[mini] = 1;
            checked[partners[mini]] = 1;
            max = min;
            min = 0;
        }
    }
    else if(order == 2) { //Sort all items in RANDOM order, only used when creating initial population
        while(itemOrder.size() < numItem) {
            int r = rand() % numScores;
            while(checked[r] == 1) {
                r = rand() % numScores;
            }
            if(r < partners[r]) {
                itemOrder.push_back(r);
                checked[r] = 1;
                checked[partners[r]] = 1;
            }
            else {
                itemOrder.push_back(partners[r]);
                checked[r] = 1;
                checked[partners[r]] = 1;
            }
        }
    }
    else if(order == 3) { //Sort the items in partialItem vector in DECREASING order (non-increasing), partialItem vector has indices of score widths
        //Used during EA, after local search and xOver operator
        int mini;
        int min = 0;
        int max = 1000; //1000 = maxItemWidth
        while(itemOrder.size() < partialItem.size() / 2) {
            for(int i = 0; i < partialItem.size(); ++i) {
                if(checked[partialItem[i]] == 1) {
                    continue;
                }
                if(itemWidths[partialItem[i]][partners[partialItem[i]]] > min
                   && itemWidths[partialItem[i]][partners[partialItem[i]]] <= max) {
                    min = itemWidths[partialItem[i]][partners[partialItem[i]]];
                    mini = partialItem[i];
                }
            }
            itemOrder.push_back(mini);
            checked[mini] = 1;
            checked[partners[mini]] = 1;
            max = min;
            min = 0;
        }
    }
    else {
        std::cerr << "[ERROR]: No order specified for MFFPlus." << std::endl;
        exit(1);
    }

    strip[0].push_back(itemOrder[0]);
    strip[0].push_back(partners[itemOrder[0]]);
    stripSum[0] += itemWidths[itemOrder[0]][partners[itemOrder[0]]];

    for(int j = 1; j < itemOrder.size(); ++j) {
        for(int i = 0; i < strip.size(); ++i) {
            if(!strip[i].empty()) {
                if(stripSum[i] + itemWidths[itemOrder[j]][partners[itemOrder[j]]] <= stripWidth) {
                    if(strip[i].size() == 2) { //If the strip only contains one item, don't run AHC, just do checks instead
                        if(adjMatrix[strip[i].back()][itemOrder[j]] == 1) {
                            strip[i].push_back(itemOrder[j]);
                            strip[i].push_back(partners[itemOrder[j]]);
                            stripSum[i] += itemWidths[itemOrder[j]][partners[itemOrder[j]]];
                            break;
                        }
                        else if(adjMatrix[strip[i].back()][partners[itemOrder[j]]] == 1) {
                            strip[i].push_back(partners[itemOrder[j]]);
                            strip[i].push_back(itemOrder[j]);
                            stripSum[i] += itemWidths[itemOrder[j]][partners[itemOrder[j]]];
                            break;
                        }
                        else if(adjMatrix[strip[i].front()][itemOrder[j]] == 1) {
                            strip[i].insert(strip[i].begin(), itemOrder[j]);
                            strip[i].insert(strip[i].begin(), partners[itemOrder[j]]);
                            stripSum[i] += itemWidths[itemOrder[j]][partners[itemOrder[j]]];
                            break;
                        }
                        else if(adjMatrix[strip[i].front()][partners[itemOrder[j]]] == 1) {
                            strip[i].insert(strip[i].begin(), partners[itemOrder[j]]);
                            strip[i].insert(strip[i].begin(), itemOrder[j]);
                            stripSum[i] += itemWidths[itemOrder[j]][partners[itemOrder[j]]];
                            break;
                        }
                    }
                    else { // If more than one item on the strip, run AHC
                        feasible = 0;
                        PreAHCMFFP(i, j, feasible, allScores, partners, itemOrder, itemWidths, itemNumbers, stripSum, strip);
                        if(feasible == 1) { //item has been packed
                            break;
                        }
                    }
                }
            }
            else if(strip[i].empty()) {
                strip[i].push_back(itemOrder[j]);
                strip[i].push_back(partners[itemOrder[j]]);
                stripSum[i] += itemWidths[itemOrder[j]][partners[itemOrder[j]]];
                break;
            }
        }
    }

    while(stripSum.back() == 0) {
        stripSum.pop_back();
        strip.pop_back();
    }

}

void PreAHCMFFP(int i1, int j1, int& feasible, std::vector<int>& allScores, std::vector<int>& partners, std::vector<int>& itemOrder,
                std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<int>& stripSum,
                std::vector<std::vector<int> >& strip) {

    feasible = 0;
    std::vector<int> scores, original, final, checkScores;

    for(int k = 0; k < strip[i1].size(); ++k) {
        scores.push_back(allScores[strip[i1][k]]);
        original.push_back(strip[i1][k]);
    }
    scores.push_back(allScores[itemOrder[j1]]);
    original.push_back(itemOrder[j1]);
    scores.push_back(allScores[partners[itemOrder[j1]]]);
    original.push_back(partners[itemOrder[j1]]);

    checkScores = scores;
    std::sort(checkScores.begin(), checkScores.end());
    if(checkScores[2] + checkScores[checkScores.size() - 1] < tau) { //failed prelim checks, no need to run AHC
        feasible = 0;
        return;
    }
    scores.push_back(tau);
    scores.push_back(tau);

    AHC(feasible, scores, original, final);

    if(feasible == 1) {
        stripSum[i1] += itemWidths[itemOrder[j1]][partners[itemOrder[j1]]];
        strip[i1] = final;
        return;
    }
    else if(feasible == 0) {
        return;
    }
}

void Mutation(int numScores, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
              std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<int>& stripSum,
              std::vector<std::vector<int> >& strip) {

    std::vector<int> stripSumX;
    std::vector<std::vector<int> > stripX;
    std::vector<int> stripSumY;
    std::vector<std::vector<int> > stripY;

    PermuteStrips(strip, stripSum);


    int r = rand() % (strip.size() - 1) + 1;

    for(int i = 0; i < r; ++i) {
        stripX.push_back(strip[i]);
        stripSumX.push_back(stripSum[i]);
    }

    for(int i = r; i < strip.size(); ++i) {
        stripY.push_back(strip[i]);
        stripSumY.push_back(stripSum[i]);
    }

    strip.clear();
    stripSum.clear();

    LocalSearch(numScores, allScores, partners, adjMatrix, itemWidths, itemNumbers, stripSum, strip, stripSumX, stripX, stripSumY, stripY);

}

void LocalSearch(int numScores, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
                 std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<int>& stripSum,
                 std::vector<std::vector<int> >& strip, std::vector<int>& stripSumX, std::vector<std::vector<int> >& stripX,
                 std::vector<int>& stripSumY, std::vector<std::vector<int> >& stripY) {

    int feasible = 0;

    PermuteStrips(stripX, stripSumX);
    PermuteStrips(stripY, stripSumY);

    do {
        feasible = 0;

        PairPair(allScores, itemWidths, itemNumbers, stripSumX, stripX, stripSumY, stripY);

        PairSin(allScores, itemWidths, itemNumbers, stripSumX, stripX, stripSumY, stripY);

        SinSin(allScores, itemWidths, itemNumbers, stripSumX, stripX, stripSumY, stripY);

        MoveSin(feasible, allScores, itemWidths, itemNumbers, stripSumX, stripX, stripSumY, stripY);
    } while(feasible == 1);

    //region End
    if(feasible == 2) {
        strip = stripX;
        stripSum = stripSumX;
    }
    else { //feasible == 0 at the end of MoveSin
        //Do FFD on stripY
        std::vector<int> partialItem;
        for(int i = 0; i < stripY.size(); ++i) {
            for(int j = 0; j < stripY[i].size(); ++j) {
                partialItem.push_back(stripY[i][j]);
            }
        }
        std::sort(partialItem.begin(), partialItem.end());

        stripY.clear();
        stripY.resize(partialItem.size() / 2);
        stripSumY.clear();
        stripSumY.resize(partialItem.size() / 2, 0);

        MFFPlus(numScores, 3, allScores, partners, partialItem, adjMatrix, itemWidths, itemNumbers, stripSumY, stripY);

        //join sets stripX and stripY together back into vector<vector<int> > strip
        strip.clear();
        stripSum.clear();

        for(int i = 0; i < stripX.size(); ++i) {
            strip.push_back(stripX[i]);
            stripSum.push_back(stripSumX[i]);
        }
        for(int i = 0; i < stripY.size(); ++i) {
            strip.push_back(stripY[i]);
            stripSum.push_back(stripSumY[i]);
        }
    }
    //endregion

}

void PairPair(std::vector<int>& allScores, std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers,
              std::vector<int>& stripSumX, std::vector<std::vector<int> >& stripX, std::vector<int>& stripSumY, std::vector<std::vector<int> >& stripY) {

    int swapType, moveType, feasible = 0;

    /*SWAPPING A PAIR OF BOXES FROM EACH SET*/
    for(int i = 0; i < stripX.size(); ++i) { //For each strip in the set stripX
        if(stripX[i].size() >= 4) { //If there are at least 2 boxes on stripX[i] (note that each element represents a score, so 4 elements = 2 boxes)
            //Go through each pair of boxes on stripX[i]
            for(int a = 0; a < stripX[i].size() - 3; a += 2) { //Starting from the first score on the first box until the first score on the penultimate box
                for(int b = a + 2; b < stripX[i].size() - 1; b += 2) { //Starting from the first score on the second box until the first score on the last box
                    int pairSizeX = itemWidths[stripX[i][a]][stripX[i][a + 1]] + itemWidths[stripX[i][b]][stripX[i][b + 1]]; //Sum box widths
                    //Check if there exists a pair of boxes on a strip in set stripY that have a combined width larger than pairSizeX
                    for(int j = 0; j < stripY.size(); ++j) { //For each strip in the set stripY
                        if(stripY[j].size() >= 4) { //If there are at least 2 boxes on stripY[j]
                            //Go through each pair of boxes on stripY[j]
                            for(int c = 0; c < stripY[j].size() - 3; c += 2) { //Starting from the first score on the first box until the first score on the penultimate box
                                for(int d = c + 2; d < stripY[j].size() - 1; d += 2) { //Starting from the first score on the second box until the first score on the last box
                                    int pairSizeY = itemWidths[stripY[j][c]][stripY[j][c + 1]] + itemWidths[stripY[j][d]][stripY[j][d + 1]]; //Sum box widths
                                    //Check if pairSizeX < pairSizeY and that boxes can fit onto strip
                                    if(pairSizeX < pairSizeY && stripSumX[i] - pairSizeX + pairSizeY <= stripWidth) {
                                        swapType = 1;
                                        //cout << "i: " << i << " j: " << j << " a: " << a << " b: " << b << " c: " << c << " d: " << d << endl;
                                        if(stripX[i].size() == 4) { //If stripX[i] only contains 2 boxes
                                            if(stripY[j].size() == 4) { //If stripY[j] only contains 2 boxes
                                                //Do a straight swap, no need for AHC
                                                stripX[i].swap(stripY[j]);
                                                Swap(stripSumX[i], stripSumY[j]);
                                                feasible = 1;
                                            }
                                            else if(d == c + 2) { //If stripY[j] contains more than 2 boxes & the two chosen boxes in stripY[j] are adjacent
                                                //Only perform AHC on stripY[j]
                                                moveType = 1;
                                                PreAHCLS(i, a, b, j, c, d, feasible, swapType, moveType, allScores, itemWidths, itemNumbers,
                                                         stripSumX, stripX, stripSumY, stripY);
                                            }
                                            else { //IIf stripY[j] contains more than 2 boxes & boxes c and d are not adjacent
                                                moveType = 0;
                                                PreAHCLS(i, a, b, j, c, d, feasible, swapType, moveType, allScores, itemWidths, itemNumbers,
                                                         stripSumX, stripX, stripSumY, stripY);
                                            }
                                        }
                                        else if(stripY[j].size() == 4) { //If stripY[j] only contains 2 boxes but stripX[i] contains > 2 boxes
                                            if(b == a + 2) { //If the two chosen boxes in stripX[i] are adjacent to one another
                                                //Only perform AHC on stripX[i]
                                                moveType = 2;
                                                PreAHCLS(i, a, b, j, c, d, feasible, swapType, moveType, allScores, itemWidths, itemNumbers,
                                                         stripSumX, stripX, stripSumY, stripY);
                                            }
                                            else { //If boxes a and b are not adjacent
                                                moveType = 0;
                                                PreAHCLS(i, a, b, j, c, d, feasible, swapType, moveType, allScores, itemWidths, itemNumbers,
                                                         stripSumX, stripX, stripSumY, stripY);
                                            }
                                        }
                                        else { //If stripX[i].size() > 4 && stripY[j].size() > 4
                                            moveType = 0;
                                            PreAHCLS(i, a, b, j, c, d, feasible, swapType, moveType, allScores, itemWidths, itemNumbers,
                                                     stripSumX, stripX, stripSumY, stripY);
                                        }
                                        if(feasible == 1) {
                                            return;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

}

void PairSin(std::vector<int>& allScores, std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers,
             std::vector<int>& stripSumX, std::vector<std::vector<int> >& stripX, std::vector<int>& stripSumY, std::vector<std::vector<int> >& stripY) {

    int d = 0;
    int swapType, moveType, feasible = 0;

    /*SWAPPING A PAIR OF BOXES FROM SET STRIPX WITH ONE BOX FROM SET STRIPY*/
    for(int i = 0; i < stripX.size(); ++i) { //For each strip in the set stripX
        if(stripX[i].size() >= 4) { //If there are at least 2 boxes on stripX[i]
            for(int a = 0; a < stripX[i].size() - 3; a += 2) { //Starting from the first score on the first box until the first score on the penultimate box
                for(int b = a + 2; b < stripX[i].size() - 1; b += 2) { //Starting from the first score on the second box until the first score on the last box
                    int pairSizeX = itemWidths[stripX[i][a]][stripX[i][a + 1]] + itemWidths[stripX[i][b]][stripX[i][b + 1]]; //Sum box widths
                    //Check if there exists a box on a strip in set stripY whose width is larger than pairSizeX
                    for(int j = 0; j < stripY.size(); ++j) {
                        //Go through each box on stripY[j]
                        for(int c = 0; c < stripY[j].size() - 1; c += 2) { //Starting from the first score on the first box unil the first score on the last box
                            //Check if pairSizeX < width of box in stripY, and that box can fit onto strip
                            if(pairSizeX <= itemWidths[stripY[j][c]][stripY[j][c + 1]] && stripSumX[i] - pairSizeX + itemWidths[stripY[j][c]][stripY[j][c + 1]] <= stripWidth) {
                                swapType = 2;
                                if(stripX[i].size() == 4) { //If stripX[i] only contains 2 boxes
                                    if(stripY[j].size() == 2) { //If stripY[j] only contains 1 box
                                        //straight swap
                                        stripX[i].swap(stripY[j]);
                                        Swap(stripSumX[i], stripSumY[j]);
                                        feasible = 1;
                                    }
                                    else { //If stripY[j] contains more than 1 box
                                        //Only perform AHC on stripY[j]
                                        moveType = 1;
                                        PreAHCLS(i, a, b, j, c, d, feasible, swapType, moveType, allScores, itemWidths, itemNumbers,
                                                 stripSumX, stripX, stripSumY, stripY);
                                    }
                                }
                                else if(stripY[j].size() == 2) { //If stripY[j] only contains 1 box, but stripX[i] contains > 2 boxes
                                    if(b == a + 2) { //If the two boxes chosen from stripX[i] are adjacent to one another
                                        //Only perform AHC on stripX[i]
                                        moveType = 2;
                                        PreAHCLS(i, a, b, j, c, d, feasible, swapType, moveType, allScores, itemWidths, itemNumbers,
                                                 stripSumX, stripX, stripSumY, stripY);
                                    }
                                    else { //If the two boxes chosen from stripX[i] are not adjacent to one another
                                        moveType = 0;
                                        PreAHCLS(i, a, b, j, c, d, feasible, swapType, moveType, allScores, itemWidths, itemNumbers,
                                                 stripSumX, stripX, stripSumY, stripY);
                                    }
                                }
                                else { //If stripX[i].size() > 4 && stripY[j].size() > 2
                                    moveType = 0;
                                    PreAHCLS(i, a, b, j, c, d, feasible, swapType, moveType, allScores, itemWidths, itemNumbers,
                                             stripSumX, stripX, stripSumY, stripY);
                                }
                                if(feasible == 1) {
                                    return;
                                }
                            }
                        }
                    }
                }
            }
        }
    }


}

void SinSin(std::vector<int>& allScores, std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers,
            std::vector<int>& stripSumX, std::vector<std::vector<int> >& stripX, std::vector<int>& stripSumY, std::vector<std::vector<int> >& stripY) {

    int b = 0;
    int d = 0;
    int swapType, moveType, feasible = 0;

    /*SWAPPING ONE BOX FROM SET STRIPX WITH ONE BOX FROM SET STRIPY*/
    for(int i = 0; i < stripX.size(); ++i) { //For each strip in the set stripX
        for(int a = 0; a < stripX[i].size() - 1; a += 2) { // Starting from the first score on the first box until the first score on the last box
            for(int j = 0; j < stripY.size(); ++j) { // For each strip in the set stripY
                for(int c = 0; c < stripY[j].size() - 1; c += 2) { //Starting from the first score on the first box until the first score on the last box
                    //Check if boxwidth[a] < boxWidth[c] and that box can fit on strip
                    if(itemWidths[stripX[i][a]][stripX[i][a + 1]] < itemWidths[stripY[j][c]][stripY[j][c + 1]]
                       && stripSumX[i] - itemWidths[stripX[i][a]][stripX[i][a + 1]] + itemWidths[stripY[j][c]][stripY[j][c + 1]] <= stripWidth) {
                        swapType = 3;
                        if(stripX[i].size() == 2) { //If stripX[i] only contains 1 box
                            if(stripY[j].size() == 2) { //If stripY[j] only contains 1 box
                                //straight swap
                                stripX[i].swap(stripY[j]);
                                Swap(stripSumX[i], stripSumY[j]);
                                feasible = 1;
                            }
                            else { //If stripY[j] contains more than 1 box
                                //Only peform AHC on stripY[j]
                                moveType = 1;
                                PreAHCLS(i, a, b, j, c, d, feasible, swapType, moveType, allScores, itemWidths, itemNumbers,
                                         stripSumX, stripX, stripSumY, stripY);
                            }
                        }
                        else if(stripY[j].size() == 2) { //If stripY[j] only contains 1 box but stripX[i].size() > 2
                            //Only perform AHC on stripX[i]
                            moveType = 2;
                            PreAHCLS(i, a, b, j, c, d, feasible, swapType, moveType, allScores, itemWidths, itemNumbers,
                                     stripSumX, stripX, stripSumY, stripY);
                        }
                        else { //stripX[i].size() > 2 && stripY[j].size() > 2
                            moveType = 0;
                            PreAHCLS(i, a, b, j, c, d, feasible, swapType, moveType, allScores, itemWidths, itemNumbers,
                                     stripSumX, stripX, stripSumY, stripY);
                        }
                        if(feasible == 1) {
                            return;
                        }
                    }
                }
            }
        }
    }

}

void MoveSin(int& feasible, std::vector<int>& allScores, std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers,
             std::vector<int>& stripSumX, std::vector<std::vector<int> >& stripX, std::vector<int>& stripSumY, std::vector<std::vector<int> >& stripY) {

    int a = 0;
    int b = 0;
    int d = 0;
    int swapType, moveType;

    /*MOVING ONE BOX FROM SET STRIPY TO SET STRIPX*/
    for(int j = 0; j < stripY.size(); ++j) { //For each strip in the set stripY
        for(int c = 0; c < stripY[j].size() - 1; c += 2) { //Starting from the first score on the first box until the first score on the last box
            for(int i = 0; i < stripX.size(); ++i) { //For each strip in the set stripX
                if(stripSumX[i] + itemWidths[stripY[j][c]][stripY[j][c + 1]] <= stripWidth) {
                    swapType = 4;
                    if(stripY[j].size() == 2) { //If stripY[j] only contains one box
                        moveType = 1;
                        PreAHCLS(i, a, b, j, c, d, feasible, swapType, moveType, allScores, itemWidths, itemNumbers,
                                 stripSumX, stripX, stripSumY, stripY);
                    }
                    else if(stripY[j].size() == 4) {
                        moveType = 2;
                        PreAHCLS(i, a, b, j, c, d, feasible, swapType, moveType, allScores, itemWidths, itemNumbers,
                                 stripSumX, stripX, stripSumY, stripY);
                    }
                    else {
                        moveType = 0;
                        PreAHCLS(i, a, b, j, c, d, feasible, swapType, moveType, allScores, itemWidths, itemNumbers,
                                 stripSumX, stripX, stripSumY, stripY);
                    }
                    if(feasible == 1 || feasible == 2) {
                        return;
                    }
                }
            }
        }
    }
}

void PreAHCLS(int i1, int a1, int b1, int j1, int c1, int d1, int& feasible, int swapType, int moveType, std::vector<int>& allScores,
              std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<int>& stripSumX,
              std::vector<std::vector<int> >& stripX, std::vector<int>& stripSumY, std::vector<std::vector<int> >& stripY) {


    //region PairPair
    /**PAIRPAIR**/
    if(swapType == 1) {
        feasible = 0;
        if(moveType == 1) { //AHC on Y[j] only, X[i] = 2 items, Y[j] > 2 items, c and d adjacent
            ++numSubSCP;
            std::vector<int> finalY, scoresY, originalY, checkScoresY;
            for(int k = 0; k < stripY[j1].size(); ++k) {
                if(k == c1 || k == c1 + 1 || k == d1 || k == d1 + 1) {
                    continue;
                }
                scoresY.push_back(allScores[stripY[j1][k]]);
                originalY.push_back(stripY[j1][k]);
            }
            scoresY.push_back(allScores[stripX[i1][a1]]);
            originalY.push_back(stripX[i1][a1]);
            scoresY.push_back(allScores[stripX[i1][a1 + 1]]);
            originalY.push_back(stripX[i1][a1 + 1]);
            scoresY.push_back(allScores[stripX[i1][b1]]);
            originalY.push_back(stripX[i1][b1]);
            scoresY.push_back(allScores[stripX[i1][b1 + 1]]);
            originalY.push_back(stripX[i1][b1 + 1]);
            checkScoresY = scoresY;
            std::sort(checkScoresY.begin(), checkScoresY.end());
            if(checkScoresY[2] + checkScoresY[checkScoresY.size() - 1] < tau) {
                feasible = 0;
                ++numInfeas;
                return;
            }
            scoresY.push_back(tau);
            scoresY.push_back(tau);
            //Run AHC on scoresY, originalY, finalY
            AHC(feasible, scoresY, originalY, finalY);
            if(feasible == 1) { //i.e. solution has been found using ahc
                stripSumX[i1] = stripSumX[i1] - (itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1 + 1]])
                                + (itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]] + itemWidths[stripY[j1][d1]][stripY[j1][d1 + 1]]);
                stripSumY[j1] = stripSumY[j1] - (itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]] + itemWidths[stripY[j1][d1]][stripY[j1][d1 + 1]])
                                + (itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1 + 1]]);
                stripX[i1].clear(); //clear the two items/four score widths from stripX[i]
                stripX[i1].push_back(stripY[j1][c1]); //Put the two items/four score widths from stripY into stripX, they are adjacent so already
                stripX[i1].push_back(stripY[j1][c1 + 1]); //in a feasible alignment
                stripX[i1].push_back(stripY[j1][d1]);
                stripX[i1].push_back(stripY[j1][d1 + 1]);
                stripY[j1] = finalY;
                feasible = 1;
                ++numFeas;
                return;
            }
            else if(feasible == 0) { //ahc did not find solution
                feasible = 0;
                ++numInfeas;
                return; //keep going through for loops in PairPair function
            }
        } //end moveType == 1

        else if(moveType == 2) { //AHC on X[i] only, Y[j] = 2 items, X[i] > 2 items, a and b adjacent
            ++numSubSCP;
            std::vector<int> finalX, scoresX, originalX, checkScoresX;
            for(int k = 0; k < stripX[i1].size(); ++k) {
                if(k == a1 || k == a1 + 1 || k == b1 || k == b1 + 1) {
                    continue;
                }
                scoresX.push_back(allScores[stripX[i1][k]]);
                originalX.push_back(stripX[i1][k]);
            }
            scoresX.push_back(allScores[stripY[j1][c1]]);
            originalX.push_back(stripY[j1][c1]);
            scoresX.push_back(allScores[stripY[j1][c1 + 1]]);
            originalX.push_back(stripY[j1][c1 + 1]);
            scoresX.push_back(allScores[stripY[j1][d1]]);
            originalX.push_back(stripY[j1][d1]);
            scoresX.push_back(allScores[stripY[j1][d1 + 1]]);
            originalX.push_back(stripY[j1][d1 + 1]);
            checkScoresX = scoresX;
            std::sort(checkScoresX.begin(), checkScoresX.end());
            if(checkScoresX[2] + checkScoresX[checkScoresX.size() - 1] < tau) { //failed prelim checks, no need to run AHC
                feasible = 0;
                ++numInfeas;
                return;
            }
            scoresX.push_back(tau);
            scoresX.push_back(tau);
            //Run AHC on scoresX, originalX, finalX
            AHC(feasible, scoresX, originalX, finalX);
            if(feasible == 1) { //i.e. solution has been found using ahc
                stripSumX[i1] = stripSumX[i1] - (itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1 + 1]])
                                + (itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]] + itemWidths[stripY[j1][d1]][stripY[j1][d1 + 1]]);
                stripSumY[j1] = stripSumY[j1] - (itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]] + itemWidths[stripY[j1][d1]][stripY[j1][d1 + 1]])
                                + (itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1 + 1]]);
                stripY[j1].clear(); //clear the two items/four score widths from stripY[j]
                stripY[j1].push_back(stripX[i1][a1]); //Put the two items/four score widths from stripX into stripY, they are adjacent so already
                stripY[j1].push_back(stripX[i1][a1 + 1]); //in a feasible alignment
                stripY[j1].push_back(stripX[i1][b1]);
                stripY[j1].push_back(stripX[i1][b1 + 1]);
                stripX[i1] = finalX;
                feasible = 1;
                ++numFeas;
                return;
            }
            else if(feasible == 0) { //ahc did not find solution
                feasible = 0;
                ++numInfeas;
                return; //keep going through for loops in PairPair function
            }
        }

        else { /** moveType == 0, AHC on both X and Y. **/
            ++numSubSCP;
            std::vector<int> finalX, scoresX, originalX, checkScoresX;
            for(int k = 0; k < stripX[i1].size(); ++k) {
                if(k == a1 || k == a1 + 1 || k == b1 || k == b1 + 1) {
                    continue;
                }
                scoresX.push_back(allScores[stripX[i1][k]]);
                originalX.push_back(stripX[i1][k]);
            }
            scoresX.push_back(allScores[stripY[j1][c1]]);
            originalX.push_back(stripY[j1][c1]);
            scoresX.push_back(allScores[stripY[j1][c1 + 1]]);
            originalX.push_back(stripY[j1][c1 + 1]);
            scoresX.push_back(allScores[stripY[j1][d1]]);
            originalX.push_back(stripY[j1][d1]);
            scoresX.push_back(allScores[stripY[j1][d1 + 1]]);
            originalX.push_back(stripY[j1][d1 + 1]);
            checkScoresX = scoresX;
            std::sort(checkScoresX.begin(), checkScoresX.end());
            if(checkScoresX[2] + checkScoresX[checkScoresX.size() - 1] < tau) { //failed prelim checks, no need to run AHC
                feasible = 0;
                ++numInfeas;
                return;
            }
            scoresX.push_back(tau);
            scoresX.push_back(tau);
            //Run AHC on scoresX, originalX, finalX
            AHC(feasible, scoresX, originalX, finalX);
            if(feasible == 1) { //i.e. solution has been found using ahc for SX
                ++numFeas;
                ++numSubSCP;
                feasible = 0;
                std::vector<int> finalY, scoresY, originalY, checkScoresY;
                for(int k = 0; k < stripY[j1].size(); ++k) {
                    if(k == c1 || k == c1 + 1 || k == d1 || k == d1 + 1) {
                        continue;
                    }
                    scoresY.push_back(allScores[stripY[j1][k]]);
                    originalY.push_back(stripY[j1][k]);
                }
                scoresY.push_back(allScores[stripX[i1][a1]]);
                originalY.push_back(stripX[i1][a1]);
                scoresY.push_back(allScores[stripX[i1][a1 + 1]]);
                originalY.push_back(stripX[i1][a1 + 1]);
                scoresY.push_back(allScores[stripX[i1][b1]]);
                originalY.push_back(stripX[i1][b1]);
                scoresY.push_back(allScores[stripX[i1][b1 + 1]]);
                originalY.push_back(stripX[i1][b1 + 1]);
                checkScoresY = scoresY;
                std::sort(checkScoresY.begin(), checkScoresY.end());
                if(checkScoresY[2] + checkScoresY[checkScoresY.size() - 1] < tau) {
                    feasible = 0;
                    ++numInfeas;
                    return;
                }
                scoresY.push_back(tau);
                scoresY.push_back(tau);
                //Run AHC on scoresY, originalY, finalY
                AHC(feasible, scoresY, originalY, finalY);
                if(feasible == 1) { //i.e. solution has been found using ahc
                    stripSumX[i1] = stripSumX[i1] - (itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1 + 1]])
                                    + (itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]] + itemWidths[stripY[j1][d1]][stripY[j1][d1 + 1]]);
                    stripSumY[j1] = stripSumY[j1] - (itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]] + itemWidths[stripY[j1][d1]][stripY[j1][d1 + 1]])
                                    + (itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1 + 1]]);
                    stripX[i1] = finalX;
                    stripY[j1] = finalY;
                    ++numFeas;
                    feasible = 1;
                    return;
                }
                else if(feasible == 0) { //ahc did not find solution for SY
                    feasible = 0;
                    ++numInfeas;
                    return;
                }
            }
            else if(feasible == 0) { //ahc did not find solution for SX
                feasible = 0;
                ++numInfeas;
                return;
            }
        }
    } //End swapType = 1
    //endregion

    //region PairSin
    /**PAIRSIN**/
    if(swapType == 2) {
        feasible = 0;
        if(moveType == 1) { //AHC on Y[j] only.
            ++numSubSCP;
            std::vector<int> finalY, scoresY, originalY, checkScoresY;
            for(int k = 0; k < stripY[j1].size(); ++k) {
                if(k == c1 || k == c1 + 1) {
                    continue;
                }
                scoresY.push_back(allScores[stripY[j1][k]]);
                originalY.push_back(stripY[j1][k]);
            }
            scoresY.push_back(allScores[stripX[i1][a1]]);
            originalY.push_back(stripX[i1][a1]);
            scoresY.push_back(allScores[stripX[i1][a1 + 1]]);
            originalY.push_back(stripX[i1][a1 + 1]);
            scoresY.push_back(allScores[stripX[i1][b1]]);
            originalY.push_back(stripX[i1][b1]);
            scoresY.push_back(allScores[stripX[i1][b1 + 1]]);
            originalY.push_back(stripX[i1][b1 + 1]);
            checkScoresY = scoresY;
            std::sort(checkScoresY.begin(), checkScoresY.end());
            if(checkScoresY[2] + checkScoresY[checkScoresY.size() - 1] < tau) {
                feasible = 0;
                ++numInfeas;
                return;
            }
            scoresY.push_back(tau);
            scoresY.push_back(tau);
            //Run AHC on scoresY, originalY, finalY
            AHC(feasible, scoresY, originalY, finalY);
            if(feasible == 1) { //i.e. solution has been found using ahc
                stripSumX[i1] = stripSumX[i1] - (itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1 + 1]])
                                + itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]];
                stripSumY[j1] = stripSumY[j1] - itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]] + (itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]]
                                                                                                  + itemWidths[stripX[i1][b1]][stripX[i1][b1 + 1]]);
                stripX[i1].clear(); //clear the two items/four score widths from stripX[i]
                stripX[i1].push_back(stripY[j1][c1]); //Put the two items/four score widths from stripY into stripX, they are adjacent so already
                stripX[i1].push_back(stripY[j1][c1 + 1]); //in a feasible alignment
                stripY[j1] = finalY;
                feasible = 1;
                ++numFeas;
                return;
            }
            else if(feasible == 0) { //ahc did not find solution
                feasible = 0;
                ++numInfeas;
                return; //keep going through for loops in PairPair function
            }
        } //end moveType == 1

        else if(moveType == 2) { //AHC on X[i] only.
            ++numSubSCP;
            std::vector<int> finalX, scoresX, originalX, checkScoresX;
            for(int k = 0; k < stripX[i1].size(); ++k) {
                if(k == a1 || k == a1 + 1 || k == b1 || k == b1 + 1) {
                    continue;
                }
                scoresX.push_back(allScores[stripX[i1][k]]);
                originalX.push_back(stripX[i1][k]);
            }
            scoresX.push_back(allScores[stripY[j1][c1]]);
            originalX.push_back(stripY[j1][c1]);
            scoresX.push_back(allScores[stripY[j1][c1 + 1]]);
            originalX.push_back(stripY[j1][c1 + 1]);
            checkScoresX = scoresX;
            std::sort(checkScoresX.begin(), checkScoresX.end());
            if(checkScoresX[2] + checkScoresX[checkScoresX.size() - 1] < tau) {
                feasible = 0;
                ++numInfeas;
                return;
            }
            scoresX.push_back(tau);
            scoresX.push_back(tau);
            //Run AHC on scoresX, originalX, finalX
            AHC(feasible, scoresX, originalX, finalX);
            if(feasible == 1) { //i.e. solution has been found using ahc
                stripSumX[i1] = stripSumX[i1] - (itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1 + 1]])
                                + itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]];
                stripSumY[j1] = stripSumY[j1] - itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]] + (itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]]
                                                                                                  + itemWidths[stripX[i1][b1]][stripX[i1][b1 + 1]]);
                stripY[j1].clear(); //clear the two items/four score widths from stripY[j]
                stripY[j1].push_back(stripX[i1][a1]); //Put the two items/four score widths from stripX into stripY, they are adjacent so already
                stripY[j1].push_back(stripX[i1][a1 + 1]); //in a feasible alignment
                stripY[j1].push_back(stripX[i1][b1]);
                stripY[j1].push_back(stripX[i1][b1 + 1]);
                stripX[i1] = finalX;
                feasible = 1;
                ++numFeas;
                return;
            }
            else if(feasible == 0) { //ahc did not find solution
                feasible = 0;
                ++numInfeas;
                return; //keep going through for loops in PairPair function
            }
        }

        else { /** moveType == 0, AHC on both X and Y. **/
            ++numSubSCP;
            std::vector<int> finalX, scoresX, originalX, checkScoresX;
            for(int k = 0; k < stripX[i1].size(); ++k) {
                if(k == a1 || k == a1 + 1 || k == b1 || k == b1 + 1) {
                    continue;
                }
                scoresX.push_back(allScores[stripX[i1][k]]);
                originalX.push_back(stripX[i1][k]);
            }
            scoresX.push_back(allScores[stripY[j1][c1]]);
            originalX.push_back(stripY[j1][c1]);
            scoresX.push_back(allScores[stripY[j1][c1 + 1]]);
            originalX.push_back(stripY[j1][c1 + 1]);
            checkScoresX = scoresX;
            std::sort(checkScoresX.begin(), checkScoresX.end());
            if(checkScoresX[2] + checkScoresX[checkScoresX.size() - 1] < tau) {
                feasible = 0;
                ++numInfeas;
                return;
            }
            scoresX.push_back(tau);
            scoresX.push_back(tau);
            //Run AHC on scoresX, originalX, finalX
            AHC(feasible, scoresX, originalX, finalX);
            if(feasible == 1) { //i.e. solution has been found using ahc
                ++numFeas;
                ++numSubSCP;
                feasible = 0;
                std::vector<int> finalY, scoresY, originalY, checkScoresY;
                for(int k = 0; k < stripY[j1].size(); ++k) {
                    if(k == c1 || k == c1 + 1) {
                        continue;
                    }
                    scoresY.push_back(allScores[stripY[j1][k]]);
                    originalY.push_back(stripY[j1][k]);
                }
                scoresY.push_back(allScores[stripX[i1][a1]]);
                originalY.push_back(stripX[i1][a1]);
                scoresY.push_back(allScores[stripX[i1][a1 + 1]]);
                originalY.push_back(stripX[i1][a1 + 1]);
                scoresY.push_back(allScores[stripX[i1][b1]]);
                originalY.push_back(stripX[i1][b1]);
                scoresY.push_back(allScores[stripX[i1][b1 + 1]]);
                originalY.push_back(stripX[i1][b1 + 1]);
                checkScoresY = scoresY;
                std::sort(checkScoresY.begin(), checkScoresY.end());
                if(checkScoresY[2] + checkScoresY[checkScoresY.size() - 1] < tau) {
                    feasible = 0;
                    ++numInfeas;
                    return;
                }
                scoresY.push_back(tau);
                scoresY.push_back(tau);
                //Run AHC on scoresY, originalY, finalY
                AHC(feasible, scoresY, originalY, finalY);
                if(feasible == 1) { //i.e. solution has been found using ahc
                    stripSumX[i1] = stripSumX[i1] - (itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1 + 1]])
                                    + itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]];
                    stripSumY[j1] = stripSumY[j1] - itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]] + (itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]]
                                    + itemWidths[stripX[i1][b1]][stripX[i1][b1 + 1]]);
                    stripX[i1] = finalX;
                    stripY[j1] = finalY;
                    feasible = 1;
                    ++numFeas;
                    return;
                }
                else if(feasible == 0) { //ahc did not find solution for SY
                    feasible = 0;
                    ++numInfeas;
                    return;
                }
            }
            else if(feasible == 0) { //ahc did not find solution for SX
                feasible = 0;
                ++numInfeas;
                return;
            }
        }
    } //End swapType = 2
    //endregion

    //region SinSin
    /**SINSIN**/
    if(swapType == 3) {
        feasible = 0;
        if(moveType == 1) { //AHC on Y[j] only.
            ++numSubSCP;
            std::vector<int> finalY, scoresY, originalY, checkScoresY;
            for(int k = 0; k < stripY[j1].size(); ++k) {
                if(k == c1 || k == c1 + 1) {
                    continue;
                }
                scoresY.push_back(allScores[stripY[j1][k]]);
                originalY.push_back(stripY[j1][k]);
            }
            scoresY.push_back(allScores[stripX[i1][a1]]);
            originalY.push_back(stripX[i1][a1]);
            scoresY.push_back(allScores[stripX[i1][a1 + 1]]);
            originalY.push_back(stripX[i1][a1 + 1]);
            checkScoresY = scoresY;
            std::sort(checkScoresY.begin(), checkScoresY.end());
            if(checkScoresY[2] + checkScoresY[checkScoresY.size() - 1] < tau) {
                feasible = 0;
                ++numInfeas;
                return;
            }
            scoresY.push_back(tau);
            scoresY.push_back(tau);
            //Run AHC on scoresY, originalY, final
            AHC(feasible, scoresY, originalY, finalY);
            if(feasible == 1) { //i.e. solution has been found using ahc
                stripSumX[i1] = stripSumX[i1] - itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]]
                                + itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]];
                stripSumY[j1] = stripSumY[j1] - itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]]
                                + itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]];
                stripX[i1].clear(); //clear the two items/four score widths from stripX[i]
                stripX[i1].push_back(stripY[j1][c1]); //Put the two items/four score widths from stripY into stripX, they are adjacent so already
                stripX[i1].push_back(stripY[j1][c1 + 1]); //in a feasible alignment
                stripY[j1] = finalY;
                feasible = 1;
                ++numFeas;
                return;
            }
            else if(feasible == 0) { //ahc did not find solution
                feasible = 0;
                ++numInfeas;
                return; //keep going through for loops in PairPair function
            }
        } //end moveType == 1

        else if(moveType == 2) { //AHC on X[i] only.
            ++numSubSCP;
            std::vector<int> finalX, scoresX, originalX, checkScoresX;
            for(int k = 0; k < stripX[i1].size(); ++k) {
                if(k == a1 || k == a1 + 1) {
                    continue;
                }
                scoresX.push_back(allScores[stripX[i1][k]]);
                originalX.push_back(stripX[i1][k]);
            }
            scoresX.push_back(allScores[stripY[j1][c1]]);
            originalX.push_back(stripY[j1][c1]);
            scoresX.push_back(allScores[stripY[j1][c1 + 1]]);
            originalX.push_back(stripY[j1][c1 + 1]);
            checkScoresX = scoresX;
            std::sort(checkScoresX.begin(), checkScoresX.end());
            if(checkScoresX[2] + checkScoresX[checkScoresX.size() - 1] < tau) {
                feasible = 0;
                ++numInfeas;
                return;
            }
            scoresX.push_back(tau);
            scoresX.push_back(tau);
            //Run AHC on scoresX, originalX, finalX
            AHC(feasible, scoresX, originalX, finalX);
            if(feasible == 1) { //i.e. solution has been found using ahc
                stripSumX[i1] = stripSumX[i1] - itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]]
                                + itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]];
                stripSumY[j1] = stripSumY[j1] - itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]]
                                + itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]];
                stripY[j1].clear(); //clear the two items/four score widths from stripY[j]
                stripY[j1].push_back(stripX[i1][a1]); //Put the two items/four score widths from stripX into stripY, they are adjacent so already
                stripY[j1].push_back(stripX[i1][a1 + 1]); //in a feasible alignment
                stripX[i1] = finalX;
                feasible = 1;
                ++numFeas;
                return;
            }
            else if(feasible == 0) { //ahc did not find solution
                feasible = 0;
                ++numInfeas;
                return; //keep going through for loops in PairPair function
            }
        }

        else { /** moveType == 0, AHC on both X and Y. **/
            ++numSubSCP;
            std::vector<int> finalX, scoresX, originalX, checkScoresX;
            for(int k = 0; k < stripX[i1].size(); ++k) {
                if(k == a1 || k == a1 + 1) {
                    continue;
                }
                scoresX.push_back(allScores[stripX[i1][k]]);
                originalX.push_back(stripX[i1][k]);
            }
            scoresX.push_back(allScores[stripY[j1][c1]]);
            originalX.push_back(stripY[j1][c1]);
            scoresX.push_back(allScores[stripY[j1][c1 + 1]]);
            originalX.push_back(stripY[j1][c1 + 1]);
            checkScoresX = scoresX;
            std::sort(checkScoresX.begin(), checkScoresX.end());
            if(checkScoresX[2] + checkScoresX[checkScoresX.size() - 1] < tau) {
                feasible = 0;
                ++numInfeas;
                return;
            }
            scoresX.push_back(tau);
            scoresX.push_back(tau);
            //Run AHC on scoresX, originalX, finalX
            AHC(feasible, scoresX, originalX, finalX);
            if(feasible == 1) { //i.e. solution has been found using ahc
                ++numFeas;
                ++numSubSCP;
                feasible = 0;
                std::vector<int> finalY, scoresY, originalY, checkScoresY;
                for(int k = 0; k < stripY[j1].size(); ++k) {
                    if(k == c1 || k == c1 + 1) {
                        continue;
                    }
                    scoresY.push_back(allScores[stripY[j1][k]]);
                    originalY.push_back(stripY[j1][k]);
                }
                scoresY.push_back(allScores[stripX[i1][a1]]);
                originalY.push_back(stripX[i1][a1]);
                scoresY.push_back(allScores[stripX[i1][a1 + 1]]);
                originalY.push_back(stripX[i1][a1 + 1]);
                checkScoresY = scoresY;
                std::sort(checkScoresY.begin(), checkScoresY.end());
                if(checkScoresY[2] + checkScoresY[checkScoresY.size() - 1] < tau) {
                    feasible = 0;
                    ++numInfeas;
                    return;
                }
                scoresY.push_back(tau);
                scoresY.push_back(tau);
                //Run AHC on scoresY, originalY, final
                AHC(feasible, scoresY, originalY, finalY);
                if(feasible == 1) { //i.e. solution has been found using ahc
                    stripSumX[i1] = stripSumX[i1] - itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]]
                                    + itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]];
                    stripSumY[j1] = stripSumY[j1] - itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]]
                                    + itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]];
                    stripX[i1] = finalX;
                    stripY[j1] = finalY;
                    feasible = 1;
                    ++numFeas;
                    return;
                }
                else if(feasible == 0) { //ahc did not find solution for SY
                    feasible = 0;
                    ++numInfeas;
                    return;
                }
            }
            else if(feasible == 0) { //ahc did not find solution for SX
                feasible = 0;
                ++numInfeas;
                return;
            }
        }
    } //End swapType = 3
    //endregion

    //region MoveSin
    /**MOVESIN**/
    if(swapType == 4) {
        feasible = 0;
        if(moveType == 1 || moveType == 2) { //AHC on X[i] only.
            ++numSubSCP;
            std::vector<int> finalX, scoresX, originalX, checkScoresX;
            for(int k = 0; k < stripX[i1].size(); ++k) {
                scoresX.push_back(allScores[stripX[i1][k]]);
                originalX.push_back(stripX[i1][k]);
            }
            scoresX.push_back(allScores[stripY[j1][c1]]);
            originalX.push_back(stripY[j1][c1]);
            scoresX.push_back(allScores[stripY[j1][c1 + 1]]);
            originalX.push_back(stripY[j1][c1 + 1]);
            checkScoresX = scoresX;
            std::sort(checkScoresX.begin(), checkScoresX.end());
            if(checkScoresX[2] + checkScoresX[checkScoresX.size() - 1] < tau) {
                feasible = 0;
                ++numInfeas;
                return;
            }
            scoresX.push_back(tau);
            scoresX.push_back(tau);
            //Run AHC on scoresX, originalX, finalX
            AHC(feasible, scoresX, originalX, finalX);
            if(feasible == 1) { //i.e. solution has been found using ahc
                stripSumX[i1] += itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]];
                stripSumY[j1] -= itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]];
                stripX[i1] = finalX; //make stripX be the strip of scores in the feasPacking struct
                ++numFeas;
                if(moveType == 1) {
                    stripY.erase(stripY.begin() + j1);
                    stripSumY.erase(stripSumY.begin() + j1);
                    if(stripY.empty()) {
                        feasible = 2;
                        return;
                    }
                    else {
                        feasible = 1;
                        return;
                    }
                }
                else if(moveType == 2) {
                    if(c1 == 0) {
                        stripY[j1].erase(stripY[j1].begin(), stripY[j1].begin() + 2);
                        feasible = 1;
                        return;
                    }
                    else if(c1 == 2) {
                        stripY[j1].pop_back();
                        stripY[j1].pop_back();
                        feasible = 1;
                        return;
                    }
                    else {
                        std::cerr << "[ERROR]: c1 in stripY[j1] is neither 0 nor 2, check that stripY[j1].size() == 4\n";
                        exit(1);
                    }
                }
            }
            else if(feasible == 0) { //ahc did not find solution
                feasible = 0;
                ++numInfeas;
                return; //keep going through for loops in PairPair function
            }
        }

        else { /** moveType == 0, AHC on both X and Y. **/
            ++numSubSCP;
            std::vector<int> finalX, scoresX, originalX, checkScoresX;
            for(int k = 0; k < stripX[i1].size(); ++k) {
                scoresX.push_back(allScores[stripX[i1][k]]);
                originalX.push_back(stripX[i1][k]);
            }
            scoresX.push_back(allScores[stripY[j1][c1]]);
            originalX.push_back(stripY[j1][c1]);
            scoresX.push_back(allScores[stripY[j1][c1 + 1]]);
            originalX.push_back(stripY[j1][c1 + 1]);
            checkScoresX = scoresX;
            std::sort(checkScoresX.begin(), checkScoresX.end());
            if(checkScoresX[2] + checkScoresX[checkScoresX.size() - 1] < tau) {
                feasible = 0;
                ++numInfeas;
                return;
            }
            scoresX.push_back(tau);
            scoresX.push_back(tau);
            //Run AHC on scoresX, originalX, finalX
            AHC(feasible, scoresX, originalX, finalX);
            if(feasible == 1) { //i.e. solution has been found using ahc
                ++numFeas;
                ++numSubSCP;
                feasible = 0;
                std::vector<int> finalY, scoresY, originalY, checkScoresY;
                for(int k = 0; k < stripY[j1].size(); ++k) {
                    if(k == c1 || k == c1 + 1) {
                        continue;
                    }
                    scoresY.push_back(allScores[stripY[j1][k]]);
                    originalY.push_back(stripY[j1][k]);
                }
                checkScoresY = scoresY;
                std::sort(checkScoresY.begin(), checkScoresY.end());
                if(checkScoresY[2] + checkScoresY[checkScoresY.size() - 1] < tau) {
                    feasible = 0;
                    ++numInfeas;
                    return;
                }
                scoresY.push_back(tau);
                scoresY.push_back(tau);
                //Run AHC on scoresY, originalY, finalY
                AHC(feasible, scoresY, originalY, finalY);
                if(feasible == 1) { //i.e. solution has been found using ahc
                    stripSumX[i1] += itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]];
                    stripSumY[j1] -= itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]];
                    stripX[i1] = finalX;
                    stripY[j1] = finalY;
                    feasible = 1;
                    ++numFeas;
                    return;
                }
                else if(feasible == 0) { //ahc did not find solution for SY
                    feasible = 0;
                    ++numInfeas;
                    return;
                }
            }
            else if(feasible == 0) { //ahc did not find solution for SX
                feasible = 0;
                ++numInfeas;
                return;
            }
        }
    } //End swapType = 4
    //endregion

}

void AHC(int& feasible, std::vector<int>& scores, std::vector<int>& original, std::vector<int>& final) {

    feasible = 0;
    int nScores = scores.size();
    int nItem = scores.size() / 2;
    std::vector<int> order;
    std::vector<int> altHam;
    std::vector<int> B(nScores, INT_MAX);
    std::vector<std::vector<int> > adjMat(nScores, std::vector<int>(nScores, 0));

    InitInstance(nScores, adjMat, scores, order, B);

    int matchSize = 0;
    std::vector<int> matchList(nScores, INT_MAX);
    std::vector<int> cycleVertex(nScores, 1);

    MCM(nScores, matchSize, adjMat, B, matchList, cycleVertex);
    if(matchSize < nItem) {
        feasible = 0;
        return;
    }

    int nCycles = 0;
    std::vector<std::vector<int> > mpStructure;
    MPS(nScores, nCycles, B, matchList, mpStructure);
    if(mpStructure[0].size() == nScores) {
        for(int j = 2; j < mpStructure[0].size(); ++j) {
            altHam.push_back(mpStructure[0][j]);
        }

        for(int i = 0; i < altHam.size(); ++i) {
            final.push_back(original[order[altHam[i]]]);
        }
        feasible = 1;
        return;
    }

    BCR(nScores, feasible, matchSize, nCycles, B, matchList, cycleVertex, mpStructure, adjMat, altHam);
    if(feasible == 0){
        return;
    }
    else if(feasible == 1){
        for(int i = 0; i < altHam.size(); ++i) {
            final.push_back(original[order[altHam[i]]]);
        }
        return;
    }


}

void InitInstance(int nScores, std::vector<std::vector<int> >& adjMat, std::vector<int>& scores, std::vector<int>& order, std::vector<int>& B) {

    std::vector<int> invOrder(nScores);

    for(int i = 0; i < nScores; ++i) {
        order.push_back(i);
    }

    for(int i = 1; i < nScores; ++i) {
        for(int j = i - 1; j >= 0; --j) {
            if(scores[i] < scores[order[j]]) {
                order[j + 1] = order[j];
                order[j] = i;
            }
        }
    }

    for(int i = 0; i < nScores; ++i) {
        invOrder[order[i]] = i;
    }

    for(int i = 0; i < nScores - 1; i += 2) {
        adjMat[invOrder[i]][invOrder[i + 1]] = 2;
        adjMat[invOrder[i + 1]][invOrder[i]] = 2;
    }

    std::sort(scores.begin(), scores.end());

    for(int i = 0; i < scores.size() - 1; ++i) {
        for(int j = i + 1; j < scores.size(); ++j) {
            if(scores[i] + scores[j] >= tau && adjMat[i][j] != 2) {
                adjMat[i][j] = 1;
                adjMat[j][i] = 1;
            }
        }
    }

    for(int i = 0; i < nScores; ++i) {
        for(int j = 0; j < nScores; ++j) {
            if(adjMat[i][j] == 2) {
                B[i] = j;
            }
        }
    }

}

void MCM(int nScores, int& matchSize, std::vector<std::vector<int> >& adjMat, std::vector<int>& B, std::vector<int>& matchList,
         std::vector<int>& cycleVertex) {

    int vacantFlag = 0;
    int lastMatch = INT_MAX;

    for(int i = 0; i < nScores; ++i) {
        vacantFlag = 0;
        if(matchList[i] == INT_MAX) {
            for(int j = nScores - 1; j > i; --j) {
                if(adjMat[i][j] == 1 && matchList[j] == INT_MAX) {
                    matchList[i] = j;
                    matchList[j] = i;
                    lastMatch = i;
                    ++matchSize;
                    if(vacantFlag == 1) {
                        cycleVertex[i] = INT_MAX;
                        cycleVertex[j] = INT_MAX;
                    }
                    break;
                }
                else if(adjMat[i][j] == 2 && matchList[j] == INT_MAX) {
                    vacantFlag = 1;
                }
            }
            if(matchList[i] == INT_MAX) {
                if((matchList[B[i]] == INT_MAX) && (lastMatch != INT_MAX)
                   && (B[i] > i) && (adjMat[lastMatch][B[i]] == 1)) {
                    matchList[i] = matchList[lastMatch];
                    matchList[lastMatch] = B[i];
                    matchList[B[i]] = lastMatch;
                    matchList[matchList[i]] = i;
                    cycleVertex[lastMatch] = INT_MAX;
                    cycleVertex[B[i]] = INT_MAX;
                    lastMatch = i;
                    ++matchSize;
                }
            }
        }
    }

}

void MPS(int nScores, int& nCycles, std::vector<int>& B, std::vector<int>& matchList, std::vector<std::vector<int> >& mpStructure) {

    int current;
    int smallest = nScores - 2;
    std::vector<int> temp;
    std::vector<int> checked(nScores, 0);

    do {
        current = smallest;
        temp.clear();
        do {
            temp.push_back(current);
            checked[current] = 1;
            temp.push_back(B[current]);
            checked[B[current]] = 1;
            current = matchList[B[current]];
        } while(current != smallest);

        mpStructure.push_back(temp);
        temp.clear();

        for(int i = 0; i < nScores; ++i) {
            if(checked[i] == 0) {
                smallest = i;
                break;
            }
        }
    } while(smallest != current);

    nCycles = mpStructure.size();

}

void BCR(int nScores, int& feasible, int matchSize, int nCycles, std::vector<int>& B, std::vector<int>& matchList,
         std::vector<int>& cycleVertex, std::vector<std::vector<int> >& mpStructure, std::vector<std::vector<int> >& adjMat, std::vector<int>& altHam) {

    feasible = 0;

    for(int i = 0; i < mpStructure.size(); ++i) {
        for(int j = 0; j < mpStructure[i].size(); ++j) {
            if(cycleVertex[mpStructure[i][j]] != INT_MAX) {
                cycleVertex[mpStructure[i][j]] = i;
            }
        }
    }

    std::vector<int> edge;
    for(int i = 0; i < matchSize; ++i) {
        while(cycleVertex[i] == INT_MAX) {
            ++i;
        }
        edge.push_back(i);
    }
    int nEdges = edge.size();

    //To find one set only
    int k = 0;
    std::vector<int> temp;
    std::vector<int> S(nCycles, 0);
    std::vector<std::vector<int> > R;

    while(k < nEdges - 2 && (adjMat[edge[k]][matchList[edge[k + 1]]] != 1 || cycleVertex[edge[k]] == cycleVertex[edge[k + 1]])) {
        ++k;
    }
    if(adjMat[edge[k]][matchList[edge[k + 1]]] == 1 && cycleVertex[edge[k]] != cycleVertex[edge[k + 1]]) {
        temp.push_back(edge[k]);
        S[cycleVertex[edge[k]]] = 1;
        while(k < nEdges - 1 && adjMat[edge[k]][matchList[edge[k + 1]]] == 1 && S[cycleVertex[edge[k + 1]]] == 0) {
            ++k;
            temp.push_back(edge[k]);
            S[cycleVertex[edge[k]]] = 1;
        }
        R.push_back(temp);
        temp.clear();
    }

    //If no set found
    if(R.size() == 0){
        feasible = 0;
        return;
    }
    else if(R[0].size() == nCycles){ //If the set found covers all cycles, |R1| = l
        CP(nScores, 0, false, matchList, B, altHam, R);
        feasible = 1;
        return;
    }

    //Otherwise, use this single set R1 to find overlapping sets:
    int SSum = 0;
    std::vector<int> S2;
    S2 = S;

    for(int i = 0; i < S.size(); ++i){
        SSum = SSum + S[i];
    }
    int SSum2 = SSum;

    //Remove edges in R1 from edge list
    for(int i = 0; i < R[0].size(); ++i){
        for(int j = 0; j < edge.size(); ++j){
            if(edge[j] == R[0][i]){
                edge.erase(edge.begin() + j);
                break;
            }
        }
    }
    nEdges = edge.size();

    //Now find overlapping sets:
    int lastRow = R.size();
    bool added = false;
    std::vector<int> temp2;
    temp.clear();


    while(!added){
        int k = 0;
        do {
            while(k < nEdges - 2 && (adjMat[edge[k]][matchList[edge[k+1]]] != 1 || cycleVertex[edge[k]] == cycleVertex[edge[k+1]])){
                ++k;
            }
            if(adjMat[edge[k]][matchList[edge[k + 1]]] == 1 && cycleVertex[edge[k]] != cycleVertex[edge[k + 1]]
               && ((S[cycleVertex[edge[k]]] == 0 && S[cycleVertex[edge[k + 1]]] == 1)
                   || (S[cycleVertex[edge[k]]] == 1 && S[cycleVertex[edge[k + 1]]] == 0))) {
                temp.push_back(edge[k]);
                temp.push_back(edge[k + 1]);
                S[cycleVertex[edge[k]]] = 1;
                S[cycleVertex[edge[k + 1]]] = 1;
                ++SSum;
                if(SSum < nCycles) {
                    ++k;
                    while(k < nEdges - 1 && S[cycleVertex[edge[k + 1]]] == 0 && adjMat[edge[k]][matchList[edge[k + 1]]] == 1) {
                        ++k;
                        temp.push_back(edge[k]);
                        S[cycleVertex[edge[k]]] = 1;
                        ++SSum;
                    }
                }
                R.push_back(temp);
                added = true;
                temp.clear();
                S2 = S;
                SSum2 = SSum;
            }

            else if(adjMat[edge[k]][matchList[edge[k+1]]] == 1 && cycleVertex[edge[k]] != cycleVertex[edge[k+1]]
                    && (S[cycleVertex[edge[k]]] == 0 && S[cycleVertex[edge[k+1]]] == 0)){
                bool overlap = false;
                temp2.push_back(edge[k]);
                temp2.push_back(edge[k+1]);
                S2[cycleVertex[edge[k]]] = 1;
                S2[cycleVertex[edge[k+1]]] = 1;
                SSum2 += 2;
                ++k;
                while(k < nEdges - 1 && S2[cycleVertex[edge[k+1]]] == 0 && adjMat[edge[k]][matchList[edge[k+1]]] == 1){
                    ++k;
                    temp2.push_back(edge[k]);
                    S2[cycleVertex[edge[k]]] = 1;
                    ++SSum2;
                }
                if(k < nEdges - 1 && S2[cycleVertex[edge[k+1]]] == 1 && adjMat[edge[k]][matchList[edge[k+1]]] == 1){
                    if(S[cycleVertex[edge[k+1]]] == 1){
                        overlap = true;
                        ++k;
                        temp2.push_back(edge[k]);
                        S2[cycleVertex[edge[k]]] = 1;
                        //++SSum2;
                        if(SSum2 < nCycles){
                            ++k;
                            while(k < nEdges - 1 && S2[cycleVertex[edge[k+1]]] == 0 && adjMat[edge[k]][matchList[edge[k+1]]] == 1){
                                ++k;
                                temp2.push_back(edge[k]);
                                S2[cycleVertex[edge[k]]] = 1;
                                ++SSum2;
                            } //end while loop
                        } //end if SSum2 < nCycles
                    }
                    else if (S[cycleVertex[edge[k+1]]] == 0){
                        //no overlap, edge from cycle that already has edge in current set, exit.
                        overlap == false;
                    }

                } //end if SSet2 == 1 overlap
                if(overlap){
                    R.push_back(temp2);
                    temp2.clear();
                    temp.clear();
                    S = S2;
                    SSum = SSum2;
                    added = true;
                }
                else{
                    temp2.clear();
                    S2 = S;
                    SSum2 = SSum;
                    temp.clear();
                }
            }//end else if two edges found == 0
            ++k;
        } while(k < nEdges - 1 && SSum < nCycles); //End do while loop


        if(SSum == nCycles){
            CP(nScores, 0, true, matchList, B, altHam, R);
            feasible = 1;
            return;
        }
        else if(added) { // &&SSum < nCycles
            for(int i = lastRow; i < R.size(); ++i) {
                for(int j = 0; j < R[i].size(); ++j) {
                    for(int l = 0; l < edge.size(); ++l) {
                        if(edge[l] == R[i][j]) {
                            edge.erase(edge.begin() + l);
                            break;
                        }
                    }
                }
            }
            lastRow = R.size();
            added = false;
            nEdges = edge.size();
            if(nEdges < 2){ //not enough edges to look for another set, and the sets currently found do not cover all cycles, infeasible.
                feasible = 0;
                return;
            }
        }
        else if(!added){
            feasible = 0;
            return;
        }
    }// End while !added


}//end BCR

void CP(int nScores, int full, bool multiple, std::vector<int>& matchList, std::vector<int>& B, std::vector<int>& altHam, std::vector<std::vector<int> >& RStar){

    std::vector<int> connectML;

    if(!multiple){ //type 0, only one set needed
        connectML = matchList;
        for(int v = 0; v < RStar[full].size() - 1; ++v) {
            connectML[RStar[full][v]] = matchList[RStar[full][v + 1]];
            connectML[matchList[RStar[full][v + 1]]] = RStar[full][v];
        }
        connectML[RStar[full][RStar[full].size() - 1]] = matchList[RStar[full][0]];
        connectML[matchList[RStar[full][0]]] = RStar[full][RStar[full].size() - 1];

        int current = nScores - 2;
        do {
            altHam.push_back(current);
            altHam.push_back(B[current]);
            current = connectML[B[current]];
        } while(altHam.size() < nScores);

        altHam.erase(altHam.begin(), altHam.begin() + 2);
        return;
    }
    else if(multiple){
        connectML = matchList;
        for(int u = 0; u < RStar.size(); ++u) {
            for(int v = 0; v < RStar[u].size() - 1; ++v) {
                connectML[RStar[u][v]] = matchList[RStar[u][v + 1]];
                connectML[matchList[RStar[u][v + 1]]] = RStar[u][v];
            }
            connectML[RStar[u][RStar[u].size() - 1]] = matchList[RStar[u][0]];
            connectML[matchList[RStar[u][0]]] = RStar[u][RStar[u].size() - 1];
        }

        int current = nScores - 2;
        do {
            altHam.push_back(current);
            altHam.push_back(B[current]);
            current = connectML[B[current]];
        } while(altHam.size() < nScores);

        altHam.erase(altHam.begin(), altHam.begin() + 2);
        return;
    }
}

void EA(const Timer& timer, int numScores, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
        std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<std::vector<int> >& populationSum,
        std::vector<std::vector<std::vector<int> > >& population) {

    std::vector<std::vector<int> > stripX;
    std::vector<std::vector<int> > stripY;
    std::vector<std::vector<int> > offspring;
    std::vector<int> stripSumX;
    std::vector<int> stripSumY;
    std::vector<int> offspringSum;

    //Choose two solutions from population at random
    int k = rand() % population.size();
    int l = rand() % population.size();
    while(k == l) {
        l = rand() % population.size();
    }

    //population[k] - parent 1
    stripX = population[k];
    stripSumX = populationSum[k];

    //population[l] - parent 2
    stripY = population[l];
    stripSumY = populationSum[l];

    double popKFit = Fitness(stripSumX, stripX);
    double popLFit = Fitness(stripSumY, stripY);

    //If GGA operator is chosen
    if(xOver == 1) {
        GGA(numScores, allScores, partners, adjMatrix, itemWidths, itemNumbers, offspringSum, offspring, stripSumX, stripX, stripSumY, stripY);
    }
        //If AGX' operator chosen
    else if(xOver == 2) {
        AGX(numScores, allScores, partners, adjMatrix, itemWidths, itemNumbers, offspringSum, offspring, stripSumX, stripX, stripSumY, stripY);
    }
    //If GPF operator chosen
    else if(xOver == 3){
        AGXDash(numScores, allScores, partners, adjMatrix, itemWidths, itemNumbers, offspringSum, offspring, stripSumX, stripX, stripSumY, stripY);
    }

    double offspringFit = Fitness(offspringSum, offspring);

    //1. If pop[k] has more strips than pop[l], replace pop[k]
    if(population[k].size() > population[l].size()){
        population[k] = offspring;
        populationSum[k] = offspringSum;
        if(offspring.size() < bestSize){
            bestSize = offspring.size();
            bestFit = offspringFit;
            bestTime = timer.elapsed();
            timeLogStream << bestSize << "\t" << bestFit << "\t" << timer.elapsed() << std::endl;
        }
        else if(offspring.size() == bestSize){
            if(offspringFit > bestFit){
                bestSize = offspring.size();
                bestFit = offspringFit;
                timeLogStream << bestSize << "\t" << bestFit << "\t" << timer.elapsed() << std::endl;
            }
        }
    }
    //2. If pop[l] has more strips than pop[k], replace pop[l]
    else if(population[l].size() > population[k].size()){
        population[l] = offspring;
        populationSum[l] = offspringSum;
        if(offspring.size() < bestSize){
            bestSize = offspring.size();
            bestFit = offspringFit;
            bestTime = timer.elapsed();
            timeLogStream << bestSize << "\t" << bestFit << "\t" << timer.elapsed() << std::endl;
        }
        else if(offspring.size() == bestSize){
            if(offspringFit > bestFit){
                bestSize = offspring.size();
                bestFit = offspringFit;
                timeLogStream << bestSize << "\t" << bestFit << "\t" << timer.elapsed() << std::endl;
            }
        }
    }

    //3. If pop[k] and pop[l] have same number of strips, compare fitness values
    else if(population[k].size() == population[l].size()){
        //3a. If popKFit is lower than popLFit, replace pop[k]
        if(popKFit < popLFit){
            population[k] = offspring;
            populationSum[k] = offspringSum;
            if(offspring.size() < bestSize){
                bestSize = offspring.size();
                bestFit = offspringFit;
                bestTime = timer.elapsed();
                timeLogStream << bestSize << "\t" << bestFit << "\t" << timer.elapsed() << std::endl;
            }
            else if(offspring.size() == bestSize){
                if(offspringFit > bestFit){
                    bestSize = offspring.size();
                    bestFit = offspringFit;
                    timeLogStream << bestSize << "\t" << bestFit << "\t" << timer.elapsed() << std::endl;
                }
            }
        }
        //3b. If popLFit is lower than popKFit, replace pop[l]
        else if(popLFit < popKFit){
            population[l] = offspring;
            populationSum[l] = offspringSum;
            if(offspring.size() < bestSize){
                bestSize = offspring.size();
                bestFit = offspringFit;
                bestTime = timer.elapsed();
                timeLogStream << bestSize << "\t" << bestFit << "\t" << timer.elapsed() << std::endl;
            }
            else if(offspring.size() == bestSize){
                if(offspringFit > bestFit){
                    bestSize = offspring.size();
                    bestFit = offspringFit;
                    timeLogStream << bestSize << "\t" << bestFit << "\t" << timer.elapsed() << std::endl;
                }
            }
        }
        //3c. If popKFit and popLFit are equal, choose at random.
        else if(popKFit == popLFit){
            double r = double(rand()) / double(RAND_MAX);
            if(r < 0.5){ //replace pop[k]
                population[k] = offspring;
                populationSum[k] = offspringSum;
                if(offspring.size() < bestSize){
                    bestSize = offspring.size();
                    bestFit = offspringFit;
                    bestTime = timer.elapsed();
                    timeLogStream << bestSize << "\t" << bestFit << "\t" << timer.elapsed() << std::endl;
                }
                else if(offspring.size() == bestSize){
                    if(offspringFit > bestFit){
                        bestSize = offspring.size();
                        bestFit = offspringFit;
                        timeLogStream << bestSize << "\t" << bestFit << "\t" << timer.elapsed() << std::endl;
                    }
                }
            }
            else if(r >= 0.5){ //replace pop[l]
                population[l] = offspring;
                populationSum[l] = offspringSum;
                if(offspring.size() < bestSize){
                    bestSize = offspring.size();
                    bestFit = offspringFit;
                    bestTime = timer.elapsed();
                    timeLogStream << bestSize << "\t" << bestFit << "\t" << timer.elapsed() << std::endl;
                }
                else if(offspring.size() == bestSize){
                    if(offspringFit > bestFit){
                        bestSize = offspring.size();
                        bestFit = offspringFit;
                        timeLogStream << bestSize << "\t" << bestFit << "\t" << timer.elapsed() << std::endl;
                    }
                }
            }
        }

    }



}

void GGA(int numScores, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
         std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<int>& offspringSum,
         std::vector<std::vector<int> >& offspring, std::vector<int>& stripSumX, std::vector<std::vector<int> >& stripX,
         std::vector<int>& stripSumY, std::vector<std::vector<int> >& stripY) {

    int k, l;
    std::vector<int> checked(numScores, 0);
    std::vector<int> partialItem;

    PermuteStrips(stripY, stripSumY);


    //choose these randomly
    if(stripY.size() == 2) {
        k = 0;
        l = 1;
    }
    else {
        k = rand() % stripY.size();
        l = rand() % stripY.size();
        while(k >= l || (k == 0 && l == stripY.size() - 1)) {
            k = rand() % stripY.size();
            l = rand() % stripY.size();
        }
    }


    for(int i = k; i <= l; ++i) {
        for(int j = 0; j < stripY[i].size(); ++j) {
            checked[stripY[i][j]] = 1;
        }
    }

    for(int i = 0; i < stripX.size(); ++i) {
        for(int j = 0; j < stripX[i].size(); ++j) {
            if(checked[stripX[i][j]] == 1) {
                stripX.erase(stripX.begin() + i);
                stripSumX.erase(stripSumX.begin() + i);
                --i;
                break;
            }
        }
    }

    for(int i = k; i <= l; ++i) {
        stripX.push_back(stripY[i]);
        stripSumX.push_back(stripSumY[i]);
    }

    for(int i = 0; i < stripX.size(); ++i) {
        for(int j = 0; j < stripX[i].size(); ++j) {
            checked[stripX[i][j]] = 1;
        }
    }

    for(int i = 0; i < checked.size(); ++i) {
        if(checked[i] == 0) {
            partialItem.push_back(i);
        }
    }

    if(partialItem.empty()) {
        //std::cout << "partialItem vector is empty - no boxes are missing\n";
        offspring = stripX;
        offspringSum = stripSumX;
        Mutation(numScores, allScores, partners, adjMatrix, itemWidths, itemNumbers, offspringSum, offspring);
    }
    else if(partialItem.size() % 2 == 0) {
        stripY.clear();
        stripY.resize(partialItem.size() / 2);
        stripSumY.clear();
        stripSumY.resize(partialItem.size() / 2, 0); //make stripSumY have size partialItem.size()/2, where each element
        //is initialized as 0.

        MFFPlus(numScores, 3, allScores, partners, partialItem, adjMatrix, itemWidths, itemNumbers, stripSumY, stripY);

        LocalSearch(numScores, allScores, partners, adjMatrix, itemWidths, itemNumbers, offspringSum, offspring, stripSumX, stripX, stripSumY, stripY);

        Mutation(numScores, allScores, partners, adjMatrix, itemWidths, itemNumbers, offspringSum, offspring);

    }
    else {
        std::cerr << "[ERROR]: partialItem.size() is odd, not valid.\n";
        exit(1);
    }
}


void AGX(int numScores, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
         std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<int>& offspringSum,
         std::vector<std::vector<int> >& offspring, std::vector<int>& stripSumX, std::vector<std::vector<int> >& stripX,
         std::vector<int>& stripSumY, std::vector<std::vector<int> >& stripY) {

    /* Alternates between sets X and Y, chooses the strip that has the largest stripSum (i.e. that is the "fullest") */

    int k;
    double r = 2.0; //must be 2.0
    std::vector<int> checked(numScores, 0);
    std::vector<int> partialItem;

    int minSetSize = std::min(stripX.size(), stripY.size());
    --minSetSize;


    //If stripX and stripY each have a strip whose stripSums are equal and the largest of all other strips, choose between them randomly.
    if(*std::max_element(stripSumX.begin(), stripSumX.end()) == *std::max_element(stripSumY.begin(), stripSumY.end())) {
        r = double(rand()) / double(RAND_MAX);
    }

    /**Fullest strip is in stripX**/
    if(*std::max_element(stripSumX.begin(), stripSumX.end()) > *std::max_element(stripSumY.begin(), stripSumY.end()) || r < 0.5) {
        while(offspring.size() < minSetSize && !stripX.empty()) {
            k = std::distance(stripSumX.begin(), std::max_element(stripSumX.begin(), stripSumX.end()));

            //Mark the items in the chosen strip from stripX in the checked vector
            for(int j = 0; j < stripX[k].size(); ++j) {
                checked[stripX[k][j]] = 1;
            }

            //Put the chosen strip from stripX into offspring
            offspring.push_back(stripX[k]);
            offspringSum.push_back(stripSumX[k]);
            stripX.erase(stripX.begin() + k);
            stripSumX.erase(stripSumX.begin() + k);

            //Go through strips in stripY, delete the strips that contain any items that have been checked (i.e. that are in offspring)
            for(int i = 0; i < stripY.size(); ++i) {
                for(int j = 0; j < stripY[i].size(); ++j) {
                    if(checked[stripY[i][j]] == 1) {
                        stripY.erase(stripY.begin() + i);
                        stripSumY.erase(stripSumY.begin() + i);
                        --i;
                        break;
                    }
                }
            }

            if(offspring.size() >= minSetSize || stripY.empty()) {
                break;
            }

            //Now go to stripY and find the fullest strip
            k = std::distance(stripSumY.begin(), std::max_element(stripSumY.begin(), stripSumY.end()));
            for(int j = 0; j < stripY[k].size(); ++j) {
                checked[stripY[k][j]] = 1;
            }
            offspring.push_back(stripY[k]);
            offspringSum.push_back(stripSumY[k]);
            stripY.erase(stripY.begin() + k);
            stripSumY.erase(stripSumY.begin() + k);

            for(int i = 0; i < stripX.size(); ++i) {
                for(int j = 0; j < stripX[i].size(); ++j) {
                    if(checked[stripX[i][j]] == 1) {
                        stripX.erase(stripX.begin() + i);
                        stripSumX.erase(stripSumX.begin() + i);
                        --i;
                        break;
                    }
                }
            }
        } // End while
    } //End If fullest strip is in stripX

        /**Fullest strip is in stripY**/
    else if(*std::max_element(stripSumX.begin(), stripSumX.end()) < *std::max_element(stripSumY.begin(), stripSumY.end()) || r >= 0.5) {
        while(offspring.size() < minSetSize && !stripY.empty()) {
            k = std::distance(stripSumY.begin(), std::max_element(stripSumY.begin(), stripSumY.end()));

            //Mark the items in the chosen strip from stripY in the checked vector
            for(int j = 0; j < stripY[k].size(); ++j) {
                checked[stripY[k][j]] = 1;
            }

            //Put the chosen strip from stripY into offspring
            offspring.push_back(stripY[k]);
            offspringSum.push_back(stripSumY[k]);
            stripY.erase(stripY.begin() + k);
            stripSumY.erase(stripSumY.begin() + k);

            //Go through strips in stripX, delete the strips that contain any items that have been checked (i.e. that are in offspring)
            for(int i = 0; i < stripX.size(); ++i) {
                for(int j = 0; j < stripX[i].size(); ++j) {
                    if(checked[stripX[i][j]] == 1) {
                        stripX.erase(stripX.begin() + i);
                        stripSumX.erase(stripSumX.begin() + i);
                        --i;
                        break;
                    }
                }
            }

            if(offspring.size() >= minSetSize || stripX.empty()) {
                break;
            }

            //Now go to stripX and find the fullest strip
            k = std::distance(stripSumX.begin(), std::max_element(stripSumX.begin(), stripSumX.end()));
            for(int j = 0; j < stripX[k].size(); ++j) {
                checked[stripX[k][j]] = 1;
            }

            offspring.push_back(stripX[k]);
            offspringSum.push_back(stripSumX[k]);
            stripX.erase(stripX.begin() + k);
            stripSumX.erase(stripSumX.begin() + k);

            for(int i = 0; i < stripY.size(); ++i) {
                for(int j = 0; j < stripY[i].size(); ++j) {
                    if(checked[stripY[i][j]] == 1) {
                        stripY.erase(stripY.begin() + i);
                        stripSumY.erase(stripSumY.begin() + i);
                        --i;
                        break;
                    }
                }
            }

        } // End while

    } //End If fullest strip is in stripY


    for(int i = 0; i < checked.size(); ++i) {
        if(checked[i] == 0) {
            partialItem.push_back(i);
        }
    }

    if(partialItem.empty()) {
        Mutation(numScores, allScores, partners, adjMatrix, itemWidths, itemNumbers, offspringSum, offspring);
    }
    else if(partialItem.size() % 2 == 0) {
        //sort(partialItem.begin(), partialItem.end());
        stripX = offspring;
        stripSumX = offspringSum;
        offspring.clear();
        offspringSum.clear();
        stripY.clear();
        stripY.resize(partialItem.size() / 2);
        stripSumY.clear();
        stripSumY.resize(partialItem.size() / 2, 0);

        MFFPlus(numScores, 3, allScores, partners, partialItem, adjMatrix, itemWidths, itemNumbers, stripSumY, stripY);

        LocalSearch(numScores, allScores, partners, adjMatrix, itemWidths, itemNumbers, offspringSum, offspring, stripSumX, stripX, stripSumY, stripY);

        Mutation(numScores, allScores, partners, adjMatrix, itemWidths, itemNumbers, offspringSum, offspring);

    }
    else {
        std::cerr << "[ERROR]: partialItem.size() is odd, not valid.\n";
        exit(1);
    }

}

void AGXDash(int numScores, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
             std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<int>& offspringSum,
             std::vector<std::vector<int> >& offspring, std::vector<int>& stripSumX, std::vector<std::vector<int> >& stripX,
             std::vector<int>& stripSumY, std::vector<std::vector<int> >& stripY) {

    /* Alternates between sets X and Y, chooses the strip that has the most items packed onto it */
    int a, b;
    double r = 2.0; //must be 2.0
    bool chooseSX = false;
    std::vector<int> stripNumX, stripNumY;
    std::vector<int> checked(numScores, 0);
    std::vector<int> partialItem;

    int minSetSize = std::min(stripX.size(), stripY.size());
    --minSetSize;


    //Create new vectors that contain the size(number of score widths/items*2) of each strip in each set X and Y
    for(int i = 0; i < stripX.size(); ++i){
        stripNumX.push_back(stripX[i].size());
    }

    for(int i = 0; i < stripY.size(); ++i){
        stripNumY.push_back(stripY[i].size());
    }


    //If stripX and stripY each have a strip whose stripNums are equal and the largest of all other strips, compare their fullness (stripSums).
    if(*std::max_element(stripNumX.begin(), stripNumX.end()) == *std::max_element(stripNumY.begin(), stripNumY.end())) {
        a = std::distance(stripNumX.begin(), std::max_element(stripNumX.begin(), stripNumX.end()));
        b = std::distance(stripNumY.begin(), std::max_element(stripNumY.begin(), stripNumY.end()));
        //If the strips also have the same stripSums (fullness), choose between them randomly
        if(stripSumX[a] == stripSumY[b]){
            r = double(rand()) / double(RAND_MAX);
        }
        else if(stripSumX[a] > stripSumY[b]){
            chooseSX = true;
        }
        else if(stripSumX[a] < stripSumY[b]){
            chooseSX = false;
        }

    }

    /**Largest strip is in stripX**/
    if(*std::max_element(stripNumX.begin(), stripNumX.end()) > *std::max_element(stripNumY.begin(), stripNumY.end()) || r < 0.5 || (chooseSX == true)) {
        while(offspring.size() < minSetSize && !stripX.empty()) {
            a = std::distance(stripNumX.begin(), std::max_element(stripNumX.begin(), stripNumX.end()));

            //Mark the items in the chosen strip from stripX in the checked vector
            for(int j = 0; j < stripX[a].size(); ++j) {
                checked[stripX[a][j]] = 1;
            }

            //Put the chosen strip from stripX into offspring
            offspring.push_back(stripX[a]);
            offspringSum.push_back(stripSumX[a]);
            stripX.erase(stripX.begin() + a);
            stripSumX.erase(stripSumX.begin() + a);
            stripNumX.erase(stripNumX.begin() + a);

            //Go through strips in stripY, delete the strips that contain any items that have been checked (i.e. that are in offspring)
            for(int i = 0; i < stripY.size(); ++i) {
                for(int j = 0; j < stripY[i].size(); ++j) {
                    if(checked[stripY[i][j]] == 1) {
                        stripY.erase(stripY.begin() + i);
                        stripSumY.erase(stripSumY.begin() + i);
                        stripNumY.erase(stripNumY.begin() + i);
                        --i;
                        break;
                    }
                }
            }

            if(offspring.size() >= minSetSize || stripY.empty()) {
                break;
            }

            //Now go to stripY and find the largest strip
            b = std::distance(stripNumY.begin(), std::max_element(stripNumY.begin(), stripNumY.end()));

            for(int j = 0; j < stripY[b].size(); ++j) {
                checked[stripY[b][j]] = 1;
            }
            offspring.push_back(stripY[b]);
            offspringSum.push_back(stripSumY[b]);
            stripY.erase(stripY.begin() + b);
            stripSumY.erase(stripSumY.begin() + b);
            stripNumY.erase(stripNumY.begin() + b);

            for(int i = 0; i < stripX.size(); ++i) {
                for(int j = 0; j < stripX[i].size(); ++j) {
                    if(checked[stripX[i][j]] == 1) {
                        stripX.erase(stripX.begin() + i);
                        stripSumX.erase(stripSumX.begin() + i);
                        stripNumX.erase(stripNumX.begin() + i);
                        --i;
                        break;
                    }
                }
            }
        } // End while
    } //End If largest strip is in stripX

        /**Largest strip is in stripY**/
    else if(*std::max_element(stripNumX.begin(), stripNumX.end()) < *std::max_element(stripNumY.begin(), stripNumY.end()) || r >= 0.5) {
        while(offspring.size() < minSetSize && !stripY.empty()) {
            b = std::distance(stripNumY.begin(), std::max_element(stripNumY.begin(), stripNumY.end()));

            //Mark the items in the chosen strip from stripY in the checked vector
            for(int j = 0; j < stripY[b].size(); ++j) {
                checked[stripY[b][j]] = 1;
            }

            //Put the chosen strip from stripY into offspring
            offspring.push_back(stripY[b]);
            offspringSum.push_back(stripSumY[b]);
            stripY.erase(stripY.begin() + b);
            stripSumY.erase(stripSumY.begin() + b);
            stripNumY.erase(stripNumY.begin() + b);

            //Go through strips in stripX, delete the strips that contain any items that have been checked (i.e. that are in offspring)
            for(int i = 0; i < stripX.size(); ++i) {
                for(int j = 0; j < stripX[i].size(); ++j) {
                    if(checked[stripX[i][j]] == 1) {
                        stripX.erase(stripX.begin() + i);
                        stripSumX.erase(stripSumX.begin() + i);
                        stripNumX.erase(stripNumX.begin() + i);
                        --i;
                        break;
                    }
                }
            }

            if(offspring.size() >= minSetSize || stripX.empty()) {
                break;
            }

            //Now go to stripX and find the largest strip
            a = std::distance(stripNumX.begin(), std::max_element(stripNumX.begin(), stripNumX.end()));

            for(int j = 0; j < stripX[a].size(); ++j) {
                checked[stripX[a][j]] = 1;
            }

            offspring.push_back(stripX[a]);
            offspringSum.push_back(stripSumX[a]);
            stripX.erase(stripX.begin() + a);
            stripSumX.erase(stripSumX.begin() + a);
            stripNumX.erase(stripNumX.begin() + a);

            for(int i = 0; i < stripY.size(); ++i) {
                for(int j = 0; j < stripY[i].size(); ++j) {
                    if(checked[stripY[i][j]] == 1) {
                        stripY.erase(stripY.begin() + i);
                        stripSumY.erase(stripSumY.begin() + i);
                        stripNumY.erase(stripNumY.begin() + i);
                        --i;
                        break;
                    }
                }
            }

        } // End while

    } //End If largest strip is in stripY


    for(int i = 0; i < checked.size(); ++i) {
        if(checked[i] == 0) {
            partialItem.push_back(i);
        }
    }

    if(partialItem.empty()) {
        Mutation(numScores, allScores, partners, adjMatrix, itemWidths, itemNumbers, offspringSum, offspring);
    }
    else if(partialItem.size() % 2 == 0) {
        //sort(partialItem.begin(), partialItem.end());
        stripX = offspring;
        stripSumX = offspringSum;
        offspring.clear();
        offspringSum.clear();
        stripY.clear();
        stripY.resize(partialItem.size() / 2);
        stripSumY.clear();
        stripSumY.resize(partialItem.size() / 2, 0);

        MFFPlus(numScores, 3, allScores, partners, partialItem, adjMatrix, itemWidths, itemNumbers, stripSumY, stripY);

        LocalSearch(numScores, allScores, partners, adjMatrix, itemWidths, itemNumbers, offspringSum, offspring, stripSumX, stripX, stripSumY, stripY);

        Mutation(numScores, allScores, partners, adjMatrix, itemWidths, itemNumbers, offspringSum, offspring);

    }
    else {
        std::cerr << "[ERROR]: partialItem.size() is odd, not valid.\n";
        exit(1);
    }

}






















































































