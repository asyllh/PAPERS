/*--------------/
ALH
func.cpp
CMSASCPP: Construct, Merge, Solve & Adapt (CMSA) algorithm for the SCPP
Article: Evolutionary Methods for the Score-Constrained Packing Problem, A.L. Hawa, R. Lewis, J.M. Thompson
17/05/2019
/--------------*/

#include <algorithm>
#include <cmath>
#include <climits>
#include "func.h"

std::ostream& operator<<(std::ostream& stream, const Component& component) {

    stream << component.m_width << " ";
    for(int i = 0; i < component.m_itemFeas.size(); ++i) {
        stream << component.m_itemFeas[i] << " ";
    }
    return stream;
}

double bsfFit = 0.000000;
std::set<Component, Order> BStar;
std::vector<std::vector<int> > Sbsf;
std::vector<std::vector<int> > SbsfScores;
std::vector<int> SbsfSum;
std::ofstream timeLogStream; //stream to open file to log timestamps when better soln found

int LowerBound(double totalItemWidth) {
    int lBound = ceil(totalItemWidth / stripWidth);
    return lBound;
}

double Fitness(std::vector<std::vector<int> >& strip, std::vector<int>& stripSum) {
    double total = 0.0;

    for(int i = 0; i < strip.size(); ++i) {
        double a = static_cast<double>(stripSum[i]) / static_cast<double>(stripWidth);
        total += pow(a, 2);
    }
    double finalLong = total / static_cast<double>(strip.size());

    double final = (trunc(finalLong * 1000000)) / 1000000;

    return final;
}

void Construct(int numScores, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
                std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers) {


    for(int p = 0; p < numPop; ++p) {
        std::vector<std::vector<int> > strip(numItem);
        std::vector<int> stripSum(numItem, 0);
        /*strip.clear();
        strip.resize(numItem);
        stripSum.clear();
        stripSum.resize(numItem, 0);*/

        MFFPlus(numScores, allScores, partners, adjMatrix, itemWidths, itemNumbers, stripSum, strip);

        for(int i = 0; i < strip.size(); ++i){
            std::vector<int> vec;
            for(int j = 0; j < strip[i].size()-1; j +=2){
                vec.push_back(itemNumbers[strip[i][j]][strip[i][j+1]]);
            }
            std::sort(vec.begin(), vec.end());
            Component component(0, stripSum[i], false, vec, strip[i]);
            BStar.insert(component);
        }

    }

}

void Adapt(){

    std::set<Component, Order>::iterator it;
    for(it = BStar.begin(); it != BStar.end(); ){ //No ++it in for loop conditions, this is done in the body of the for loop
        if(it->m_inSolution == true){
            it->m_age = 0;
            it->m_inSolution = false;
            ++it;
        }
        else if(it->m_inSolution == false){
            if(it->m_age >= maxAge-1){
                it = BStar.erase(it);
            }
            else if(it->m_age >= 0){
                it->m_age++;
                ++it;
            }
        }
        else{
            std::cerr << "ERROR: m_age < 0." << std::endl;
            exit(1);
        }
    }

}

void MFFPlus(int numScores, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix, std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers, std::vector<int>& stripSum, std::vector<std::vector<int> >& strip) {

    //Modified First-Fit Plus (shell)

    int feasible;
    std::vector<int> itemOrder;
    std::vector<int> checked(numScores, 0);

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
                        overlap = false;
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

void CP(int nScores, int full, bool multiple, std::vector<int>& matchList, std::vector<int>& B, std::vector<int>& altHam,
        std::vector<std::vector<int> >& RStar){

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


void CMSA(const Timer& timer, int numScores, std::vector<int>& allScores, std::vector<int>& partners, std::vector<std::vector<int> >& adjMatrix,
          std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers){

    while(timer.elapsed() < timeLimit){
        //Construct function - creates numPop solutions using MFFR+, adds all strips to BStar with age = 0.
        Construct(numScores, allScores, partners, adjMatrix, itemWidths, itemNumbers);

        //std::cout << "BStar size: " << BStar.size() << std::endl;

        std::vector<std::vector<int> > SStar; //Optimal Solution to be produced using XC
        std::vector<std::vector<int> > SStarScores;
        std::vector<int> SStarSum;
        XC(SStar, SStarScores, SStarSum);



        if(solutionFound) {
            double starFit = Fitness(SStar, SStarSum);
            if(numIterations == 0){
                Sbsf = SStar;
                SbsfScores = SStarScores;
                SbsfSum = SStarSum;
                bsfFit = starFit;
                timeLogStream << numIterations << "\t" << Sbsf.size() << "\t" << bsfFit << "\t" << timer.elapsed() << std::endl;
            }
            //if solution found using XC better than the best solution found so far, update best solution.
            else if(SStar.size() < Sbsf.size()){
                Sbsf = SStar;
                SbsfScores = SStarScores;
                SbsfSum = SStarSum;
                bsfFit = starFit;
                timeLogStream << numIterations << "\t" << Sbsf.size() << "\t" << bsfFit << "\t" << timer.elapsed() << std::endl;
            }
            else if(SStar.size() == Sbsf.size()) {
                if(starFit > bsfFit) {
                    Sbsf = SStar;
                    SbsfScores = SStarScores;
                    SbsfSum = SStarSum;
                    bsfFit = starFit;
                    timeLogStream << numIterations << "\t" << Sbsf.size() << "\t" << bsfFit << "\t" << timer.elapsed() << std::endl;
                }
            }
        }
        //else if(!solutionFound) { do nothing }

        //Adapt function - reset age of components in optimal solution, increase age of all other components, remove components that have reached maxAge
        Adapt();

        //std::cout << "CMSA Time: " << timer.elapsed() << std::endl;

        ++numIterations;
    }

}//end CMSA



