/*--------------/
ALH
dlx.cpp
CMSASCPP: Construct, Merge, Solve & Adapt (CMSA) algorithm for the SCPP
Article: Evolutionary Methods for the Score-Constrained Packing Problem, A.L. Hawa, R. Lewis, J.M. Thompson
17/05/2019
/--------------*/

#include "dlx.h"

int numSolutions = 0;
int bestLevel = RAND_MAX;
double bestFitXC = 0.000000;
Column* colRoot = nullptr;
std::vector<Node*> choice;
std::vector<std::vector<std::string> > xSolution;
bool timeExceeded = false;
bool solutionFound = false;
std::ofstream xcLogStream;
int xcTimeLimit = 600;

void XC(std::vector<std::vector<int> >& SStar, std::vector<std::vector<int> >& SStarScores, std::vector<int>& SStarSum){

    Column* currentCol;
    Column* prevCol = new Column;
    colRoot = prevCol;
    for(int i = 0; i < numItem; ++i){
        currentCol = new Column;
        currentCol->name = std::to_string(i);
        currentCol->head.up = currentCol->head.down = &currentCol->head;
        currentCol->len = 0;
        currentCol->prev = prevCol;
        prevCol->next = currentCol;
        prevCol = currentCol;
    }
    prevCol->next = colRoot;
    colRoot->prev = prevCol;
    currentCol = nullptr;
    prevCol = nullptr;

    int count = 1;
    Node* currentNode;
    std::set<Component, Order>::iterator it;
    for(it = BStar.begin(); it != BStar.end(); ++it){
        Node* rowStart = nullptr;
        Column* ccol = nullptr;
        Node* prevNode = nullptr;
        Row* currentRow = new Row;
        std::vector<int> tempIN = it->m_itemFeas;
        currentRow->width = it->m_width;
        currentRow->number = count;
        for(int j = 0; j < tempIN.size(); ++j){
            std::string tempString = std::to_string(tempIN[j]);
            for(ccol = colRoot->next; tempString != ccol->name && ccol != colRoot; ccol = ccol->next);
            currentNode = new Node;
            if(!rowStart) {
                rowStart = currentNode;
                currentRow->node = currentNode;
            }
            else {
                currentNode->left = prevNode;
                prevNode->right = currentNode;
            }
            currentNode->col = ccol;

            currentNode->up = ccol->head.up, ccol->head.up->down = currentNode;
            ccol->head.up = currentNode, currentNode->down = &ccol->head;
            currentNode->row = currentRow;
            ++ccol->len;
            ++currentRow->len;
            prevNode = currentNode;
        }
        rowStart->left = prevNode;
        prevNode->right = rowStart;
        ++count;
    }
    currentNode = nullptr;

    Column* bestCol = nullptr;

    int level = 0;
    solutionFound = false;

    Timer xctimer;
    RecursiveSearch(xctimer, level, currentNode, bestCol);

    xcLogStream << numIterations << "\t" << BStar.size() << "\t" << bestLevel << "\t" << bestFitXC << "\t" << numSolutions << "\t" << xctimer.elapsed() << std::endl;

    //itt = iterator for 'true' (solutionFound == true), itf = iterator for 'false' (solutionFound == false).
    if(solutionFound) {
        std::vector<int> starStrip;
        for(auto& vec : xSolution) {
            for(auto& element : vec) {
                int num = std::stoi(element);
                starStrip.push_back(num);
            }
            std::sort(starStrip.begin(), starStrip.end());
            std::set<Component, Order>::iterator itt = BStar.find(Component(starStrip));
            if(itt != BStar.end()) {
                itt->m_inSolution = true;
                SStar.push_back(starStrip);
                SStarScores.push_back(itt->m_scoresFeas);
                SStarSum.push_back(itt->m_width);
            }
            else if(itt == BStar.end()) {
                std::cerr << "ERROR: itemFeas vector not found in BStar." << std::endl;
                exit(1);
            }
            starStrip.clear();
        }
    }
    else if(!solutionFound){ //If no solution found, Sbsf becomes the current solution
        //Need to mark strips in Sbsf as inSolution
        for(int i = 0; i < Sbsf.size(); ++i){
            std::set<Component, Order>::iterator itf = BStar.find(Component(Sbsf[i]));
            if(itf != BStar.end()){
                itf->m_inSolution = true;
            }
            else if(itf == BStar.end()){ //If the strip isn't in BStar, that means it must have been removed from BStar at some point due to reaching maxAge
                //Therefore we must reinsert the strip back into BStar, and set it as being in the current best so far solution
                Component component(0, SbsfSum[i], true, Sbsf[i], SbsfScores[i]);

            }
        }
    }

    //Deleting circularly linked lists:
    //For each column, starting at the head node of the column, delete every node in the column (excluding head node)
    for(Column* currCol = colRoot->next; currCol != colRoot; currCol = currCol->next){
        Node* currNode = &(currCol->head);
        while(currNode->down != &(currCol->head)){
            delete currNode->down;
        }
    }

    //Delete every column to the right of colRoot until only colRoot is left
    while(colRoot->next != colRoot){
        delete colRoot->next;
    }

    numSolutions = 0;
    bestLevel = RAND_MAX;
    bestFitXC = 0.000000;
    colRoot = nullptr;
    choice.clear();
    xSolution.clear();
    timeExceeded = false;


} //End of XC function

void StoreSolution(){

    xSolution.clear();
    std::vector<std::string> temp;
    std::vector<Node*>::iterator it;

    for(it = choice.begin(); it != choice.end(); ++it){
        Node* first = *it;
        Node* n = first;
        do{
            temp.push_back(n->col->name);
            n = n->right;
        } while (n != first);
        xSolution.push_back(temp);
        temp.clear();
    }
}

double FitnessXC(){

    double total = 0.0;
    std::vector<Node*>::iterator it;
    for(it = choice.begin(); it != choice.end(); ++it){
        Node* n = *it;
        double a = static_cast<double>(n->row->width) / static_cast<double>(stripWidth);
        total += pow(a, 2);
    }

    double finalLong = total / static_cast<double>(choice.size());

    double final = (trunc(finalLong * 1000000)) / 1000000;

    return final;
}

//Cover column, block row, leaves all columns except column that is being covered, so a node is never removed from a list twice.
void Cover(Column* c){

    Node* nn = nullptr;
    Node* rr = nullptr;
    Node* uu = nullptr;
    Node* dd = nullptr;
    Column* l = nullptr;
    Column* r = nullptr;

    l = c->prev;
    r = c->next;
    l->next = r;
    r->prev = l;

    for(rr = c->head.down; rr != &(c->head); rr = rr->down){
        for(nn = rr->right; nn != rr; nn = nn->right){
            uu = nn->up;
            dd = nn->down;
            uu->down = dd;
            dd->up = uu;
            --nn->col->len;
        }
    }
}

//Uncovering a column, done in exact reverse order of covering.
void Uncover(Column* c){

    Node* nn = nullptr;
    Node* rr = nullptr;
    Node* uu = nullptr;
    Node* dd = nullptr;
    Column* l = nullptr;
    Column* r = nullptr;

    for(rr = c->head.up; rr!= &(c->head); rr = rr->up){
        for(nn = rr->left; nn != rr; nn = nn->left){
            uu = nn->up;
            dd = nn->down;
            uu->down = dd->up = nn;
            ++nn->col->len;
        }
    }
    l = c->prev;
    r = c->next;
    l->next = r->prev = c;
}

//Function to choose column object c = 'bestCol', should return pointer to bestCol only.
void SelectBestColumn(Column*& bestCol){

    long int minLen = LONG_MAX; //is this large enough?
    for(Column* currCol = colRoot->next; currCol != colRoot; currCol = currCol->next){
        if(currCol->len < minLen){
            bestCol = currCol, minLen = currCol->len;
        }
    }
}

void RecursiveSearch(const Timer& xctimer, int& level, Node*& currentNode, Column*& bestCol){
    SelectBestColumn(bestCol); // Returns bestCol pointer (line 2)
    Cover(bestCol); // Cover bestCol column (line 3)
    if(numSolutions == 0 && choice.size() <= level){ //Set r <- D[c] and O_k <- r, starting from first node below head node of column (line 4/5)
        choice.push_back(bestCol->head.down);
        currentNode = bestCol->head.down;
    }
    else{ currentNode = choice[level] = bestCol->head.down; }

    while(currentNode != &(bestCol->head)){ // While r != c, continue going down column until head node is reached (line 4)
        for(Node* rowNode = currentNode->right; rowNode != currentNode; rowNode = rowNode->right){ // For each j <- R[r] ... while j!=r, cover column j (line 6/7)
            Cover(rowNode->col);
        }
        if(xctimer.elapsed() > xcTimeLimit){
            if(numIterations == 0 && solutionFound == false){
                timeExceeded = false;
            }
            else {
                timeExceeded = true;
            }
        }
        if(colRoot->next != colRoot && timeExceeded == false){ // Do search(k+1) if root is not the only column left
            if(level < bestLevel){
                ++level;
                RecursiveSearch(xctimer, level, currentNode, bestCol);
            }
        }
        else if(colRoot->next == colRoot){
            if(level < bestLevel){
                ++numSolutions;
                solutionFound = true;
                bestLevel = level;
                choice.resize(bestLevel+1);
                StoreSolution();
                bestFitXC = FitnessXC();
            }
            else if(level == bestLevel){
                double currentFit = FitnessXC();
                if(currentFit >= bestFitXC){
                    ++numSolutions;
                    solutionFound = true;
                    StoreSolution();
                    bestFitXC = currentFit;
                }
            }
            else{ std::cerr << "[ERROR]: level > bestLevel." << std::endl; exit(1); }
        }
        for(Node* rowNode = currentNode->left; rowNode != currentNode; rowNode = rowNode->left){ // For each j <- L[r],... while j!=r, uncover column j (line 10/11)
            Uncover(rowNode->col);
        }
        if(timeExceeded){
            break; //break out of while loop
        }
        else { currentNode = choice[level] = currentNode->down; }// Set currentNode to the next node down in the column
    }
    //Exit while loop when all nodes in column have been assessed, currentNode == head node of bestCol
    Uncover(bestCol);
    if(level > 0){ // Go up a level, exit this search function and go back into the previous search function
        --level;
        currentNode = choice[level];
        bestCol = currentNode->col;
    }
    else{ // level == 0, top level reached, end function
        return; //is this needed??
    } // This code will end here if level == 0, or will go back into previous recursiveSearch function if level > 0.

} //End recursive search function
