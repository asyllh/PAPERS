/*--------------/
ALH
dlx.h
CMSASCPP: Construct, Merge, Solve & Adapt (CMSA) algorithm for the SCPP
Article: Evolutionary Methods for the Score-Constrained Packing Problem, A.L. Hawa, R. Lewis, J.M. Thompson
17/05/2019
/--------------*/

#ifndef DLX_H
#define DLX_H

#include <iostream>
#include <cmath>
#include <climits>
#include <string>
#include <vector>

#include "consts.h"

class Node; // Declare Node struct so that it can be used in the Row struct below
class Column; //declare Column struct so that it can be used in the Node struct below

class Row {
public:
    int width; // The total width of all items in the row
    int number; // The number of the row in order of rows from text file, 0 to r-1
    int len; // The number of nodes (items) in the row
    Node* node; //The first node in the row, use it to find/list all nodes in the row

    Row() : width(0), number(0), len(0), node(nullptr) {}

    ~Row(){
        node = nullptr;
    }

};

class Node {
public:
    Node* left; // Predecessor in row (node before this current node in row)
    Node* right; // Successor in row (node after this current node in row)
    Node* up; // Predecessor in column (node above this current node in column)
    Node* down; // Successor in column (node below this current node in column)
    Column* col; // The column containing this current node
    Row* row; // The row containing this node

    Node() : left(nullptr), right(nullptr), up(nullptr), down(nullptr), col(nullptr), row(nullptr) {}

    ~Node(){
        if(left){ //using if(left) because some nodes (Node head in Column class) have left set to nullptr, so left->right cannot be accessed.
            left->right = right;
        }
        if(right){ //using if(right) because some nodes (Node head in Column class) have right set to nullptr, so right->left cannot be accessed.
            right->left = left;
        }
        up->down = down;
        down->up = up;
        left = nullptr;
        right = nullptr;
        up = nullptr;
        down = nullptr;
        col = nullptr;
        if(row){
            row->len--; //decrement the number of nodes in the row
            if(row->len == 0){ //if there are no more nodes in the row, delete the row
                delete row;
            }
        }
        row = nullptr;
    }

};

class Column {
public:
    int len; // The number of non-header items currently in this column's list
    std::string name; // The name of the column (name of item)
    Node head; // The list header
    Column* prev; // The column before this current column
    Column* next; // The column after this current column

    Column() : len(0), name(""), head(), prev(nullptr), next(nullptr) {}

    ~Column(){
        prev->next = next;
        next->prev = prev;
        prev = nullptr;
        next = nullptr;
    }

};


extern int numSolutions;
extern int bestLevel;
extern double bestFitXC;
extern Column* colRoot;
extern std::vector<Node*> choice;
extern std::vector<std::vector<std::string> > xSolution;

extern bool timeExceeded;
extern bool solutionFound;


void XC(std::vector<std::vector<int> >& SStar, std::vector<std::vector<int> >& SStarScores, std::vector<int>& SStarSum);

void StoreSolution();

double FitnessXC();

void Cover(Column* c);

void Uncover(Column* c);

void SelectBestColumn(Column*& bestCol);

void RecursiveSearch(const Timer& xctimer, int& level, Node*& currentNode, Column*& bestCol);




#endif
