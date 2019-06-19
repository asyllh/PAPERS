/*--------------/
ALH
gfunc.cpp
InstGen: Problem Instance Generator for the SCPP
Article: Evolutionary Methods for the Score-Constrained Packing Problem, A.L. Hawa, R. Lewis, J.M. Thompson
11/01/2019
/--------------*/
#include "gfunc.h"

int instance; //Number of instances
int tau = 70; //Minimum scoring distance
int numItem = 100; //Number of items n in the set I
int minWidth = 1; //Minimum width of the scores
int maxWidth = 70; //Maximum width of the scores
int minItemWidth = 150; //Minimum width of the items
int maxItemWidth = 1000; //Maximum width of the items
int seed = 1;
int instType = 1; //Type of instance, 1 = artificial, 2 = real
int numTypes = 0;

void ProgramInfo() {

    std::cout << "Problem Instance Generator for the SCPP:\n-------------\n"
              << "PARAMETERS:\n"
              << "       -i <int>    [Instance number.]\n"
              << "       -t <int>    [Constraint value. Default = 70.]\n"
              << "       -n <int>    [Number of items. Default = 100.]\n"
              << "       -a <int>    [Minimum score width. Default = 1.]\n"
              << "       -b <int>    [Maximum score width. Default = 70.]\n"
              << "       -m <int>    [Minimum item width. Default = 150.]\n"
              << "       -M <int>    [Maximum item width. Default = 1000.]\n"
              << "       -s <int>    [Seed. Default = 1.]\n"
              << "       -c <int>    [Instance type. 1: Artificial. 2: Real. Default = 1.]\n"
              << "---------------\n\n";
}


void CreateAInstance(int numScores, double& totalItemWidth, std::vector<int>& allScores, std::vector<int>& partners,
                     std::vector<std::vector<int> >& adjMatrix, std::vector<std::vector<int> >& itemWidths,
                     std::vector<std::vector<int> >& itemNumbers) {

    std::vector<int> randOrder;
    std::vector<int> checkItem(numScores, 0);
    totalItemWidth = 0.0;

    //Create random values to be used as score widths, put in allScores vector (except last two elements)
    for(int i = 0; i < numScores; ++i) {
        allScores.push_back(rand() % (maxWidth - minWidth + 1) + minWidth);
    }

    //Sort all of the scores in the allScores vector in ascending order
    std::sort(allScores.begin(), allScores.end()); //sorts elements of vector in ascending order

    //Filling in adjacency matrix - if sum of two scores >= tau (70), then insert 1 into the matrix, else leave as 0
    for(int i = 0; i < allScores.size() - 1; ++i) {
        for(int j = i + 1; j < allScores.size(); ++j) {
            if(allScores[i] + allScores[j] >= tau) {
                adjMatrix[i][j] = 1;
                adjMatrix[j][i] = 1;
            }
        }
    }

    //Initially, randOrder vector will contain elements in the order 0, ..., numScores -2, numScores -1
    for(int i = 0; i < numScores; ++i) {
        randOrder.push_back(i);
    }

    //Randomly shuffle all values in randOrder vector
    //std::random_shuffle(randOrder.begin(), randOrder.end());
    std::shuffle(randOrder.begin(), randOrder.end(), std::default_random_engine(seed));

    //Assign partners to each score (i.e. pair up scores to define which scores are either side of the same box)
    //In the adjacency matrix, this will be represented by value 2
    //Therefore there will be a value of 2 in every row and every column, non repeating
    for(int i = 0; i < numItem; ++i) {
        adjMatrix[randOrder[2 * i]][randOrder[2 * i + 1]] = 2;
        adjMatrix[randOrder[2 * i + 1]][randOrder[2 * i]] = 2;
    }

    for(int i = 0; i < numScores; ++i) {
        for(int j = 0; j < numScores; ++j) {
            if(adjMatrix[i][j] == 2) {
                partners[i] = j;
                break;
            }
        }
    }

    int k = 0;
    for(int i = 0; i < numScores - 1; ++i) {
        for(int j = i + 1; j < numScores; ++j) {
            if(adjMatrix[i][j] == 2) {
                itemWidths[i][j] = rand() % (maxItemWidth - minItemWidth + 1) + minItemWidth;
                itemWidths[j][i] = itemWidths[i][j];
                itemNumbers[i][j] = k;
                itemNumbers[j][i] = k;
                //allWidths[k] = itemWidths[i][j];
                ++k;
                break;
            }
        }
    }


    //std::cout << std::right << std::setw(5) << "Item#" << std::setw(10) << "Scores" << std::setw(12) << "Partners" << std::setw(10) << "Width\n";
    int count = 1;
    for(int i = 0; i < numScores; ++i) {
        if(checkItem[i] == 1) {
            continue;
        }
        /*std::cout << count << std::setw(9) << allScores[i] << "-" << allScores[partners[i]]
        << std::setw(10) << i  << "-" << partners[i] << std::setw(10) << itemWidths[i][partners[i]] << std::endl;*/
        totalItemWidth += itemWidths[i][partners[i]];
        checkItem[i] = 1;
        checkItem[partners[i]] = 1;
        ++count;
    }
    //std::cout << "Total Item Widths: " << totalItemWidth << std::endl << std::endl;


}

void CreateRInstance(int numScores, double& totalItemWidth, std::vector<int>& allScores, std::vector<int>& partners,
                     std::vector<std::vector<int> >& adjMatrix, std::vector<std::vector<int> >& itemWidths,
                     std::vector<std::vector<int> >& itemNumbers, std::vector<std::vector<int> >& itemTypes, std::vector<int>& typeNumber) {

    totalItemWidth = 0.0;

    numTypes = rand() % 21 + 10; //Create number of item types between 10 and 30

    double sumRand = 0.0;
    std::vector<double> randNums;

    for(int i = 0; i < numTypes; ++i) {
        double r = ((double) rand()) / ((double) RAND_MAX) / 2.0 + 0.1;
        randNums.push_back(r);
        sumRand += r;
    }

    int sum = 0;
    double coeff = static_cast<double>(numItem) / sumRand;
    for(int i = 0; i < randNums.size() - 1; ++i) {
        int num = floor(coeff * randNums[i]);
        sum += num;
        itemTypes[0].push_back(num);
    }

    int last = numItem - sum;
    itemTypes[0].push_back(last);

    std::vector<int> allWidthsUnsort;
    std::vector<int> typeNumberUnsort;

    for(int j = 0; j < numTypes; ++j) {
        int r1 = rand() % (maxWidth - minWidth + 1) + minWidth;
        int r2 = rand() % (maxWidth - minWidth + 1) + minWidth;
        if(r1 <= r2) {
            itemTypes[1].push_back(r1);
            itemTypes[2].push_back(r2);
        }
        else {
            itemTypes[1].push_back(r2);
            itemTypes[2].push_back(r1);
        }
        int w = rand() % (maxItemWidth - minItemWidth + 1) + minItemWidth;
        itemTypes[3].push_back(w);
        for(int k = 0; k < itemTypes[0][j]; ++k) {
            allScores.push_back(itemTypes[1][j]);
            allScores.push_back(itemTypes[2][j]);
            allWidthsUnsort.push_back(itemTypes[3][j]);
            allWidthsUnsort.push_back(itemTypes[3][j]);
            typeNumberUnsort.push_back(j);
            typeNumberUnsort.push_back(j);
        }
    }

    std::vector<int> partUnsort(numScores);

    for(int i = 0; i < partUnsort.size() - 1; i += 2) {
        partUnsort[i] = i + 1;
        partUnsort[i + 1] = i;
    }

    std::vector<int> order;
    std::vector<int> invOrder(numScores);

    for(int i = 0; i < numScores; ++i) {
        order.push_back(i);
    }

    for(int i = 1; i < numScores; ++i) {
        for(int j = i - 1; j >= 0; --j) {
            if(allScores[i] < allScores[order[j]]) {
                order[j + 1] = order[j];
                order[j] = i;
            }
        }
    }

    for(int i = 0; i < numScores; ++i) {
        invOrder[order[i]] = i;
    }

    std::sort(allScores.begin(), allScores.end());

    std::vector<int> checked(numScores, 0);

    for(int a = 0; a < numScores; ++a) {
        if(checked[a] == 1) {
            continue;
        }
        partners[a] = invOrder[partUnsort[order[a]]];
        partners[partners[a]] = a;
        checked[a] = 1;
        checked[partners[a]] = 1;
    }
    checked.clear();
    checked.resize(numScores, 0);

    std::vector<int> allWidths(numScores);

    for(int i = 0; i < numScores; ++i) {
        allWidths[i] = allWidthsUnsort[order[i]];
        typeNumber[i] = typeNumberUnsort[order[i]];
    }

    for(int i = 0; i < allScores.size() - 1; ++i) {
        for(int j = i + 1; j < allScores.size(); ++j) {
            if(allScores[i] + allScores[j] >= tau) {
                adjMatrix[i][j] = 1;
                adjMatrix[j][i] = 1;
            }
        }
    }

    for(int i = 0; i < partners.size(); ++i) {
        if(checked[i] == 1) {
            continue;
        }
        adjMatrix[i][partners[i]] = 2;
        adjMatrix[partners[i]][i] = 2;
        checked[i] = 1;
        checked[partners[i]] = 1;
    }
    checked.clear();
    checked.resize(numScores, 0);

    int k = 0;
    for(int i = 0; i < numScores; ++i) {
        if(checked[i] == 1) {
            continue;
        }
        itemWidths[i][partners[i]] = allWidths[i];
        itemWidths[partners[i]][i] = allWidths[i];
        itemNumbers[i][partners[i]] = k;
        itemNumbers[partners[i]][i] = k;
        checked[i] = 1;
        checked[partners[i]] = 1;
        ++k;
    }

    for(int j = 0; j < itemTypes[0].size(); ++j) {
        int sum = itemTypes[0][j] * itemTypes[3][j];
        totalItemWidth += sum;
    }

}

void OutputProbInst(int numScores, double totalItemWidth, std::vector<int>& allScores, std::vector<int>& partners,
                    std::vector<std::vector<int> >& itemWidths, std::vector<std::vector<int> >& itemNumbers,
                    std::vector<std::vector<int> >& itemTypes, std::vector<int>& typeNumber) {

    std::string filename;
    int n = numItem / 100;
    if(instType == 1) {
        filename = "PA" + std::to_string(n) + "_" + std::to_string(instance) + ".txt";
    }
    else {
        filename = "PR" + std::to_string(n) + "_" + std::to_string(instance) + ".txt";
    }

    std::ofstream ofs(filename.c_str());
    if(!ofs) {
        std::cerr << "[ERROR]: Cannot write to file." << std::endl;
        exit(1);
    }


    ofs << instance << std::endl; //Instance number
    ofs << numItem << std::endl; //Number of items
    ofs << instType << std::endl; //Artificial or real instance type

    for(const auto& v : allScores) {
        ofs << v << " ";
    }
    ofs << std::endl;

    for(const auto& v : partners) {
        ofs << v << " ";
    }
    ofs << std::endl;

    for(int i = 0; i < numScores; ++i) {
        ofs << itemWidths[i][partners[i]] << " ";
    }
    ofs << std::endl;

    for(int i = 0; i < numScores; ++i) {
        ofs << itemNumbers[i][partners[i]] << " ";
    }
    ofs << std::endl;

    ofs << totalItemWidth << std::endl;

    //For real instances only
    if(instType == 2) {
        ofs << numTypes << std::endl;

        for(const auto& subVec : itemTypes) {
            for(const auto& v : subVec) {
                ofs << v << " ";
            }
            ofs << std::endl;
        }

        for(const auto& v : typeNumber) {
            ofs << v << " ";
        }
    }

    ofs.close();

}

