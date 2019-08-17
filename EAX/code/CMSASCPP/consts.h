/*--------------/
ALH
consts.h
CMSASCPP: Construct, Merge, Solve & Adapt (CMSA) algorithm for the SCPP
Article: Evolutionary Methods for the Score-Constrained Packing Problem, A.L. Hawa, R. Lewis, J.M. Thompson
17/05/2019
/--------------*/

#ifndef CONSTS_H
#define CONSTS_H

#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <set>
#include <algorithm>

class Timer {
private:
    std::chrono::high_resolution_clock::time_point m_start;

public:
    Timer()
            : m_start(std::chrono::high_resolution_clock::now())
    {

    }

    void reset()
    {
        m_start = std::chrono::high_resolution_clock::now();
    }

    double elapsed() const
    {
        return std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - m_start).count();
    }

    /*~Timer()
    {
        std::cout << "Total Time: "
                  << std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(std::chrono::high_resolution_clock::now() - m_start).count()
                  << "seconds." << std::endl;
    }*/
};

struct Component {
    mutable int m_age; // Age of strip
    int m_width; //Total width of strip (sum of all itemWidths on strip)
    mutable bool m_inSolution; //mark as true if component (strip) is in the solution found by XC, otherwise false.
    std::vector<int> m_itemFeas; //will compare using this, itemNumbers
    std::vector<int> m_scoresFeas; //score numbers, indicies

    Component() { }

    Component(int a, int w, bool s, std::vector<int>& vec1, std::vector<int>& vec2)
            : m_age(a), m_width(w), m_inSolution(s), m_itemFeas(vec1), m_scoresFeas(vec2) {

    }

    Component(std::vector<int>& vec3)
            : m_itemFeas(vec3) {

    }

    friend std::ostream& operator<<(std::ostream& stream, const Component& component);
};

struct Order {
    bool operator()(const Component& lhs, const Component& rhs) const {
        return lhs.m_itemFeas < rhs.m_itemFeas;
    }
};

extern int instance;
extern int timeLimit;
extern int tau;
extern int numItem;
extern int stripWidth;
extern int numPop;
extern int maxAge;
extern int seed;
extern int instType;
extern int numIterations;
extern std::ofstream outputStream;

extern double bsfFit; //best so far Fitness (fitness of the best so far solution Sbsf)
extern std::set<Component, Order> BStar;
extern std::vector<std::vector<int> > Sbsf; //Solution (best so far)
extern std::vector<std::vector<int> > SbsfScores;
extern std::vector<int> SbsfSum;
extern std::ofstream timeLogStream;

extern std::ofstream xcLogStream;

extern int xcTimeLimit;

#endif
