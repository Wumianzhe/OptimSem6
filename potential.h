#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include "matrix.h"
#include "vector.h"

struct transpTask {
    Matrix C;
    Matrix X = Matrix(0,0);
    vector_t a;
    vector_t b;
    std::pair<int,std::pair<int,int>> dMaxCell = {0, {0,0}};
    transpTask(Matrix& _C, vector_t& _a, vector_t& _b) : C(_C), a(_a), b(_b),X(_a.size(),_b.size(), -1) {}
    void buildInit();
    void solvePot();
    bool isOptimal();
    std::vector<std::pair<int,int>> buildCycle();
};

transpTask readTransport(std::string filename);


#endif // POTENTIAL_H_
