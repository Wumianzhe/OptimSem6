#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include "matrix.h"
#include "vector.h"

struct transpTask {
    Matrix C;
    Matrix X = Matrix(0,0);
    vector_t a;
    vector_t b;
    int diff = 0;
    std::pair<int,std::pair<int,int>> dMaxCell = {0, {0,0}};
    transpTask(Matrix& _C, vector_t& _a, vector_t& _b, int _diff) : C(_C), a(_a), b(_b),X(_a.size(),_b.size(), -1),diff(_diff) {}
    void buildInit();
    void solvePot();
    bool isOptimal();
    int results(const Matrix& Xr);
    void prettyPrint();
    std::vector<std::pair<int,int>> buildCycle();
};

transpTask readTransport(std::string filename);


#endif // POTENTIAL_H_
