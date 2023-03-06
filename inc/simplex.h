#ifndef SIMPLEX_H_
#define SIMPLEX_H_
#include "matrix.h"
#include "vector.h"

struct task_t {
    vector_t C;
    Matrix A;
    vector_t b;
};

void genTest();
vector_t simplex(task_t task, vector_t x0);
vector_t initBasic(task_t task);
std::set<int> comb(std::set<int> L, int count, bool reset);
vector_t enumerate(task_t task);
double det(Matrix A);
#endif // SIMPLEX_H_
