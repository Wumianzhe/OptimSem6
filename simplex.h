#ifndef SIMPLEX_H_
#define SIMPLEX_H_
#include "matrix.h"
#include "vector.h"

enum extrem { min = -1, max = 1 };
struct task_t {
    vector_t C;
    Matrix A;
    vector_t b;
};

task_t genCanon(task_t task, std::set<int>& unbound, std::set<int>& noneq, int sign);
std::tuple<task_t, std::set<int>, std::set<int>> genDual(task_t task, std::set<int>& unbound, std::set<int>& noneq);
vector_t restore(vector_t& xCan, int size, std::set<int>& unbound);
vector_t simplex(task_t task, vector_t x0, std::set<int> Nk0);
vector_t initBasic(task_t task);
std::set<int> comb(std::set<int> L, int count, bool reset);
vector_t enumerate(task_t task);
Matrix solve(const Matrix& A,const Matrix& b);
double det(const Matrix& A);
#endif // SIMPLEX_H_
