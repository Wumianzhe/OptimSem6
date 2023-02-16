#ifndef VECTOR_H_
#define VECTOR_H_
#include "matrix.h"

struct vector_t : public matrix_t {
    vector_t(int h) : matrix_t(h, 1){};
    inline double& operator[](const int i) { return el(i, 1); }
    inline double operator[](const int i) const { return el(i, 1); }
    double eqNorm() const;
};
double norm(const vector_t& v);
#endif // VECTOR_H_
