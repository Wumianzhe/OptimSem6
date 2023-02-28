#ifndef VECTOR_H_
#define VECTOR_H_
#include "matrix.h"

struct vector_t : public Matrix {
    vector_t(int h) : Matrix(h, 1){};
    vector_t(Matrix& mat) : Matrix(mat) {
        rows *= cols;
        cols = 1;
    }
    vector_t(Matrix&& mat) : Matrix(mat) {
        rows *= cols;
        cols = 1;
    }
    inline double& operator[](const int i) { return el(i, 0); }
    inline double operator[](const int i) const { return el(i, 0); }
    bool operator>(const double R) const;
    bool operator<(const double R) const;
    bool operator<=(const double R) const;
    bool operator>=(const double R) const;
    void push_back(double el) {
        _data.push_back(el);
        rows++;
    }
    double eqNorm() const;
};
double norm(const vector_t& v);
double dot(const vector_t& r, const vector_t& l);
#endif // VECTOR_H_
