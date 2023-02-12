#include <cmath>
#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>
#ifndef MATRIX_H_
#define MATRIX_H_

namespace mat {
double random(double a, double b);
}

struct matrix_t {
    int rows;
    int cols;

    ~matrix_t() {}
    matrix_t(const matrix_t& M);
    matrix_t(int h, int w = 1);
    inline double& el(int i, int j) { return _data[i * cols + j]; }
    inline double el(int i, int j) const { return _data[i * cols + j]; }
    inline double& operator()(int i, int j) { return el(i, j); }
    inline double operator()(int i, int j) const { return el(i, j); }
    inline double& operator[](const int i) { return _data[i]; }
    inline double operator[](const int i) const { return _data[i]; }
    matrix_t operator*(const matrix_t R) const;
    matrix_t operator*(const double n) const;
    matrix_t operator/(const double n) const;
    matrix_t operator-(const matrix_t R) const;
    matrix_t operator+(const matrix_t R) const;
    void operator=(const matrix_t& R);
    matrix_t T();
    double norm() const;
    double infnorm() const;
    inline void write(std::ofstream& fstr) const {
        fstr.write(reinterpret_cast<const char*>(_data.data()), _data.size() * sizeof(double));
    }

  private:
    std::vector<double> _data;
};
matrix_t eyes(int size);
std::ostream& operator<<(std::ostream& os, const matrix_t& M);
double norm(const matrix_t& M);
double infnorm(const matrix_t& M);
void matPrint(const matrix_t& M);
matrix_t ThomasAlg(const matrix_t& B);

#endif // MATRIX_H_
