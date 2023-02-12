#include "matrix.h"
#include <cstdio>
#include <cstdlib>

matrix_t::matrix_t(const matrix_t& M) {
    rows = M.rows;
    cols = M.cols;
    _data = M._data;
}
matrix_t::matrix_t(int h, int w) : _data(h * w, 0) {
    rows = h;
    cols = w;
    // _data = std::vector<double>(h * w);
}

matrix_t matrix_t::operator*(const matrix_t R) const {
    if (R.rows != cols) {
        throw std::runtime_error("matrix size mismatch");
    }
    matrix_t res(this->rows, R.cols);
    for (int i = 0; i < res.rows; i++) {
        for (int j = 0; j < res.cols; j++) {
            res(i, j) = 0;
            for (int k = 0; k < R.rows; k++) {
                res(i, j) += el(i, k) * R(k, j);
            }
        }
    }
    return res;
}

matrix_t matrix_t::operator*(const double n) const {
    matrix_t res(this->rows, this->cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            res(i, j) = n * el(i, j);
        }
    }
    return res;
}
matrix_t matrix_t::operator/(const double n) const {
    matrix_t res(this->rows, this->cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            res(i, j) = el(i, j) / n;
        }
    }
    return res;
}

matrix_t matrix_t::operator-(const matrix_t R) const {
    if (R.rows != rows || R.cols != cols) {
        throw std::runtime_error("matrix size mismatch");
    }
    matrix_t res(this->rows, this->cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            res(i, j) = el(i, j) - R(i, j);
        }
    }
    return res;
}
matrix_t matrix_t::operator+(const matrix_t R) const {
    if (R.rows != rows || R.cols != cols) {
        throw std::runtime_error("matrix size mismatch");
    }
    matrix_t res(this->rows, this->cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            res(i, j) = el(i, j) + R(i, j);
        }
    }
    return res;
}

void matrix_t::operator=(const matrix_t& R) {
    if (this == &R) {
        return;
    }
    this->_data = R._data;
}

matrix_t matrix_t::T() {
    matrix_t res(this->cols, this->rows);
    for (int i = 0; i < this->cols; i++) {
        for (int j = 0; j < this->rows; j++) {
            res(i, j) = (*this)(j, i);
        }
    }
    return res;
}

double matrix_t::norm() const {
    double res = 0;
    if (this->rows == 1) {
        for (int i = 0; i < this->cols; i++) {
            res += el(0, i) * el(0, i);
        }
    } else {
        for (int i = 0; i < this->rows; i++) {
            res += el(i, 0) * el(i, 0);
        }
    }
    res = sqrt(res);
    return res;
}

double matrix_t::infnorm() const {
    double res = 0;
    for (int i = 0; i < rows; i++) {
        double s = 0;
        for (int j = 0; j < cols; j++) {
            s += std::abs(el(i, j));
        }
        if (s > res) {
            res = s;
        }
    }
    return res;
}

matrix_t eyes(int size) {
    matrix_t mat(size, size);
    for (int i = 0; i < mat.rows; i++) {
        mat(i, i) = 1;
    }
    return mat;
}

std::ostream& operator<<(std::ostream& os, const matrix_t& M) {
    for (int i = 0; i < M.rows; i++) {
        for (int j = 0; j < M.cols; j++) {
            os << M.el(i, j) << ' ';
        }
        os << std::endl;
    }
    return os;
}

double mat::random(double a, double b) { return a + (double)rand() * (b - a) / RAND_MAX; }

double norm(const matrix_t& M) { return M.norm(); }
double infnorm(const matrix_t& M) { return M.infnorm(); }

void matPrint(const matrix_t& M) {
    fflush(stdout);
    std::printf("\n");
    for (int i = 0; i < M.rows; i++) {
        for (int j = 0; j < M.cols; j++) {
            std::printf("%0.6e ", M(i, j));
        }
        std::printf("\n");
    }
    fflush(stdout);
}

matrix_t ThomasAlg(const matrix_t& B) {
    int size = B.rows;
    matrix_t delta(size), lambda(size);
    delta[0] = -B(0, 2) / B(0, 1);
    lambda[0] = B(0, 3) / B(0, 1);
    for (int i = 1; i < size; i++) {
        delta[i] = -B(i, 2) / (B(i, 1) + B(i, 0) * delta[i - 1]);
        lambda[i] = (B(i, 3) - B(i, 0) * lambda[i - 1]) / (B(i, 1) + B(i, 0) * delta[i - 1]);
    }
    matrix_t M(size);
    M[size - 1] = lambda[size - 1];
    for (int i = size - 2; i >= 0; i--) {
        M[i] = delta[i] * M[i + 1] + lambda[i];
    }
    return M;
}
