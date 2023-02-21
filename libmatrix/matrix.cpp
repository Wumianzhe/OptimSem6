#include "matrix.h"
#include "vector.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>

Matrix::Matrix(const Matrix& M) {
    rows = M.rows;
    cols = M.cols;
    _data = M._data;
}
Matrix::Matrix(int h, int w) : _data(h * w, 0) {
    rows = h;
    cols = w;
    // _data = std::vector<double>(h * w);
}
auto Matrix::colSpan(int index) const {
    auto begin = _data.cbegin() + index * rows;
    auto end = begin + rows;
    return std::ranges::subrange(begin, end);
}
Matrix::Matrix(const Matrix& M, std::set<int> colIndices) {
    // whether all indices are inbound
    if ((*colIndices.begin() < 0) || (*colIndices.rbegin() >= M.cols)) {
        throw std::out_of_range("column index out of range");
    }
    cols = colIndices.size();
    rows = M.rows;
    for (int index : colIndices) {
        auto range = M.colSpan(index);
        std::ranges::copy(range, std::back_inserter(_data));
    }
}
Matrix Matrix::operator*(const Matrix R) const {
    if (R.rows != cols) {
        throw std::runtime_error("matrix size mismatch");
    }
    Matrix res(this->rows, R.cols);
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

Matrix Matrix::operator*(const double n) const {
    Matrix res(this->rows, this->cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            res(i, j) = n * el(i, j);
        }
    }
    return res;
}
Matrix Matrix::operator/(const double n) const {
    Matrix res(this->rows, this->cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            res(i, j) = el(i, j) / n;
        }
    }
    return res;
}

Matrix Matrix::operator-(const Matrix R) const {
    if (R.rows != rows || R.cols != cols) {
        throw std::runtime_error("matrix size mismatch");
    }
    Matrix res(this->rows, this->cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            res(i, j) = el(i, j) - R(i, j);
        }
    }
    return res;
}
Matrix Matrix::operator+(const Matrix R) const {
    if (R.rows != rows || R.cols != cols) {
        throw std::runtime_error("matrix size mismatch");
    }
    Matrix res(this->rows, this->cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            res(i, j) = el(i, j) + R(i, j);
        }
    }
    return res;
}

void Matrix::operator=(const Matrix& R) {
    if (this == &R) {
        return;
    }
    this->rows = R.rows;
    this->cols = R.cols;
    this->_data = R._data;
}

Matrix Matrix::operator[](const std::set<int> N) const { return Matrix(*this, N); }
Matrix Matrix::T() {
    Matrix res(this->cols, this->rows);
    for (int i = 0; i < this->cols; i++) {
        for (int j = 0; j < this->rows; j++) {
            res(i, j) = el(j, i);
        }
    }
    return res;
}

double Matrix::infnorm() const {
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

Matrix Matrix::eyes(int size) {
    Matrix mat(size, size);
    for (int i = 0; i < mat.rows; i++) {
        mat(i, i) = 1;
    }
    return mat;
}

std::ostream& operator<<(std::ostream& os, const Matrix& M) {
    for (int i = 0; i < M.rows; i++) {
        for (int j = 0; j < M.cols; j++) {
            os << M.el(i, j) << ' ';
        }
        os << std::endl;
    }
    return os;
}

double mat::random(double a, double b) { return a + (double)rand() * (b - a) / RAND_MAX; }

double infnorm(const Matrix& M) { return M.infnorm(); }

void matPrint(const Matrix& M) {
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

Matrix Matrix::ThomasAlg(const Matrix& B) {
    int size = B.rows;
    vector_t delta(size), lambda(size);
    delta[0] = -B(0, 2) / B(0, 1);
    lambda[0] = B(0, 3) / B(0, 1);
    for (int i = 1; i < size; i++) {
        delta[i] = -B(i, 2) / (B(i, 1) + B(i, 0) * delta[i - 1]);
        lambda[i] = (B(i, 3) - B(i, 0) * lambda[i - 1]) / (B(i, 1) + B(i, 0) * delta[i - 1]);
    }
    vector_t M(size);
    M[size - 1] = lambda[size - 1];
    for (int i = size - 2; i >= 0; i--) {
        M[i] = delta[i] * M[i + 1] + lambda[i];
    }
    return M;
}
