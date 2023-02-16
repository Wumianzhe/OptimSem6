#include "vector.h"

double vector_t::eqNorm() const {
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

double norm(const vector_t& v) { return v.eqNorm(); }
