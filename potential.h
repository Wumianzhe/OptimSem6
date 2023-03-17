#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include "matrix.h"
#include "vector.h"

struct transpTask {
    Matrix C;
    vector_t a;
    vector_t b;
    transpTask(Matrix& _C, vector_t& _a, vector_t& _b) : C(_C), a(_a), b(_b) {}
};

transpTask readTransport(std::string filename);

#endif // POTENTIAL_H_
