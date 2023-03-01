#include "simplex.h"
#include <algorithm>

using namespace std;

double eps = 1e-8;
struct {
    set<int> supSet;
    int size;
    vector<bool> mask;
    bool first;
    set<int> next() {
        if (first) {
            first = !first;
        } else {
            if (!prev_permutation(mask.begin(), mask.end())) {
                return {};
            }
        }
        set<int> perm;
        auto it = supSet.begin();
        for (int i = 0; i < mask.size(); i++) {
            if (mask[i]) {
                perm.insert(*it);
            }
            it = std::next(it);
        }
        return perm;
    }
    void reset(set<int> L, int size) {
        supSet = L;
        mask.resize(size, true);
        mask.resize(supSet.size(), false);
        first = true;
    }
} generator;

// should not alter initial objects
set<int> merge(set<int> N, set<int> L) {
    N.merge(L);
    return N;
}
vector_t genNext(task_t task, Matrix B, vector_t xk, int jk, set<int> Nk);
set<int> fill(Matrix& A, set<int> N);
double det(Matrix A);
Matrix solve(Matrix A, Matrix b, Matrix& Ai);

vector_t simplex(task_t task, vector_t x0) {
    int m = task.A.rows;
    int n = task.A.cols;
    set<int> Nkp;
    set<int> Nk;
    vector_t xk = x0;
    Matrix B = Matrix::eyes(m); // A^-1
    // 1.
    for (int i = 0; i < n; i++) {
        if (xk[i] != 0) {
            Nkp.insert(i);
        }
    }
    set<int> Lkp;
    for (int i = 0; i < n; i++) {
        if (Nkp.find(i) == Nkp.end()) {
            Lkp.insert(i);
        }
    }
    if (Nkp.size() < m) {
        generator.reset(Lkp, m - Nkp.size());
        Nk = fill(task.A, Nkp);
    } else {
        Nk = Nkp;
    }
    Matrix ykT = solve(task.A[Nk], task.b, B);
    // 2.
    while (true) {
        // Lk = N \ Nk (not Nkp!)
        set<int> Lk;
        for (int i = 0; i < n; i++) {
            if (Nk.find(i) == Nk.end()) {
                Lk.insert(i);
            }
        }
        // whole vector to avoid index confusion like in 3.b
        vector_t dk = task.C.T() - ykT * task.A;
        // 3.a
        auto it = find_if(dk.begin(), dk.end(), [&dk](int j) {
            // to avoid -0
            return (abs(dk[j]) > 1e-12) && (dk[j] < 0);
        });
        // dk >= 0
        if (it == dk.end()) {
            return xk;
        }
        int jk = *it;
        // 3.b
        vector_t uknk = B * task.A[{jk}];
        vector_t uk(n);
        // copy using indices I need
        int i = 0;
        for (int index : Nk) {
            uk[index] = uk[i++];
        }
        uk[jk] = -1; // rest is 0 by construction
        // 4.a
        set<int> P;
        // equivalent of P = { i in N | uk[i] > 0}
        copy_if(Nk.begin(), Nk.end(), inserter(P, P.end()), [&uk](int i) { return uk[i] > 0; });
        if (P.size() == 0) {
            throw logic_error("no lower bound");
        }
        // 4.b
        int ik = *P.begin();
        double thk = xk[ik] / uk[ik];
        for (int i : P) {
            if (xk[i] / uk[i] < thk) {
                ik = i;
                thk = xk[i] / uk[i];
            }
        }
        xk = xk - uk * thk;
        // 5
        if (abs(thk) > 1e-12) { // thk != 0
            Matrix F = Matrix::eyes(n);
            double div = uk[ik];
            for (int i = 0; i < n; i++) {
                if (i != ik) {
                    F(i, ik) = -uk[i] / div;
                } else {
                    F(i, i) = 1 / div;
                }
            }
            B = F * B;
            ykT = task.C.T() * B;
            // swapping indices and resetting generator
            Nk.erase(jk);
            Nk.insert(ik);
            Nkp.erase(jk);
            Nkp.insert(ik);
            generator.reset(Nkp, m - Nkp.size());
        } else {
            Nk = fill(task.A, Nkp);
            ykT = solve(task.A[Nk], task.b, B);
        }
    }
}

// finding starting basic vector
vector_t initBasic(task_t task) {
    int n = task.A.cols;
    int m = task.A.rows;
    vector_t secC(n + m);
    for (int i = 0; i < m; i++) {
        secC[n + i] = 1;
    }
    Matrix secA(m, n + m);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            secA(i, j) = task.A(i, j);
        }
        secA(i, n + i) = 1;
    }
    vector_t x0(n + m);
    for (int i = 0; i < m; i++) {
        if (task.b[i] >= 0) {
            x0[i] = task.b[i];
        } else {
            x0[i] = -task.b[i];
            for (int j = 0; j < n + m; j++) {
                secA(i, j) *= -1;
            }
        }
    }
    return simplex({secC, secA, task.b}, x0);
}

vector_t genNext(task_t task, Matrix B, vector_t xk, int jk, set<int> Nk) {
    int n = task.A.cols;
    int m = task.A.rows;
    vector_t uk(task.A.cols); // initialized with zeroes
    auto uknk = B * task.A[{jk}];
    bool lt = true;
    // there's no Nk[i] so I have to iterate
    for (int i : Nk) {
        uk[i] = uknk(i, 0);
        if (uk[i] > 0) {
            lt = false;
        }
    }
    uk[jk] = -1;
    if (lt) {
        throw logic_error("No lower bound");
    }
    double thk = 0;
    for (int i : Nk) {
        if (uk[i] > 0 && xk[i] / uk[i] < thk) {
            thk = xk[i] / uk[i];
        }
    }
    Matrix xk1 = xk - uk * thk; // I didn't make a move constructor
    return xk1;
}

set<int> fill(Matrix& A, set<int> Np) {
    // initial fill, using leftmost vectors
    set<int> perm = generator.next();
    // no exact comparison for doubles, should be det == 0
    // std::prev_permutation is used to generate sequences, so there are no repeats
    while (abs(det(A[merge(Np, perm)])) < eps) {
        perm = generator.next();
        if (perm == set<int>{}) {
            throw logic_error("looping");
        }
    }

    return merge(Np, perm);
}

// requires det!=0 due to gauss
Matrix solve(Matrix A, Matrix b, Matrix& Ainv) {
    int m = A.rows;
    auto [Qt, R] = Matrix::QtRdecomp(A);
    b = Qt * b;
    // gaussBack
    vector_t res(m);
    for (int i = m - 1; i >= 0; i--) {
        double s = 0;
        for (int k = i + 1; k < m; k++) {
            s += R(i, k) * res[k];
        }
        double x = (b(i, 0) - s) / R(i, i);
        res[i] = x;
    }
    Ainv = Matrix::RInverse(R) * Qt;
    return res;
}

// no requirements for A other than squareness, it can always be transformed
double det(Matrix A) {
    if (A.cols != A.rows) {
        throw std::invalid_argument("determinant is not defined for non-square matrices");
    }
    auto [Qt, R] = Matrix::QtRdecomp(A);
    double det = 1;
    for (int i = 0; i < R.rows; i++) {
        det *= R(i, i);
    }
    return det;
}

vector_t enumerate(task_t task) {
    int m = task.A.rows;
    int n = task.A.cols;
    Matrix B(m, m);
    vector_t x(n);
    // just build all matrices and test solutions
    set<int> N;
    for (int i = 0; i < n; i++) {
        N.insert(i);
    }
    generator.reset(N, m);
    set<int> Nk = generator.next();
    // final combination is empty
    while (Nk != set<int>{}) {
        if (abs(det(task.A[Nk])) > eps) {
            // x_t is only using m indices, while n are needed for sizes to match
            vector_t x_t = solve(task.A[Nk], task.b, B);
            if (!(x_t >= 0)) {
                Nk = generator.next();
                continue;
            }
            vector_t x_tFull(n);

            auto it = Nk.begin();
            for (int i = 0; i < m; i++) {
                x_tFull[*it] = x_t(i, 0);
                it++;
            }
            // compiler doesn't know it's 1x1
            if (dot(task.C, x_tFull) > dot(task.C, x)) {
                copy(x_tFull.begin(), x_tFull.end(), x.begin());
            }
        }
        Nk = generator.next();
    }
    return x;
}
