#include "simplex.h"
#include <algorithm>
#include <functional>

using namespace std;

double eps = 1e-8;
double eps0 = 1e-14;
struct {
    set<int> supSet;
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
        mask.resize(0);
        mask.resize(size, true);
        mask.resize(supSet.size(), false);
        first = true;
    }
} generator;

void genTest() {
    generator.reset({1, 2, 3, 4, 5}, 3);
    set<int> perm = generator.next();
    while (perm != set<int>{}) {
        cout << "{";
        for (int i : perm) {
            cout << " " << i;
        }
        cout << "}\n";
        perm = generator.next();
    }
}
void trim(Matrix& v) {
    for (auto& el : v) {
        if (abs(el) < eps0) {
            el = 0;
        }
    }
}

// should not alter initial objects
set<int> merge(set<int> N, set<int> L) {
    N.merge(L);
    return N;
}
vector_t genNext(task_t task, Matrix B, vector_t xk, int jk, set<int> Nk);
set<int> fill(Matrix& A, set<int> N);
double det(Matrix A);
Matrix solve(Matrix A, Matrix b, Matrix& Ai);
Matrix inverse(Matrix A);
void set_filter(const set<int>& supSet, set<int>& resSet, const function<bool(int)>& pred) {
    copy_if(supSet.begin(), supSet.end(), inserter(resSet, resSet.end()), pred);
}

vector_t simplex(task_t task, vector_t x0, set<int> Nk0) {
    int m = task.A.rows;
    int n = task.A.cols;
    set<int> N; // so I can replace most of for loops with set_filter
    for (int i = 0; i < n; i++) {
        N.insert(i);
    }
    set<int> Nkp;
    set<int> Nk;
    vector_t xk = x0;
    Matrix B = Matrix::eyes(m); // A^-1
    // 1.
    set_filter(N, Nkp, [&](int i) { return (xk[i] != 0); });
    set<int> Lkp;
    set_filter(N, Lkp, [&](int i) { return (Nkp.find(i) == Nkp.end()); });
    generator.reset(Lkp, m - Nkp.size());
    if (Nk0 != set<int>{}) {
        Nk = Nk0;
    } else {
        if (Nkp.size() < m) {
            Nk = fill(task.A, Nkp);
        } else {
            Nk = Nkp;
        }
    }
    B = inverse(task.A[Nk]);
    // 2.
    while (true) {
        // Lk = N \ Nk (not Nkp!)
        set<int> Lk;
        set_filter(Lkp, Lk, [&](int i) { return (Nk.find(i) == Nk.end()); });
        // whole vector to avoid index confusion like in 3.b
        vector_t dk = task.C.T() - task.C.T()[Nk] * B * task.A;
        // cout << "\nIter xk: " << xk.T();
        // cout << "dk: " << dk.T();
        trim(dk);
        // 3.a
        auto it = find_if(Lk.begin(), Lk.end(), [&dk](int j) { return (dk[j] < 0); });
        // dk >= 0
        if (it == Lk.end()) {
            return xk;
        }
        int jk = *it;
        // 3.b
        vector_t uknk = B * task.A[{jk}];
        trim(uknk);
        vector_t uk(n);
        // copy using indices I need
        int i = 0;
        for (int index : Nk) {
            uk[index] = uknk[i++];
        }
        uk[jk] = -1; // rest is 0 by construction
        // cout << "uk: " << uk.T();
        // 4.a
        set<int> P;
        // equivalent of P = { i in N | uk[i] > 0}
        set_filter(Nk, P, [&](int i) { return uk[i] > 0; });
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
        // to prevent -0. Limits maximum precision, but should prevent many bugs
        trim(xk);
        // cout << "theta: " << thk << ", C: " << dot(task.C, xk) << ", x_n: " << xk.T();
        // 5
        if (abs(thk) > eps0) { // thk != 0
            // Matrix F = Matrix::eyes(m);
            // // F dimensions are m,m. So it's possible to go out of bounds
            // int dist = distance(Nk.begin(), find(Nk.begin(), Nk.end(), ik));
            // double div = uk[ik];
            // for (int i = 0; i < m; i++) {
            //     if (i != dist) {
            //         F(i, dist) = -uknk[i] / div;
            //     } else {
            //         F(i, i) = 1 / div;
            //     }
            //     it++;
            // }
            // B = F * B;
            // recalculating sets
            Nk.erase(ik);
            Nk.insert(jk);
            Lk.erase(jk);
            Lk.insert(ik);

            Nkp.clear();
            set_filter(N, Nkp, [&](int i) { return xk[i] > 0; });
            Lkp.clear();
            set_filter(N, Lkp, [&](int i) { return xk[i] == 0; });
            // using F can swap indices (A[Nk] is always in ascending order, which is not guaranteed with F), which leads to incorrect answers
            B = inverse(task.A[Nk]);
            generator.reset(Lk, m - Nkp.size());
        } else {
            Nk = fill(task.A, Nkp);
            auto Ak = task.A[Nk];
            B = inverse(Ak);
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
            x0[n + i] = task.b[i];
        } else {
            x0[n + i] = -task.b[i];
            for (int j = 0; j < n + m; j++) {
                secA(i, j) *= -1;
            }
        }
    }
    set<int> Nk0;
    for (int i = 0; i < m; i++) {
        Nk0.insert(i + n);
    }
    vector_t res = simplex({secC, secA, task.b}, x0, Nk0);
    for (int i = 0; i < m; i++) {
        if (res[n + i] > 0) {
            throw logic_error("Empty solution space");
        }
    }
    vector_t xInit(n);
    copy(res.begin(), res.begin() + n, xInit.begin());
    return xInit;
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

Matrix inverse(Matrix A) {
    auto [Qt, R] = Matrix::QtRdecomp(A);
    return Matrix::RInverse(R) * Qt;
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
    bool first = true;
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
            trim(x_t);
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
            // cout << x_tFull.T();
            // compiler doesn't know it's 1x1
            if (!first) {
                if (dot(task.C, x_tFull) > dot(task.C, x)) {
                    copy(x_tFull.begin(), x_tFull.end(), x.begin());
                }
            } else {
                copy(x_tFull.begin(), x_tFull.end(), x.begin());
                first = false;
            }
        }
        Nk = generator.next();
    }
    return x;
}
