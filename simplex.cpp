#include "simplex.h"
#include <algorithm>

using namespace std;

double eps = 1e-8;
inline int sign(double n) { return (n >= 0); }
// should not alter initial objects
set<int> merge(set<int> N, set<int> L) {
    N.merge(L);
    return N;
}
set<int> comb(set<int> L, int count, bool reset) {
    static std::vector<bool> bitmask;
    if (reset) {
        bitmask.resize(count, true);
        bitmask.resize(L.size(), false);
    }
    if (!reset) {
        if (!prev_permutation(bitmask.begin(), bitmask.end())) {
            return {};
        }
    }
    set<int> perm;
    auto it = L.begin();
    for (int i = 0; i < bitmask.size(); i++) {
        if (bitmask[i]) {
            perm.insert(*it);
        }
        it = next(it);
    }
    return perm;
}
vector_t genNext(task_t task, Matrix B, vector_t xk, vector_t dk, set<int> N, set<int> L);
set<int> fill(task_t t, set<int> N);
double det(Matrix A);
Matrix solve(Matrix A, Matrix b, Matrix& Ai);

vector_t simplex(task_t task, vector_t x0) {
    int m = task.A.rows;
    int n = task.A.cols;
    set<int> Nkp;
    set<int> Nk;
    vector_t xk = x0;
    Matrix B = Matrix::eyes(m); // A^-1
    bool gt = true;
    do {
        // refill each time, but should not be actually long
        Nkp.clear();
        for (int i = 0; i < n; i++) {
            if (xk[i] != 0) {
                Nkp.insert(i);
            }
        }
        if (Nkp.size() < m) {
            Nk = fill(task, Nkp);
        } else {
            Nk = Nkp;
        }
        // Lk = N \ Nk (not Nkp!)
        set<int> Lk;
        for (int i = 0; i < n; i++) {
            if (Nk.find(i) == Nk.end()) {
                Lk.insert(i);
            }
        }
        // always solving for now, shortcut is WIP
        Matrix ykT = solve(task.A[Nk], task.b, B);
        Matrix dk = task.C.T() - ykT * task.A; // using full vector so there's no index confusion later
        // dk[Lk] >= 0
        for (int i = 0; i < Lk.size(); i++) {
            if (dk(i, 0) < 0) {
                gt = false;
                break;
            }
        }
        xk = genNext(task, B, xk, dk, Nk, Lk);
    } while (!gt);
    return xk;
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

vector_t genNext(task_t task, Matrix B, vector_t xk, vector_t dk, set<int> Nk, set<int> L) {
    int n = task.A.cols;
    int m = task.A.rows;
    int jk;
    for (int j : L) {
        if (dk[j] < 0) {
            jk = j;
            break;
        }
    }
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

set<int> fill(task_t t, set<int> Np) {
    // initial fill, using leftmost vectors
    int req = t.A.rows - Np.size();
    set<int> L;
    set<int> N0;
    for (int i = 0; i < t.A.rows; i++) {
        // if not already in set
        if (Np.find(i) == Np.end()) {
            L.insert(i);
            // first rows not in set
            if (N0.size() < req) {
                N0.insert(i);
            }
        }
    }
    // no exact comparison for doubles, should be det == 0
    if (abs(det(t.A[merge(Np, N0)])) < eps) {
        return merge(Np, N0);
    }
    // shuffle around (first permutation should be same as N0, but it's not a big optimization)
    // std::prev_permutation is used to generate sequences, so there are no repeats
    set<int> perm = comb(L, req, 1);
    while (abs(det(t.A[merge(Np, perm)])) < eps) {
        perm = comb(L, req, 0);
    }

    return merge(Np, perm);
}

// requires det!=0 due to gauss
Matrix solve(Matrix A, Matrix b, Matrix& Ai) {
    int m = A.rows;
    // for some reason does not work when i multiply A and b separately. Need a lot more debugging time for that
    Matrix B(A.rows, A.cols + 1);
    copy(A.begin(), A.end(), B.begin());
    copy(b.begin(), b.end(), B.begin() + A.size());

    for (int i = 0; i < m; i++) {
        double s = 0;
        for (int k = i; k < m; k++) {
            s += A(k, i) * A(k, i);
        }

        // w construction
        vector_t W(m);
        for (int k = 0; k < i; k++) {
            W[k] = 0;
        }
        W(i, 0) = A(i, i) + sign(A(i, i)) * sqrt(s);
        for (int k = i + 1; k < m; k++) {
            W[k] = A(k, i);
        }
        double beta = 1 / (s + abs(A(i, i)) * sqrt(s));

        Matrix H = Matrix::eyes(A.rows) - W * W.T() * beta;
        B = H * B;
        // unsure about it
        Ai = H * Ai;
    }

    // gaussBack
    vector_t res(m);
    for (int i = m - 1; i >= 0; i--) {
        double s = 0;
        for (int k = i + 1; k < m; k++) {
            s += A(i, k) * res[k];
        }
        double x = (b(i, 0) - s) / (A(i, i));
        res[i] = x;
    }
    return res;
}

// no requirements for A, it can always be transformed
double det(Matrix A) {
    int m = A.rows;
    for (int i = 0; i < m; i++) {
        double s = 0;
        for (int k = i; k < m; k++) {
            s += A(k, i) * A(k, i);
        }

        // w construction
        vector_t W(m);
        for (int k = 0; k < i; k++) {
            W[k] = 0;
        }
        W(i, 0) = A(i, i) + sign(A(i, i)) * sqrt(s);
        for (int k = i + 1; k < m; k++) {
            W[k] = A(k, i);
        }
        double beta = 1 / (s + abs(A(i, i)) * sqrt(s));

        Matrix H = Matrix::eyes(A.rows) - W * W.T() * beta;
        A = H * A;
    }
    double det = 1;
    for (int i = 0; i < m; i++) {
        det *= A(i, i);
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
    set<int> Nk = comb(N, m, 1);
    // final combination is empty
    while (Nk != set<int>{}) {
        if (abs(det(task.A[Nk])) > eps) {
            // x_t is only using m indices, while n are needed for sizes to match
            Matrix x_t = solve(task.A[Nk], task.b, B);
            bool positive = true;
            for (int i = 0; i < m; i++) {
                if (x_t(i, 0) < 0) {
                    positive = false;
                    break;
                }
            }
            if (!positive) {
                Nk = comb(N, m, 0);
                continue;
            }
            vector_t x_tFull(n);
            auto it = Nk.begin();
            for (int i = 0; i < m; i++) {
                x_tFull[*it] = x_t(i, 0);
                it++;
            }
            // compiler doesn't know it's 1x1
            if ((task.C.T() * x_tFull)(0, 0) > (task.C.T() * x)(0, 0)) {
                copy(x_t.begin(), x_t.end(), x.begin());
            }
        }
        Nk = comb(N, m, 0);
    }
    return x;
}
