#include "simplex.h"
#include <algorithm>
#include <functional>

using namespace std;

double eps = 1e-8;
double eps0 = 1e-12;
struct {
    set<int> supSet;
    vector<bool> mask;
    bool first;
    set<int> next() {
        // https://stackoverflow.com/a/28698654
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

// округление до eps0. В основном для избавления от -0 и очень маленьких положительных значений
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
set<int> fill(Matrix& A, set<int> N);
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
    // Если задано начальное множество индексов (когда начинаем с решения метода исск. базиса)
    // используем его, иначе fill заполняет его так, чтобы определитель был не равен 0
    if (!Nk0.empty()) {
        Nk = Nk0;
    } else {
        if (Nkp.size() < m) {
            Nk = fill(task.A, Nkp);
        } else {
            Nk = Nkp;
        }
    }
    B = inverse(task.A[Nk]);
    while (true) {
        // разность множеств
        // Lk = N \ Nk (not Nkp!)
        set<int> Lk;
        set_filter(Lkp, Lk, [&](int i) { return (Nk.find(i) == Nk.end()); });
        // 2.
        // тут используем весь вектор dk, чтобы не преобразовывать индексы
        vector_t dk = task.C.T() - task.C.T()[Nk] * B * task.A;
        trim(dk);

        // 3.a
        auto it = find_if(Lk.begin(), Lk.end(), [&dk](int j) { return (dk[j] < 0); });
        // dk >= 0 (не найден элемент, меньший 0)
        if (it == Lk.end()) {
            return xk;
        }
        int jk = *it;

        // 3.b
        vector_t uknk = B * task.A[{jk}];
        trim(uknk);
        vector_t uk(n);
        // копируем нужные индексы
        int i = 0;
        for (int index : Nk) {
            uk[index] = uknk[i++];
        }
        uk[jk] = -1;
        // всё остальное в uk 0 по построению
         
        // 4.a
        set<int> P;
        // P = { i in N | uk[i] > 0}
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
        trim(xk);
        // cout << "theta: " << thk << ", C: " << dot(task.C, xk) << ", x_n: " << xk.T();
        
        // 5
        if (abs(thk) > eps0) { // thk != 0
            // recalculating sets
            Nk.erase(ik);
            Nk.insert(jk);
            Lk.erase(jk);
            Lk.insert(ik);

            Nkp.clear();
            set_filter(N, Nkp, [&](int i) { return xk[i] > 0; });
            Lkp.clear();
            set_filter(N, Lkp, [&](int i) { return xk[i] == 0; });

            // все красивые записи в теории о перестройке матрицы игнорируем. Я
            // использую множества, а они в c++ отсортированы. Из-за чего в
            // матрице A столбцы всегда по возрастанию индекса, что не всегда
            // так для B=A^-1, полученной построением. Как результат, получаем
            // значения не там, где они должны быть, и пару ночей дебага
            B = inverse(task.A[Nk]);
            generator.reset(Lk, m - Nkp.size());
        } else {
            Nk = fill(task.A, Nkp);
            auto Ak = task.A[Nk];
            B = inverse(Ak);
        }
    }
}

// Поиск начального опорного вектора
vector_t initBasic(task_t task) {
    int n = task.A.cols;
    int m = task.A.rows;
    // очень много копирования. В общем всё по теории
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
            task.b[i] *= -1;
            x0[n + i] = task.b[i];
            for (int j = 0; j < n; j++) {
                secA(i, j) *= -1;
            }
        }
    }
    set<int> Nk0;
    for (int i = 0; i < m; i++) {
        Nk0.insert(i + n);
    }
    trim(x0);
    vector_t res = simplex({secC, secA, task.b}, x0, Nk0);
    for (int i = 0; i < m; i++) {
        if (res[n + i] > 0 && abs(res[n + i]) > 3e-11) {
            cout << res.T();
            throw logic_error("Empty solution space");
        }
    }
    vector_t xInit(n);
    // копирование первых n компонент
    copy(res.begin(), res.begin() + n, xInit.begin());
    return xInit;
}

set<int> fill(Matrix& A, set<int> Np) {
    // initial fill, using leftmost vectors
    set<int> perm = generator.next();
    // no exact comparison for doubles, should be det == 0
    // std::prev_permutation is used to generate sequences, so there are no repeats
    double dcur = det(A[merge(Np, perm)]);
    while (isnan(dcur) || abs(dcur) < eps) {
        perm = generator.next();
        // если генератор вернул пустое мн-во, пройдены все комбинации.
        if (perm == set<int>{}) {
            throw logic_error("looping");
        }
        dcur = det(A[merge(Np, perm)]);
    }

    return merge(Np, perm);
}

Matrix inverse(Matrix A) {
    auto [Qt, R] = Matrix::QtRdecomp(A);
    return Matrix::RInverse(R) * Qt;
}

// requires det!=0 due to gauss
Matrix solve(const Matrix& A,const Matrix& b) {
    int m = A.rows;
    auto [Qt, R] = Matrix::QtRdecomp(A);
    auto bq = Qt * b;
    // gaussBack
    vector_t res(m);
    for (int i = m - 1; i >= 0; i--) {
        double s = 0;
        for (int k = i + 1; k < m; k++) {
            s += R(i, k) * res[k];
        }
        double x = (bq(i, 0) - s) / R(i, i);
        res[i] = x;
    }
    return res;
}

// no requirements for A other than squareness, it can always be transformed
// ну да, только иногда NaN возвращает...
double det(const Matrix& A) {
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
            vector_t x_t = solve(task.A[Nk], task.b);
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
            // dot это скалярное произведение (dot product)
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

std::tuple<task_t, set<int>, set<int>> genDual(task_t task, set<int>& unbound, set<int>& noneq) {
    Matrix A(task.A.T());
    vector_t b = task.C;
    vector_t C = task.b;
    set<int> dualUnbound;
    set<int> dualNoneq;
    // x_i >=0 => A^T[i] y <= C[i]
    for (int i = 0; i < b.rows; i++) {
        // x_i >= 0
        if (unbound.find(i) == unbound.end()) {
            dualNoneq.insert(i);
        }
    }
    // A[i] x >= b[i] => y_i >= 0
    for (int i = 0; i < C.rows; i++) {
        // A[i] x = b[i]
        if (noneq.find(i) == noneq.end()) {
            dualUnbound.insert(i);
        }
    }
    return {{C, A, b}, dualUnbound, dualNoneq};
}

task_t genCanon(task_t task, set<int>& unbound, set<int>& noneq, int sign) {
    task_t& baseline = task;
    int size = baseline.A.cols;

    Matrix A(baseline.A.rows, size + unbound.size() + noneq.size());
    // due to way matrices are stored, it will copy first n columns, which is exactly what I want
    copy(baseline.A.begin(), baseline.A.end(), A.begin());
    // x_i = x'_i - v_i, x_'i >=0, v_i >=0
    vector_t C = baseline.C;
    int col = size;
    for (int i : unbound) {
        C.push_back(-C[i]);
        for (int j = 0; j < A.rows; j++) {
            A(j, col) = -A(j, i);
        }
        col++;
    }
    // a >= b => a - w = b, w >=0 if sign = -1
    // or a <= b => a + w = b, w >=0 if sign = 1
    col = size + unbound.size();
    for (int i : noneq) {
        C.push_back(0);
        A(i, col++) = sign;
    }
    return {C, A, baseline.b};
}
