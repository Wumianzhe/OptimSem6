#include "matrix.h"
#include "simplex.h"
#include "vector.h"
#include <functional>
#include <iostream>
#include <sstream>

using namespace std;

enum extrem { min = -1, max = 1 };

std::tuple<task_t, set<int>, set<int>> read(string filename);
vector<string> split(string line, char sep);
task_t genCanon(task_t task, set<int>& unbound, set<int>& noneq, int sign);
std::tuple<task_t, set<int>, set<int>> genDual(task_t task, set<int>& unbound, set<int>& noneq);
vector_t restore(vector_t& xCan, int size, set<int>& unbound);

int main(int argc, char* argv[]) {
    auto [task, unbound, noneq] = read("task.csv");
    auto prim = genCanon(task, unbound, noneq, extrem::min);
    cout << "Прямая задача:\n";
    cout << prim.C.T() << endl;
    cout << prim.A << endl;
    cout << "b: " << prim.b.T() << endl;
    auto init = initBasic(prim);
    auto primSol = simplex(prim, init, {});
    // cout << "Решение:";
    // cout << primSol.T() << endl;
    cout << "Решение: " << restore(primSol, task.C.size(), unbound).T();
    cout << "Значение целевой функции: " << dot(prim.C, primSol) << endl;

    auto [dual, dUnbound, dNoneq] = genDual(task, unbound, noneq);
    auto dualCan = genCanon(dual, dUnbound, dNoneq, extrem::max);
    cout << "\n\nДвойственная задача:\n";
    cout << dualCan.C.T() << endl << dualCan.A << endl;
    cout << "b: " << dualCan.b.T() << endl;
    vector_t dualCanSol = enumerate(dualCan);
    // cout << "Решение:";
    // cout << dualCanSol.T() << endl;
    vector_t dualSol = restore(dualCanSol, dual.C.size(), dUnbound);
    cout << "Решение: " << dualSol.T();
    cout << "Значение целевой функции: " << dot(dualCan.C, dualCanSol) << endl;

    set<int> dualNoneq;
    for (int i = 0; i < dualSol.size(); i++) {
        if (dualSol[i] != 0) {
            dualNoneq.insert(i);
        }
    }
    cout << "Восстановленное по решению двойственной задачи решение прямой: ";
    cout << solve(task.A.T()[dualNoneq].T(), task.b.T()[dualNoneq].T()).T();
    return 0;
}

std::tuple<task_t, set<int>, set<int>> read(string filename) {
    fstream in(filename);
    int height, width;
    string line;
    // sizes
    getline(in, line);
    sscanf(line.c_str(), "%d,%d\n", &height, &width);
    vector_t C(width);
    Matrix A(width, height);
    vector_t b(height);
    // target function
    getline(in, line);
    auto funCoefs = split(line, ',');
    transform(funCoefs.begin(), funCoefs.end(), C.begin(), [](string& str) { return stod(str); });
    // reading coefs
    set<int> noneqIndices;
    for (int i = 0; i < height; i++) {
        getline(in, line);
        auto words = split(line, ',');
        transform(words.begin(), words.begin() + width, A.begin() + i * width, [](string& str) { return stod(str); });
        b[i] = stod(words[width + 1]);
        if (words[width] != "EQ") {
            noneqIndices.insert(i);
            if (words[width] == "LT") {
                transform(A.begin() + i * width, A.begin() + (i + 1) * width, A.begin() + i * width,
                          [](double d) { return -1 * d; });
                b[i] *= -1;
            }
        }
    }
    // unbounded variables
    getline(in, line);
    auto words = split(line, ',');
    set<int> unboundIndices;
    transform(words.begin(), words.end(), inserter(unboundIndices, unboundIndices.end()),
              [](string& str) { return stoi(str) - 1; });
    // I almost forgot i store in column-major, so consecutive areas are columns, not rows
    return {{C, A.T(), b}, unboundIndices, noneqIndices};
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

vector<string> split(string line, char sep) {
    stringstream iss(line);
    vector<string> ret;
    for (string word; getline(iss, word, sep);) {
        ret.push_back(word);
    }
    return ret;
}

vector_t restore(vector_t& xCan, int size, set<int>& unbound) {
    vector_t res(size);
    copy(xCan.begin(), xCan.begin() + size, res.begin());
    for (int i : unbound) {
        int ord = distance(unbound.begin(), find(unbound.begin(), unbound.end(), i));
        if (xCan[size + ord] > 0) {
            res[i] -= xCan[size + ord];
        }
    }
    return res;
}
