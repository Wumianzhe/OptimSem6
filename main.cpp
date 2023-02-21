#include "matrix.h"
#include "simplex.h"
#include "vector.h"
#include <functional>
#include <iostream>
#include <sstream>

using namespace std;

std::tuple<task_t, set<int>, set<int>> read(string filename);
vector<string> split(string line, char sep);
task_t genCanon(task_t task, set<int>& unbound, set<int>& noneq, int sign);
task_t genDual(task_t task, set<int>& unbound, set<int>& noneq);

int main(int argc, char* argv[]) {
    auto [task, unbound, noneq] = read("task.csv");
    auto prim = genCanon(task, unbound, noneq, -1);
    cout << "Прямая задача:\n";
    cout << prim.C.T() << endl;
    cout << prim.A << endl;
    cout << "b: " << prim.b.T() << endl;
    auto dual = genDual(task, unbound, noneq);
    cout << "Двойственная задача:\n";
    cout << dual.C.T() << endl << dual.A << endl;
    cout << "b: " << dual.b.T() << endl;
    vector_t dualRoot = enumerate(dual);
    cout << "Решение:";
    cout << dualRoot.T() << endl;
    cout << "Значение целевой функции: " << dual.C.T() * dualRoot;
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

task_t genDual(task_t task, set<int>& unbound, set<int>& noneq) {
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
    return genCanon({C, A, b}, dualUnbound, dualNoneq, 1);
}

vector<string> split(string line, char sep) {
    stringstream iss(line);
    vector<string> ret;
    for (string word; getline(iss, word, sep);) {
        ret.push_back(word);
    }
    return ret;
}
