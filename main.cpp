#include "matrix.h"
#include "vector.h"
#include <functional>
#include <iostream>
#include <sstream>

using namespace std;

struct task_t {
    vector_t C;
    Matrix A;
    vector_t b;
    int size;         // always size of base task, for restoration
    set<int> unbound; // these are kept to be able to restore solution
};

std::pair<task_t, set<int>> read(string filename);
vector<string> split(string line, char sep);
task_t genCanon(task_t task, set<int>& noneq, int sign);
task_t genDual(task_t task, set<int>& noneq);

int main(int argc, char* argv[]) {
    auto [task, noneq] = read("task.csv");
    auto prim = genCanon(task, noneq, -1);
    cout << prim.C.T() << endl;
    cout << prim.A << endl << prim.b.T();
    auto dual = genDual(task, noneq);
    cout << dual.C.T() << endl;
    cout << dual.A << endl << dual.b.T();
    return 0;
}

std::pair<task_t, set<int>> read(string filename) {
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
                transform(A.begin() + i * width, A.begin() + (i + 1) * width, A.begin() + i * width, [](double d) { return -1 * d; });
                b[i] *= -1;
            }
        }
    }
    // unbounded variables
    getline(in, line);
    auto words = split(line, ',');
    set<int> unboundIndices;
    transform(words.begin(), words.end(), inserter(unboundIndices, unboundIndices.end()), [](string& str) { return stoi(str) - 1; });
    // I almost forgot i store in column-major, so consecutive areas are columns, not rows
    return {{C, A.T(), b, width, unboundIndices}, noneqIndices};
}

task_t genCanon(task_t task, set<int>& noneq, int sign) {
    task_t& baseline = task;

    Matrix A(baseline.A.rows, baseline.size + baseline.unbound.size() + noneq.size());
    // due to way matrices are stored, it will copy first n columns, which is exactly what I want
    copy(baseline.A.begin(), baseline.A.end(), A.begin());
    // x_i = x'_i - v_i, x_'i >=0, v_i >=0
    vector_t C = baseline.C;
    int col = baseline.size;
    for (int i : baseline.unbound) {
        C.push_back(-C[i]);
        for (int j = 0; j < A.rows; j++) {
            A(j, col) = -A(j, i);
        }
        col++;
    }
    // a >= b => a - w = b, w >=0 if sign = -1
    // or a <= b => a + w = b, w >=0 if sign = 1
    col = baseline.size + baseline.unbound.size();
    for (int i : noneq) {
        C.push_back(0);
        A(i, col++) = sign;
    }
    return {C, A, baseline.b, baseline.size, baseline.unbound};
}

task_t genDual(task_t task, set<int>& noneq) {
    Matrix A(task.A.T());
    vector_t b = task.C;
    vector_t C = task.b;
    set<int> dualUnbound;
    set<int> dualNoneq;
    // x_i >=0 => A^T[i] y <= C[i]
    for (int i = 0; i < b.rows; i++) {
        // x_i >= 0
        if (task.unbound.find(i) == task.unbound.end()) {
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
    task_t dual = {C, A, b, C.rows, dualUnbound};
    return genCanon(dual, dualNoneq, 1);
}

vector<string> split(string line, char sep) {
    stringstream iss(line);
    vector<string> ret;
    for (string word; getline(iss, word, sep);) {
        ret.push_back(word);
    }
    return ret;
}
