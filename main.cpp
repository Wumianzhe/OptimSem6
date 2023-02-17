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
task_t genPrimary(std::pair<task_t, set<int>>& generic);
task_t genDual(std::pair<task_t, set<int>>& generic);

int main(int argc, char* argv[]) {
    auto generic = read("task.csv");
    auto prim = genPrimary(generic);
    cout << prim.C.T() << endl;
    cout << prim.A << endl << prim.b.T();
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
                transform(A.begin() + i * width, A.begin() + (i + 1) * width + 1, A.begin() + i * width, [](double d) { return -1 * d; });
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
    A = A.T();
    task_t task = {C, A, b, width, unboundIndices};
    return {task, noneqIndices};
}

task_t genPrimary(std::pair<task_t, set<int>>& generic) {
    task_t& baseline = generic.first;
    set<int>& noneq = generic.second;

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
    // a >= b => a - w = b, w >=0
    col = baseline.size + baseline.unbound.size();
    for (int i : noneq) {
        C.push_back(0);
        A(i, col++) = -1;
    }
    return {C, A, baseline.b, baseline.size, baseline.unbound};
}

vector<string> split(string line, char sep) {
    stringstream iss(line);
    vector<string> ret;
    for (string word; getline(iss, word, sep);) {
        ret.push_back(word);
    }
    return ret;
}
