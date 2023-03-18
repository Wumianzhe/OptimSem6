#include "potential.h"
#include <functional>
#include <iostream>
#include <sstream>
#include <stack>

enum dir {
hor,vert
};

using coord_t = std::pair<int,int>;
using namespace std;

vector<string> split(string line, char sep);
transpTask readTransport(string filename) {
    fstream in(filename);
    int height, width;
    string line;
    // sizes
    getline(in, line);
    sscanf(line.c_str(), "%d,%d\n", &height, &width);
    Matrix C(width, height);
    vector_t a(height);
    vector_t b(width);
    vector_t p(width); // shorthand for penalties

    // most coefficients
    for (int i = 0; i < height; i++) {
        getline(in, line);
        auto words = split(line, ',');
        transform(words.begin(), words.begin() + width, C.begin() + i * width, [](string& str) { return stod(str); });
        a[i] = stod(words[width]);
    }

    // required amounts
    getline(in, line);
    vector<string> words = split(line, ',');
    transform(words.begin(), words.end(), b.begin(), [](string& str) { return stod(str); });

    // penalties
    getline(in, line);
    words = split(line, ',');
    if (words.empty()) {
        p = vector_t(width);
    } else {
        transform(words.begin(), words.end(), p.begin(), [](string& str) { return stod(str); });
    }

    // open to closed conversion
    int sum = 0;
    for (int i=0 ; i < height; i ++ ) {
        sum += a[i];
    }
    for (int j=0; j < width; j++) {
        sum -= b[j];
    }
    // more required than provided
    if (sum < 0) {
        C.resize(width,height+1);
        for (int i=0; i < width; i++) {
            C(i,height) = p[i];
        }
    }

    C = C.T();

    // more provided than required
    if (sum > 0) {
        C.resize(height, width+1);
        for (int i=0; i < height; i++) {
            C(i,width) = 0;
        }
    }
    return {C, a, b};
}

vector<string> split(string line, char sep) {
    stringstream iss(line);
    vector<string> ret;
    for (string word; getline(iss, word, sep);) {
        ret.push_back(word);
    }
    return ret;
}

void transpTask::buildInit() {
    int m = a.size();
    int n = b.size();
    auto tmpa = a;
    auto tmpb = b;
    int i=0,j=0;
    for (int k=0; k < m+n-1; k++) {
        int x_k = min(tmpa[i],tmpb[j]);
        tmpa[i] -= x_k;
        tmpb[j] -= x_k;
        X(i,j) = x_k;
        if (tmpa[i] == 0) {
            i++;
            continue;
        }
        j++;
        continue;
    }
}

void transpTask::solvePot() {

}

bool transpTask::isOptimal() {
    int m = a.size();
    int n = b.size();
    vector_t v(n);
    vector_t u(m);

    stack<pair<coord_t, dir>> evals;
    for (int j = 0; j < n; j++) {
        if (X(0, j) >= 0) {
            evals.push({{0, j}, dir::vert});
        }
    }
    // пока остались нерешённые ограничения
    while (!evals.empty()) {
        auto [coords, dir] = evals.top();
        auto [i, j] = coords;
        evals.pop();
        if (dir == dir::vert) {
            v[j] = C(i, j) - u[i];
            for (int i = 0; i < m; i++) {
                if (X(i, j) >= 0) {
                    evals.push({{i, j}, dir::hor});
                }
            }
        } else {
            u[i] = C(i,j) - v[j];
            for (int j=0; j < n; j++) {
                if (X(i,j) >= 0) {
                    evals.push({{i,j}, dir::vert});
                }
            }
        }
    }

    bool optim = true;
    for (int i=0; i < m; i++) {
        for (int j=0; j < n; j++) {
            int delta = C(i,j) - v[j] + u[i];
            if (delta < 0) {
                optim = false;
            }
            if (abs(delta) > dMaxCell.first) {
                dMaxCell.first = abs(delta);
                dMaxCell.second = {i,j};
            }
        }
    }
    return optim;
}
