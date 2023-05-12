#include <iostream>
#include "potential.h"
#include "simplex.h"

using namespace std;

task_t constructSimplex(const transpTask& task);
Matrix transportFromSimplex(const vector_t& sol, pair<int,int> sizes);

int main(int argc, char* argv[]) {
    transpTask task = readTransport("task.csv");
    task.buildInit();
    task.solvePot();
    task.prettyPrint();
    auto taskSimple = constructSimplex(task);
    cout << "\nСимплекс:" << endl;
    // cout << taskSimple.C.T();
    // cout << taskSimple.A;
    // cout << taskSimple.b.T();
    cout << "Начальное приближение:" << endl;
    auto init = initBasic(taskSimple);
    auto transInit = transportFromSimplex(init, {task.a.size(),task.b.size()});
    cout << "Значение функции цели:" << task.results(transInit)<< endl;
    cout << transInit;
    cout << "Решение:" << endl;
    auto sol = simplex(taskSimple,init,{});
    auto transSol = transportFromSimplex(sol, {task.a.size(),task.b.size()});
    cout << "Значение функции цели:" << task.results(transSol)<< endl;
    cout << transSol;
    return 0;
}

task_t constructSimplex(const transpTask& task) {
    int m = task.a.size();
    int n = task.b.size();
    vector_t C(n*m);
    vector_t b(n+m-1);
    Matrix A(n+m-1,n*m);

    for (int i=0; i < n*m; i++) {
        A(i/n,i) = 1;
        if (i%n != n-1) {
            A(m+i%n,i) = 1;
        }
        C[i] = task.C(i/n,i%n);
    }
    copy(task.a.begin(),task.a.end(),b.begin());
    copy(task.b.begin(),task.b.end()-1,b.begin()+m);

    return {C,A,b};
}

Matrix transportFromSimplex(const vector_t& sol, pair<int,int> sizes) {
    Matrix X(sizes.first,sizes.second);
    for (int i=0; i < sizes.first; i++) {
        for (int j=0; j < sizes.second; j++) {
            X(i,j) = round(sol[i*sizes.second+j]);
        }
    }
    return X;
}
