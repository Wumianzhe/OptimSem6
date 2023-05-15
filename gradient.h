#ifndef GRADIENT_H_
#define GRADIENT_H_

#include "linmethods.h"
#include "vector.h"

class Grad1 {
  public:
    Grad1(std::function<double(vector_t)> f, std::function<vector_t(vector_t)> grad, vector_t sol = {0,0}) : m_f(f),m_grad(grad),m_sol(sol){};
    vector_t solve(std::function<double(vector_t)> f, std::function<vector_t(vector_t)> grad, const vector_t& x_0,
                      const double eps, const LinMethod& method, std::string of = "");
    void write(std::pair<double, double> eps_bounds, const LinMethod& method, const vector_t& x_0);

  private:
    // LinMethod& linMin;
    std::function<double(vector_t)> m_f;
    std::function<vector_t(vector_t)> m_grad;
    const vector_t m_sol;
};

class Grad2 {
  public:
    Grad2(std::function<double(vector_t)> f, std::function<vector_t(vector_t)> grad, std::function<Matrix(vector_t)> Hessian, vector_t sol = {0,0}) : m_f(f),m_grad(grad),m_sol(sol),m_Hessian(Hessian){};
    vector_t solve(std::function<double(vector_t)> f, std::function<vector_t(vector_t)> grad, const vector_t& x_0,
                      const double eps, std::string of = "");
    void write(std::pair<double, double> eps_bounds, const vector_t& x_0);
  private:
    // LinMethod& linMin;
    std::function<double(vector_t)> m_f;
    std::function<vector_t(vector_t)> m_grad;
    std::function<Matrix(vector_t)> m_Hessian;
    const vector_t m_sol;
};

#endif // GRADIENT_H_
