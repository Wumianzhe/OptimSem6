#ifndef GRADIENT_H_
#define GRADIENT_H_

#include "linmethods.h"
#include "vector.h"

class Grad1 {
  public:
    Grad1(LinMethod<vector_t>& _linMin) : linMin(_linMin){};
    vector_t solve(std::function<double(vector_t)> f, std::function<vector_t(vector_t)> grad, const vector_t& x_0, const double eps, std::string of = "");
    void write(std::pair<double, double> eps_bounds, std::function<double (vector_t)> f, std::function<vector_t (vector_t)> grad, const vector_t& x_0, const vector_t& x);

  private:
    LinMethod<vector_t>& linMin;
};

#endif // GRADIENT_H_
