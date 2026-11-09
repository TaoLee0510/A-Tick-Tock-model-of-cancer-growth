#pragma once

#include <vector>

namespace atcg3d::structured_pde {
double beta_cdf(double x, double alpha, double beta);
std::vector<double> beta_duration_kernel(double alpha, double beta, double cycle,
                                         double width, double maximum);
} // namespace atcg3d::structured_pde
