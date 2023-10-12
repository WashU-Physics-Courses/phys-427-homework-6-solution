#include "rk4.h"

#include <iostream>
#include <fstream>

int main() {
  // Potential function
  double k0 = 1.0e-4;
  double k1 = 1.0e-7;
  auto V = [k0, k1](double x) {
    return 0.5 * k0 * x * x + k1 * x * x * x * x;
  };

  // Stop condition. Stop when hitting the other boundary at x = a
  double a = 30.0;
  auto stop = [a](double x, const std::vector<double>& y) {
    return x >= a;
  };

  // Right hand side of the differential equation. E is captured by reference,
  // so that we can update it afterwards
  double E = 0.0;
  auto dydx = [&E, V] (double x, const std::vector<double>& y) {
    return std::vector<double>{y[1], 2.0*(V(x) - E) * y[0]};
  };

  runge_kutta_4 rk4(2);

  // Initial condition at boundary x = -a
  double dphi0 = 0.001;
  std::vector<double> y0{0.0, dphi0};

  // Constants
  double phi_a = 100.0;
  double h_a = 1e-6;
  double x0 = -a;
  double dF = 0.0;

  // Target tolerance is 10^-8
  double tol = 1e-8;
  int n = 0;

  while(std::abs(phi_a) > tol) {
    // debug message
    // std::cout << n << std::endl;

    // Try to integrate with current value of E to obtain psi(a)
    auto y = rk4.integrate(dydx, stop, 0.1, x0, y0);
    phi_a = y[0];
    // Break if we are already successful
    if (std::abs(phi_a) <= tol)
      break;

    // Increment E in order to estimate a numerical derivative
    E = E + h_a;
    auto y_delta = rk4.integrate(dydx, stop, 0.1, x0, y0);
    // decrease E again so that we are at the same E where we obtained phi_a
    E = E - h_a;
    // Numerical derivative
    dF = (y_delta[0] - phi_a) / h_a;
    // Newton iteration
    E -= phi_a / dF;

    // debug message
    // std::cout << "E is " << E << ", phi_a is " << phi_a << ", y_delta is " << y_delta[0] << std::endl;

    // Do not iterate for too many steps
    n += 1;
    if (n > 100) break;
  }
  // Print the energy eigenvalue
  std::cout << "E = " << E << std::endl;

  // Write the wave function to an output file
  std::ofstream outfile("problem1.csv");
  for (int i = 0; i < rk4.result.size(); i++) {
    outfile << rk4.xs[i] << "," << rk4.result[i][0] << std::endl;
  }
  outfile.close();

  return 0;
}
