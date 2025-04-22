#include <iostream>
#include <vector>
#include <cmath>
 
// Function f(x) = sin(ln(x)) / x^2
double f(double x) {
    return sin(log(x)) / (x * x);
}
 
// Solves the tridiagonal system Ax = b using the Thomas algorithm
std::vector<double> thomas_algorithm(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d) {
    int n = d.size();
    std::vector<double> x(n);
    
    // Forward elimination
    for (int i = 1; i < n; i++) {
        double w = a[i] / b[i - 1];
        b[i] -= w * c[i - 1];
        d[i] -= w * d[i - 1];
    }
    
    // Back substitution
    x[n - 1] = d[n - 1] / b[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
    }
    
    return x;
}
 
int main() {
    double h = 0.25;
    std::vector<double> x = {1.25, 1.5, 1.75};
    int n = x.size();
 
    // Coefficients for the tridiagonal system
    std::vector<double> a(n - 1), b(n), c(n - 1), d(n);
 
    // Constructing the linear system
    for (int i = 0; i < n; i++) {
        double xi = x[i];
        double hi = h * h;
        
        a[i] = -1 - h * (1 / xi);
        b[i] = 2 + h * h * (2 / (xi * xi));
        c[i] = -1 + h * (1 / xi);
        d[i] = hi * f(xi);
    }
 
    // Boundary conditions
    d[0] -= a[0] * 1; // y(1) = 1
    d[n - 1] -= c[n - 2] * 2; // y(2) = 2
 
    // Solve using Thomas algorithm
    std::vector<double> y = thomas_algorithm(a, b, c, d);
 
    // Output results
    for (int i = 0; i < n; i++) {
        std::cout << "y(" << x[i] << ") = " << y[i] << std::endl;
    }
 
    return 0;
}
 