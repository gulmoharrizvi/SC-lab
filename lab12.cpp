#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <bits/stdc++.h>

using namespace std;
const double PI = 3.141592653589793;
const double alpha = 1.0; 
 

double initialCondition(double x) {
    return sin(100 * PI * x) + sin(50 * PI * x);
}
 

void thomasSolver(const std::vector<double>& a, const std::vector<double>& b,
                  const std::vector<double>& c, std::vector<double>& d, std::vector<double>& x) {
    int n = d.size();
    std::vector<double> c_prime(n), d_prime(n);
 
    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];
 
    for (int i = 1; i < n; ++i) {
        double m = b[i] - a[i] * c_prime[i - 1];
        c_prime[i] = c[i] / m;
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / m;
    }
 
    x[n - 1] = d_prime[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }
}
 
int main() {
    int nx = 100;              
    int nt = 11;               
    double L = 1.0;           
    double T = 0.1;            
    double dx = 0.25;
    double dt = T / (nt - 1);
    double r = alpha * dt / (dx * dx);
 
    std::vector<std::vector<double>> U(nt, std::vector<double>(nx + 1, 0.0));
 
    
    for (int i = 1; i < nx; ++i) {
        double x = i * dx;
        U[0][i] = initialCondition(x);
    }
 
    
    std::vector<double> a(nx - 1, -r / 2);
    std::vector<double> b(nx - 1, 1 + r);
    std::vector<double> c(nx - 1, -r / 2);
    std::vector<double> rhs(nx - 1), sol(nx - 1);
 
    
    for (int n = 0; n < nt - 1; ++n) {
        for (int i = 1; i < nx; ++i) {
            rhs[i - 1] = (r / 2) * U[n][i - 1] + (1 - r) * U[n][i] + (r / 2) * U[n][i + 1];
        }
 
        thomasSolver(a, b, c, rhs, sol);
 
        for (int i = 1; i < nx; ++i) {
            U[n + 1][i] = sol[i - 1];
        }
    }

    for (int n = 0; n < nt; ++n) {
        std::cout << "t = " << n * dt << ": "<<endl;
        for (int i = 0; i <= nx; ++i) {
            std::cout << U[n][i] << " ";
        }
        std::cout << "\n";
    }
 
    return 0;
}
 