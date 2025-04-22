#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
 
using namespace std;
 
// Problem parameters
const double L = 1;  // Length of the rod
const double T = 0.05;  // Total simulation time
const int Nx = 10;     // Number of spatial steps
const int Nt = 10;     // Number of time steps
 
double alpha = 1.0;    // Thermal diffusivity
double dx = L / Nx;
double dt = T / Nt;
double r = alpha * dt / (dx * dx); // Stability parameter
 
// Initial condition function
double initial_condition(double x) {
    return (x <= 0.5) ? (2 * x) : (2 * (1 - x));
}
 
// FTCS Scheme
void solve_ftcs(vector<double>& U) {
    vector<double> U_new(Nx + 1);
    
    for (int n = 0; n <= Nt; ++n) {
        cout << "Time step " << n * dt << ": ";
        for (double u : U) cout << setw(8) << fixed << setprecision(4) << u << " ";
        cout << endl;
 
        // Time stepping
        for (int i = 1; i < Nx; ++i) {
            U_new[i] = U[i] + r * (U[i - 1] - 2 * U[i] + U[i + 1]);
        }
 
        U = U_new;
        U[0] = U[Nx] = 0.0;  // Boundary conditions
    }
}
 
// Crank-Nicolson Scheme (solving A*x = b)
void solve_crank_nicolson(vector<double>& U) {
    vector<double> a(Nx - 1, -r / 2), b(Nx - 1, 1 + r), c(Nx - 1, -r / 2), d(Nx - 1), U_new(Nx + 1);
 
    for (int n = 0; n <= Nt; ++n) {
        cout << "Time step " << n * dt << ": ";
        for (double u : U) cout << setw(8) << fixed << setprecision(4) << u << " ";
        cout << endl;
 
        // Construct RHS vector d
        for (int i = 1; i < Nx; ++i) {
            d[i - 1] = (1 - r) * U[i] + (r / 2) * (U[i - 1] + U[i + 1]);
        }
 
        // Solve tridiagonal system using Thomas algorithm
        vector<double> c_star(Nx - 2), d_star(Nx - 1);
        c_star[0] = c[0] / b[0];
        d_star[0] = d[0] / b[0];
 
        for (int i = 1; i < Nx - 1; ++i) {
            double denom = b[i] - a[i] * c_star[i - 1];
            c_star[i] = c[i] / denom;
            d_star[i] = (d[i] - a[i] * d_star[i - 1]) / denom;
        }
 
        // Back substitution
        U_new[Nx - 1] = d_star[Nx - 2];
        for (int i = Nx - 3; i >= 0; --i) {
            U_new[i + 1] = d_star[i] - c_star[i] * U_new[i + 2];
        }
 
        U = U_new;
        U[0] = U[Nx] = 0.0;  // Boundary conditions
    }
}
 
int main() {
    vector<double> U(Nx + 1);
    
    // Initialize U(x,0)
    for (int i = 0; i <= Nx; ++i) {
        U[i] = initial_condition(i * dx);
    }
 
    cout << "FTCS Method:\n";
    solve_ftcs(U);
 
    // Reinitialize U for Crank-Nicolson
    for (int i = 0; i <= Nx; ++i) {
        U[i] = initial_condition(i * dx);
    }
 
    cout << "\nCrank-Nicolson Method:\n";
    solve_crank_nicolson(U);
 
    return 0;
}
 