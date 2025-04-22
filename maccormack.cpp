#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

const double a = 1.0;          // Wave speed
const double L = 1.0;          // Domain length
const double dx = 0.1;         // Spatial step
const double dt = 0.05;        // Time step
const double T = 0.2;          // Final time
const int Nx = L / dx + 1;     // Number of spatial points
const int Nt = T / dt;         // Number of time steps
const double pi = acos(-1.0);

double initial_condition(double x) {
    return sin(2 * pi * x);
}

double exact_solution(double x, double t) {
    return sin(2 * pi * (x - a * t));
}

int main() {
    vector<double> x(Nx), u(Nx), u_pred(Nx), u_new(Nx);

    // Setup spatial grid and initial condition
    for (int i = 0; i < Nx; ++i) {
        x[i] = i * dx;
        u[i] = initial_condition(x[i]);
    }

    // Time loop
    for (int n = 0; n < Nt; ++n) {
        // Predictor step
        for (int i = 0; i < Nx - 1; ++i) {
            u_pred[i] = u[i] - a * dt / dx * (u[i + 1] - u[i]);
        }
        u_pred[Nx - 1] = u[Nx - 1]; // Neumann/right BC approximation

        // Corrector step
        for (int i = 1; i < Nx; ++i) {
            u_new[i] = 0.5 * (u[i] + u_pred[i] - a * dt / dx * (u_pred[i] - u_pred[i - 1]));
        }
        u_new[0] = 0.0; // Dirichlet boundary condition

        u = u_new; // Update for next step
    }

    // Print final solution
    cout << fixed << setprecision(5);
    cout << "\n x\tNumerical\tExact\t\tError\n";
    double max_error = 0.0;
    for (int i = 0; i < Nx; ++i) {
        double exact = exact_solution(x[i], T);
        double error = fabs(u[i] - exact);
        max_error = max(max_error, error);
        cout << x[i] << "\t" << u[i] << "\t\t" << exact << "\t\t" << error << '\n';
    }

    cout << "\nMax error: " << max_error << endl;
    return 0;
}

