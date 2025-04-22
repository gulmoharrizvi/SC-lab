#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

const double L = 1.0, T = 0.05, dx = 0.1, dt = 0.005, alpha = 1.0;
const int Nx = L / dx + 1, Nt = T / dt;
const double r = alpha * dt / (dx * dx);

double initial(double x) {
    return sin(M_PI * x);
}

double exact(double x, double t) {
    return exp(-M_PI * M_PI * t) * sin(M_PI * x);
}

void solveDufortFrankel() {
    vector<double> u_prev(Nx), u_curr(Nx), u_next(Nx), x(Nx);
    for (int i = 0; i < Nx; ++i) {
        x[i] = i * dx;
        u_prev[i] = u_curr[i] = initial(x[i]);
    }

    // First step using FTCS
    for (int i = 1; i < Nx - 1; ++i) {
        u_curr[i] = u_prev[i] + r * (u_prev[i + 1] - 2 * u_prev[i] + u_prev[i - 1]);
    }

    // Time loop
    for (int n = 1; n < Nt; ++n) {
        for (int i = 1; i < Nx - 1; ++i) {
            u_next[i] = ((1 - 2 * r) * u_prev[i] + 2 * r * (u_curr[i - 1] + u_curr[i + 1])) / (1 + 2 * r);
        }
        u_next[0] = u_next[Nx - 1] = 0;

        u_prev = u_curr;
        u_curr = u_next;
    }

    cout << "\nDufort–Frankel Scheme Results:\n";
    cout << "x\tNumerical\tExact\t\tError\n";
    double max_error = 0;
    for (int i = 0; i < Nx; ++i) {
        double ex = exact(x[i], T);
        double err = fabs(u_curr[i] - ex);
        max_error = max(max_error, err);
        cout << fixed << setprecision(5) << x[i] << "\t" << u_curr[i] << "\t\t" << ex << "\t\t" << err << '\n';
    }
    cout << "Max error: " << max_error << endl;
}
