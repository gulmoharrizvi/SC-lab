#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

const double a = 1.0, L = 1.0, dx = 0.1, dt = 0.05, T = 0.2;
const int Nx = L / dx + 1, Nt = T / dt;
const double pi = acos(-1.0);

double initial(double x) {
    return sin(2 * pi * x);
}

double exact(double x, double t) {
    return sin(2 * pi * (x - a * t));
}

void solveUpwind() {
    vector<double> x(Nx), u(Nx), u_new(Nx);

    for (int i = 0; i < Nx; ++i) {
        x[i] = i * dx;
        u[i] = initial(x[i]);
    }

    for (int n = 0; n < Nt; ++n) {
        for (int i = 1; i < Nx; ++i)
            u_new[i] = u[i] - a * dt / dx * (u[i] - u[i - 1]);
        u_new[0] = 0.0;  // Dirichlet BC
        u = u_new;
    }

    cout << "\nUpwind Scheme Results:\n";
    cout << "x\tNumerical\tExact\t\tError\n";
    double max_error = 0;
    for (int i = 0; i < Nx; ++i) {
        double ex = exact(x[i], T);
        double err = fabs(u[i] - ex);
        max_error = max(max_error, err);
        cout << fixed << setprecision(5) << x[i] << "\t" << u[i] << "\t\t" << ex << "\t\t" << err << '\n';
    }
    cout << "Max error: " << max_error << endl;
}
