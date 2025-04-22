#include<bits/stdc++.h>
using namespace std;
#define  ld long double
ld func(ld x, ld y) {
    return x * (x + y) - 2;
}
ld euler(ld h,ld x0, ld y,  ld x)
{
cout<<setprecision(10);
    ld temp = -0;
 int steps=(x-x0)/h;
    while (steps--) {
        temp = y;
        y = y + h * func(x0, y);
        x0 = x0 + h;
    }
 
   return y;
}


struct GaussLegendre {
    vector<ld> points;
    vector<ld> weights;
};


GaussLegendre getGaussLegendre(int n) {
    if (n == 2) return {{-0.5773502692, 0.5773502692}, {1.0, 1.0}};
    if (n == 3) return {{-0.7745966692, 0, 0.7745966692}, {0.5555555556, 0.8888888889, 0.5555555556}};
    if (n == 5) return {{-0.9061798459, -0.5384693101, 0, 0.5384693101, 0.9061798459},
                        {0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851}};
    return {{}, {}};
}



GaussLegendre getGaussChebyshev(int n) {
    std::vector<ld> points, weights;
    for (int i = 1; i <= n; i++) {
        ld x = cos(M_PI * (2 * i - 1) / (2 * n));  // Chebyshev nodes
        points.push_back(x);
        weights.push_back(M_PI / n);
    }
    return {points, weights};
}


ld integrand(ld x) {
    return sqrt(1 - x * x) * cos(x);
}

ld gaussQuadrature(int n, bool isChebyshev = false) {
    GaussLegendre quad = isChebyshev ? getGaussChebyshev(n) : getGaussLegendre(n);
    ld integral = 0.0;
    for (size_t i = 0; i < quad.points.size(); i++) {
        ld x = quad.points[i];
        ld w = quad.weights[i];
        ld fx = integrand(x);
        integral += w * fx;
    }
    return integral;
}
 
int main() {
    ld x0 = 0, y0 = 2, xn = 3.0;
    cout << "Euler Method Approximations for y(3.0):"<<endl;
    cout << "h = 0.3: " << euler(0.3, x0, y0, xn) << endl;
    cout << "h = 0.2: " << euler(0.2, x0, y0, xn) << endl;
    cout << "h = 0.15: " << euler(0.15, x0, y0, xn) << endl;
 
   cout<<endl;
    cout << "Gauss-Legendre Integration Results:"<<endl;;
    for (int n : {2, 3, 5})
        cout << "n = " << n << ": " << gaussQuadrature(n, false) << "\n";
 
    cout << "\nGauss-Chebyshev Integration Results:\n";
    for (int n : {2, 3, 5})
        cout << "n = " << n << ": " << gaussQuadrature(n, true) << "\n";
 
    return 0;
}
