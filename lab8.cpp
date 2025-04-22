#include <iostream>
#include <iomanip>
#include <cmath>
 
using namespace std;
 

double f1(double x, double y) {
    return -2 * x * y * y;
}
 

double f2(double x, double y) {
    return y - x * x + 1;
}
 

double modified_euler(double (*f)(double, double), double x0, double y0, double h, double x_target) {
    double x = x0, y = y0;
    while (x < x_target) {
        double y_predict = y + h * f(x, y);  
        y = y + (h / 2) * (f(x, y) + f(x + h, y_predict)); // Corrector step
        x += h;
    }
    return y;
}
 

double heun(double (*f)(double, double), double x0, double y0, double h, double x_target) {
    double x = x0, y = y0;
    while (x < x_target) {
        double k1 = f(x, y);
        double k2 = f(x + h, y + h * k1);
        y = y + (h / 2) * (k1 + k2);
        x += h;
    }
    return y;
}
 

double runge_kutta_2nd_order(double (*f)(double, double), double x0, double y0, double h, double x_target) {
    double x = x0, y = y0;
    while (x < x_target) {
        double k1 = h * f(x, y);
        double k2 = h * f(x + h / 2, y + k1 / 2);
         y = y + (1.0 / 6.0) * (k1 + 2 * k2);
        x += h;
    }
    return y;
}
 

double percentage_error(double exact, double approx) {
    return fabs((exact - approx) / exact) * 100.0;
}
 

double exact_solution1(double x) {
    return 1 / (1 + x * x);
}
 

double exact_solution2(double x) {
    return exp(x) - x * x - 2 * x - 1;
}
 
int main() {
    double h = 0.2;
    
    // Solve using Modified Euler and Heun’s methods for first IVP
    double y_modified_euler_02 = modified_euler(f1, 0, 1, h, 0.2);
    double y_modified_euler_04 = modified_euler(f1, 0, 1, h, 0.4);
 
    double y_heun_02 = heun(f1, 0, 1, h, 0.2);
    double y_heun_04 = heun(f1, 0, 1, h, 0.4);
 
    // Exact values for comparison
    double exact_y_02 = exact_solution1(0.2);
    double exact_y_04 = exact_solution1(0.4);
 
    // Compute percentage errors
    double error_modified_02 = percentage_error(exact_y_02, y_modified_euler_02);
    double error_modified_04 = percentage_error(exact_y_04, y_modified_euler_04);
 
    double error_heun_02 = percentage_error(exact_y_02, y_heun_02);
    double error_heun_04 = percentage_error(exact_y_04, y_heun_04);
 
    cout << fixed << setprecision(6);
    cout << "First IVP: dy/dx = -2xy^2, y(0) = 1\n";
    cout << "Modified Euler:\n";
    cout << "y(0.2) = " << y_modified_euler_02 << ", Percentage Error = " << error_modified_02 << "%\n";
    cout << "y(0.4) = " << y_modified_euler_04 << ", Percentage Error = " << error_modified_04 << "%\n\n";
 
    cout << "Heun's Method:\n";
    cout << "y(0.2) = " << y_heun_02 << ", Percentage Error = " << error_heun_02 << "%\n";
    cout << "y(0.4) = " << y_heun_04 << ", Percentage Error = " << error_heun_04 << "%\n\n";
 
    return 0;
}
