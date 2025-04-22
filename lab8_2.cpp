#include <bits/stdc++.h>
using namespace std;
 
long double f(long double x, long double y){
    return y - x*x + 1;
}
 
long double F(long double x){
    return x*x + 2*x + 1 - 0.5*exp(x);
}
 
long double runge_kutta_2(long double x, long double y, long double h, int steps){
    for(int i = 0; i < steps; i++){
        long double k1 = h*f(x, y);
        long double k2 = h*f(x + h, y + k1);
        y = y + (k1 + k2)/2;
        x = x + h;
    }
    return y;
}
 
int main(){
    // Runge Kutta 2nd Order Method
    cout << "Runge Kutta 2nd Order Method at x = 0.6 --> " << setprecision(8) << runge_kutta_2(0.0, 0.5, 0.2, 3) << endl;
    long double y = F(0.6);
    cout << "Relative Error at x = 0.6 --> " << abs(y - runge_kutta_2(0.0, 0.5, 0.2, 3))/abs(y) << "\n";
    cout << endl;
    cout << "Runge Kutta 2nd Order Method at x = 1.2 --> " << setprecision(8) << runge_kutta_2(0.0, 0.5, 0.2, 6) << endl;
    y = F(1.2);
    cout << "Relative Error at x = 1.2 --> " << abs(y - runge_kutta_2(0.0, 0.5, 0.2, 6))/abs(y) << "\n";

}
