#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
 
using namespace std;
 

const double L = 100.0;   
const int Nx = 101;       
const double dx = L / (Nx - 1);  
const double dt = 0.1;    
const double T = 1.0;    
const int Nt = T / dt;    
 
void initialize(vector<double>& c) {
    for (int i = 0; i < Nx; i++) {
        double x = i * dx;
        if (x >= 25 && x <= 35) {
            c[i] = 1.0 + cos(M_PI * (x - 30) / 5);
        } else {
            c[i] = 0.0;
        }
    }
}
 

void solveFTCS(vector<double>& c, double (*u_func)(double)) {
    vector<double> c_new(Nx, 0.0);
 
 
    
    for (int n = 0; n < Nt; n++) {
        double t = n * dt;
        double u = u_func(t);  
        double alpha = u * dt / (2 * dx);  
 
        for (int i = 1; i < Nx - 1; i++) {
            c_new[i] = c[i] - alpha * (c[i + 1] - c[i - 1]);
        }
 
        c_new[0] = c_new[Nx - 1] = 0.0;
 
        
        c=c_new;
        
        if (n % 1 == 0 && t==0) {
            cout << "Time: " << t << " min\n";
            for (int i = 0; i < Nx; i++) {
                if(c[i]!=0)
                cout<<i<<" ";
                // cout << c[i] << "   ";
            }
            cout << "\n\n";
        }
    }
 
   
}
 

double constantSpeed(double t) {
    return 2.0;  
}
 

double timeDependentSpeed(double t) {
    return t / 20.0; 
}
 
int main() {
    vector<double> c(Nx, 0.0);
 
    initialize(c);
    solveFTCS(c, constantSpeed);
 
    cout<<"Time dependent"<<endl;
    initialize(c);
    solveFTCS(c, timeDependentSpeed);
 
    return 0;
}
 