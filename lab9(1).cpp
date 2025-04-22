#include<bits/stdc++.h>

using namespace std;

#define ld long double
 
const ld tolerance = 1e-6;
 
ld f(ld x , ld y){

	return y - x*x + 1;

}
 
ld Runge_Kutta(ld x0 , ld y0 , ld h , ld target){

	ld xn = x0 , yn = y0;

	int max_iteration = 10000;

	for(int i = 0 ; i < max_iteration && target > xn ; i++){

		ld r1 = h*(yn + 1 - xn*xn - 2*xn*( (1.0/2.0) - (sqrt(3.0)/6.0) )*h - ((1.0/2.0) - (sqrt(3.0))/6.0)*((1.0/2.0) - (sqrt(3.0))/6.0)*h*h);

		ld r2 = h*(yn + 1 - xn*xn - 2*xn*( (1.0/2.0) + (sqrt(3.0)/6.0) )*h - ((1.0/2.0) + (sqrt(3.0))/6.0)*((1.0/2.0) + (sqrt(3.0))/6.0)*h*h);

		ld D = (1.0 - (h/4.0))*(1.0 - (h/4.0)) - h*h*((1.0/4.0) - (sqrt(3.0)/6.0))*((1.0/4.0) + (sqrt(3.0)/6.0));

		ld k1 = (r1*(1 - h/4.0) + h*r2*(1.0/4.0 - sqrt(3.0)/6.0))/D;

		ld k2 = (r1*(1 - h/4.0) + h*r2*(1.0/4.0 + sqrt(3.0)/6.0))/D;

		yn = yn + (k1 + k2)/2.0;

		xn = xn + h;

	}

	return yn;

}
 
 
int main(){
 
ld E1 = 1.649 , E2 = 3.18;

ld val1 = Runge_Kutta(0.0 , 0.5 , 0.2 , 0.6);

cout<<"runge Kutta for x = 0.6 "<<val1<<endl;

cout<<"Relative Error : "<<(fabs(val1 - E1)/E1)*100.0<<"%"<<endl;

cout<<endl;
 
ld val2 = Runge_Kutta(0.0 , 0.5 , 0.2 , 1.2);

cout<<"runge Kutta for x = 1.2 "<<val2<<endl;

cout<<"Relative Error : "<<(fabs(val2 - E2)/E2)*100.0<<"%"<<endl;

cout<<endl;

return 0;

}
 