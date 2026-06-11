#include<iostream>
#include<cmath>

const double DerivativeEpsilon=1e-5;
const double pi=std::acos(-1);

double wrap(double x) {
    double y = std::fmod(x + pi, 2.0 * pi);
    if (y <= 0.0) {
        y += 2.0 * pi;
    }
    return y - pi;
}

double Derivative(double t,double (*fc)(double),int k=1){
    const double h = DerivativeEpsilon;
    if(k==0)
        return fc(t);
    else if(k==1)
        return (fc(t+h)-fc(t-h))/(2*h);
    else if(k==2)
        return (fc(t + h) - 2.0 * fc(t) + fc(t - h)) / (h * h);
    else if(k>2){
        return (Derivative(t+h,fc,k-1)-Derivative(t-h,fc,k-1))/(2*h);
    }else
        throw std::runtime_error("invalid k");
}