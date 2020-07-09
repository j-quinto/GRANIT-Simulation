#pragma once
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

//double h = 1;
//double mu_n = 1;
const double pi = 3.1415926535897;
const double h = 1.054571817e-34;
const double m_n = 1.67492749805e-27;
const double g = 9.81;
const double zo = pow((pow(h, 2) / (2 * g * pow(m_n, 2))), 1 / 3);
const double gam_n = 183e6;
const double mu_n = gam_n*h/(4*pi);
const double v = 4.91;
const double d = 1e-3;
const double w = 2 * pi * v / d;

ArrayXd n(4);
ArrayXd z(4);
ArrayXd E(4);

int N = 2; //number of states/size of a
ArrayXcd a(N);
ArrayXXd B(3, 3);
ArrayXcd dadt(N, 1);

double func(double x);
//da/dt
double func(double x)
{
    return x;
}

ArrayXcd da_dt(double t, ArrayXcd a, ArrayXXcd B)
{
    dadt(0) = (1i / h) * mu_n * (a(0) * B(0, 2) + a(1) * B(0, 0));
    dadt(1) = (-1i / h) * mu_n * (a(1) * B(0, 2) - a(0) * B(0, 0));
    return dadt;
}
