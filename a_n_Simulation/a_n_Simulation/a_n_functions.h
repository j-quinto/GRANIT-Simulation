#pragma once
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <string>

using namespace std;
using namespace Eigen;

//double h = 1;
//double mu_n = 1;
const double pi = 3.1415926535897;
const double h = 1.054571817e-34; //Js hbar
const double m_n = 1.67492749805e-27; //kg neutron mass
const double g = 9.81; //m/s^2 gravitational constant
const double zo = pow(pow(h,2)/(2.0*g*pow(m_n, 2)),1.0/3.0);
const double gam_n = 183e6;
const double mu_n = -1.0*gam_n*h/2; //neutron magnetic moment
const double v = 5.5; //m/s
const double L = 0.16; //m/s
const double d = 1e-2; //m
const double f = 130; //Hz
const double w = 2*pi*v/d; //rad/s
const double wd = 2*pi*f;
const double p = 4*pi/4;
const double Ia = 1.4; //A
const double Ib = 3.5; //A
const double mu = 4*pi*1e-7;
const double dx = 1e-3; //mm wire cross section length
const double dz = 1e-3; //mm wire cross section height
const double Dz = 1.3e-3; //mm distance between center of wire to the mirror
const double k = -(mu*d/(pow(pi,2)*dx*dz))*exp(-2*pi*Dz/d)*sinh(pi*dz/d)*sin(pi*dx/d);
//const double w = 2910.0 - 200.0; //rad/s
double dt;
double TF;
double t;
double t_2;
double t_3;
double Bx;
double Bz;
int m;
int l;

complex<double> dadt_u;
complex<double> dadt_d;
complex<double> DaDt_u;
complex<double> DaDt_d;
complex<double> HBpp;
complex<double> HBpm;
complex<double> HBmm;
complex<double> HBmp;
complex<double> w_ml;


const int N = 8; //number of states/size of a
ArrayXcd k1(N), k2(N), k3(N), k4(N);
ArrayXd A(N);
ArrayXd n(N); //state vector
ArrayXd z(N); //m state heights 
ArrayXd E(N); //J state energies
ArrayXcd a(N); //coefficient vect
ArrayXcd dadt(N,1);
ArrayXXd b(3,3);
ArrayXXd Bo(3,3); //time dependent magnetic field coefficients
ArrayXd Bh(3); //holding magnetic field

const char* mode = "AC";

ArrayXXd B(double t, ArrayXd Bh, ArrayXXd Bo, const char* mode)
{
    if (strcmp(mode,"DC") == 0) {
        b << Bh(0) + Bo(0,0)*sin(w*t), Bh(1) + Bo(0,1), Bh(2) + Bo(0,2)*cos(w*t),
                     Bo(1,0)*sin(w*t),         Bo(1,1),         Bo(1,2)*cos(w*t),
                     Bo(2,0)*sin(w*t),         Bo(2,1),         Bo(2,2)*cos(w*t);
        return b;
    }
    else if (strcmp(mode, "AC") == 0) {
        Bx = k*cos(wd*t + p)*(Ia*((sqrt(2) - 2)*cos(w*t) + sqrt(2)*sin(w*t)) + Ib*((sqrt(2) + 2)*sin(w*t) - sqrt(2)*cos(w*t)));
        Bz = k*cos(wd*t + p)*(Ia*((sqrt(2) - 2)*sin(w*t) - sqrt(2)*cos(w*t)) - Ib*((sqrt(2) + 2)*cos(w*t) + sqrt(2)*sin(w*t)));
        b <<      Bh(0) + Bx, Bh(1),       Bh(2) + Bz,
                 (2*pi/d)*Bx,   0.0,      (2*pi/d)*Bz,
            pow(2*pi/d,2)*Bx,   0.0, pow(2*pi/d,2)*Bz;
        return b;
    }
    else{
        b << 0, 0, 0,
             0, 0, 0,
             0, 0, 0;
        return b;
    }
}

ArrayXcd da_dt(double t, ArrayXcd a, ArrayXXd B)
{
    for (m = 0; m < N; m++) 
    {
        if (m % 2 == 0) {
            DaDt_u = 0.0;
            for (l = 0; l < N; l = l + 2) {
                if (l == m) {
                    HBpp = -(B(0, 2) + B(1, 2) * (2.0 / 3.0) * z(m) + B(2, 2) * (4.0 / 15.0) * pow(z(m), 2));
                    HBpm = -(B(0, 0) - B(0, 1) * 1.0i + B(1, 0) * (2.0 / 3.0) * z(m) + B(2, 0) * (4.0 / 15.0) * pow(z(m), 2));
                }
                else {
                    HBpp = -(B(1, 2) * 2.0 * pow(-1, n(m) - n(l) - 1) * pow(zo, 3) / pow(z(m) - z(l), 2) +
                        B(2, 2) * 12.0 * pow(-1, n(m) - n(l) - 1) * pow(zo, 6) / pow(z(m) - z(l), 4));
                    HBpm = -(B(1, 0) * 2.0 * pow(-1, n(m) - n(l) - 1) * pow(zo, 3) / pow(z(m) - z(l), 2) +
                        B(2, 0) * 12.0 * pow(-1, n(m) - n(l) - 1) * pow(zo, 6) / pow(z(m) - z(l), 4));
                }
                w_ml = (E(l) - E(m)) / h;
                dadt_u = exp(-1.0i * w_ml * t) * (a(l) * HBpp + a(l + 1) * HBpm);
                DaDt_u = DaDt_u + dadt_u;
            }
            dadt(m) = (-1i * mu_n / h) * DaDt_u;
        }
        else {
            DaDt_d = 0.0;
            for (l = 1; l < N; l = l + 2) {
                if (l == m) {
                    HBmm = (B(0, 2) + B(1, 2) * (2.0 / 3.0) * z(m) + B(2, 2) * (4.0 / 15) * pow(z(m), 2));
                    HBmp = -(B(0, 0) + B(0, 1) * 1.0i + B(1, 0) * (2.0 / 3.0) * z(m) + B(2, 0) * (4.0 / 15.0) * pow(z(m), 2));
                }
                else {
                    HBmm = (B(1, 2) * 2.0 * pow(-1, n(m) - n(l) - 1) * pow(zo, 3) / pow(z(m) - z(l), 2) +
                        B(2, 2) * 12.0 * pow(-1, n(m) - n(l) - 1) * pow(zo, 6) / pow(z(m) - z(l), 4));
                    HBmp = -(B(1, 0) * 2.0 * pow(-1, n(m) - n(l) - 1) * pow(zo, 3) / pow(z(m) - z(l), 2) +
                        B(2, 0) * 12.0 * pow(-1, n(m) - n(l) - 1) * pow(zo, 6) / pow(z(m) - z(l), 4));
                }
                w_ml = (E(l) - E(m)) / h;
                dadt_d = exp(-1.0i * w_ml * t) * (a(l) * HBmm + a(l - 1) * HBmp);
                DaDt_d = DaDt_d + dadt_d;
            }
            dadt(m) = (-1i * mu_n / h)*DaDt_d;
        }
    }
    return dadt;
}

//da_dt for constant test magnetic field
/*
ArrayXcd da_dt(double t, ArrayXcd a, ArrayXXcd B)
{
    dadt(0) = (1i / h) * mu_n * (a(0) * B(0, 2) + a(1) * B(0, 0));
    dadt(1) = (-1i / h) * mu_n * (a(1) * B(0, 2) - a(0) * B(0, 0));
    return dadt;
}
*/