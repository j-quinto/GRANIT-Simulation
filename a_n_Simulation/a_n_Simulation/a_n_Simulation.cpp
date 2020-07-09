// a_n_Simulation.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <Eigen/Dense>
#include "a_n_functions.h"

using namespace std;
using namespace Eigen;

int main()
{
    double t;
    double dt;
    double TF;

    //RK4 func approx
    ArrayXcd k1(N), k2(N), k3(N), k4(N);
    ArrayXd A(N);
    //initial conditions
    n << 1.0, 1.0, 3.0, 3.0;
    z << 13.7, 13.7, 32.4, 32.4;
    E = m_n * g * z;
    dt = 1.0e-6;
    TF = 1.0e-1;
    t = 0.0;
    a << 1.0,
        0.0;
    B << 1e-3, 0, 1e-3,
        0.6, 0, 0.6,
        0, 0, 0;

    ofstream fout("sim.csv");
    if (!fout)
    {
        cout << "\n error" << endl;
    }
    fout << scientific;
    fout.precision(8);

    fout << "%Time(s)" << "," << "a1_up" << "," << "a1_down" << "\n";
    while (t <= TF)
    {
        //std::cout << "time = " << t << '\n';
        double t_2 = t + 0.5 * dt;
        double t_3 = t + dt;
        
        k1 = da_dt(t, a, B);
        k2 = da_dt(t_2, a + dt * k1 / 2, B);
        k3 = da_dt(t_2, a + dt * k2 / 2, B);
        k4 = da_dt(t_3, a + dt * k3, B);
        a = a + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        A = a.abs2();
        //output
        fout << t << "," << A(0) << "," << A(1) << "\n";
        t = t + dt;
    }
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu