// a_n_Simulation.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <Eigen/Dense>
#include "a_n_functions.h"

using namespace std;
using namespace Eigen;

int main()
{
    //RK4 func approx
    //initial conditions
    n << 1.0, 1.0, 3.0, 3.0;
    z << 13.7e-6, 13.7e-6, 32.4e-6, 32.4e-6;
    E = m_n * g * z;
    dt = 5.0e-7;
    TF = 1.0e-1;
    t = 0.0;
    a << 0.0, 0.0, 1.0, 0.0;
    Bo << 1e-3, 0.0, 1e-3,
          6e-1, 0.0, 6e-1,
           0.0, 0.0,  0.0;
    Bh << 0.0, 0.0, 15e-4;

    ofstream fout("B_z.csv");
    if (!fout)
    {
        cout << "\n error" << endl;
    }
    fout << scientific;
    fout.precision(8);

    fout << "%Time(s)" << "," << "a1_up" << "," << "a1_down" << "," << "a3_up" << "," << "a3_down" << "\n";
    while (t <= TF)
    {
        //std::cout << "time = " << t << '\n';
        double t_2 = t + 0.5 * dt;
        double t_3 = t + dt;
        
        k1 = da_dt(t, a, B(t, Bh, Bo));
        k2 = da_dt(t_2, a + dt * k1 / 2, B(t_2, Bh, Bo));
        k3 = da_dt(t_2, a + dt * k2 / 2, B(t_2, Bh, Bo));
        k4 = da_dt(t_3, a + dt * k3, B(t_3, Bh, Bo));
        a = a + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        A = a.abs2();
        //output
        fout << t << "," << A(0) << "," << A(1) << "," << A(2) << "," << A(3) << "\n";
        t = t + dt;
    }
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu