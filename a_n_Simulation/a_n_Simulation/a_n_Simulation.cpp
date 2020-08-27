// a_n_Simulation.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <Eigen/Dense>
#include "a_n_functions.h"

using namespace std;
using namespace Eigen;
int prog_len; 
int prog_size;
int i;

int main()
{
    //RK4 func approx
    //initial conditions
    //n << 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0, 5.0;
    //z << 13.7e-6, 13.7e-6, 24.0e-6, 24.0e-6, 32.4e-6, 32.4e-6, 39.8e-6, 39.8e-6, 46.6e-6, 46.6e-6;
    n << 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0;
    z << 13.7e-6, 13.7e-6, 24.0e-6, 24.0e-6, 32.4e-6, 32.4e-6, 39.8e-6, 39.8e-6;
    //n << 1.0, 1.0, 3.0, 3.0;
    //z << 13.7e-6, 13.7e-6, 32.4e-6, 32.4e-6;
    E = m_n * g * z;
    //Bh = 15e-4 dt = 5.0e-7, Bh = 50e-4 dt = 1.0e-7, Bh = 200e-4 dt - 3e-8
    dt = 1.0e-7;
    //TF = 10.0e-2;
    TF = L / v;
    t = 0.0;
    prog_len = TF / dt;
    prog_size = prog_len / 50;

    //a << 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0;
    //a << sqrt(0.02 / 2), sqrt(0.02 / 2), sqrt(0.3/2), sqrt(0.3 / 2), sqrt(0.3 / 2), sqrt(0.3 / 2), sqrt(0.3 / 2), sqrt(0.3 / 2);
    a << 0.0, 0.0, 1.0 / sqrt(2.0), 1.0 / sqrt(2.0), 0.0, 0.0, 0.0, 0.0;
    //a << 1.0 / sqrt(2.0), 1.0 / sqrt(2.0), 0.0, 0.0;
    /*
    Bo << 1e-3, 0.0, 1e-3,
          6e-1, 0.0, 6e-1,
         400.0, 0.0, 400.0;
    */
    //Bh << 0.0, 0.0, 15e-4; //DC
    
    Bo << 8.1306e-4, 0.0, -8.1306e-4,
          5.109e-1, 0.0, -5.109e-1,
          315.8273, 0.0, -315.8273;
    /*
    Bo <<  12.6e-4, 0.0, -12.6e-4,
          7.931e-1, 0.0, -7.931e-1,
          498.3077, 0.0, -498.3077;
        */
    Bh << 0.0, 0.3e-3, 0.0; //AC

    ofstream fout("B_z_AC_n1234_f130_v5.5_ppi1.csv");
    if (!fout)
    {
        cout << "\n error" << endl;
    }
    fout << scientific;
    fout.precision(8);

    //fout << "%Time(s)" << "," << "a1_up" << "," << "a1_down" << "," << "a2_up" << "," << "a2_down" << "," << "a3_up" << "," << "a3_down" << "," << "a4_up" << "," << "a4_down" << "," << "a5_up" << "," << "a5_down" << "\n";
    fout << "%Time(s)" << "," << "a1_up" << "," << "a1_down" << "," << "a2_up" << "," << "a2_down" << "," << "a3_up" << "," << "a3_down" << "," << "a4_up" << "," << "a4_down" << "\n";
    //fout << "%Time(s)" << "," << "a1_up" << "," << "a1_down" << "," << "a3_up" << "," << "a3_down" << "\n";
    cout << 0 << "%" << endl;
    i = 1;
    for(t = 0.0; t <= TF; t = t + dt)
    {
        t_2 = t + 0.5 * dt;
        t_3 = t + dt;

        k1 = da_dt(t, a, B(t, Bh, Bo, mode));
        k2 = da_dt(t_2, a + dt * k1 / 2, B(t_2, Bh, Bo, mode));
        k3 = da_dt(t_2, a + dt * k2 / 2, B(t_2, Bh, Bo, mode));
        k4 = da_dt(t_3, a + dt * k3, B(t_3, Bh, Bo, mode));
        a = a + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        A = a.abs2();
        //output
        //fout << t << "," << A(0) << "," << A(1) << "," << A(2) << "," << A(3) << "," << A(4) << "," << A(5) << "," << A(6) << "," << A(7) << "," << A(8) << "," << A(9) << "\n";
        fout << t << "," << A(0) << "," << A(1) << "," << A(2) << "," << A(3) << "," << A(4) << "," << A(5) << "," << A(6) << "," << A(7) << "\n";
        //fout << t << "," << A(0) << "," << A(1) << "," << A(2) << "," << A(3) << "\n";
        i++;
        if(i % prog_size == 0)
        {
            cout << 100 * (double) i / prog_len << "%" << endl;
            cout << "a1 = " << A(0)+A(1) << endl;

        }
    }
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu