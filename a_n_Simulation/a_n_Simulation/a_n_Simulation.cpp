// a_n_Simulation.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "a_n_functions.h"

int main()
{
    double t;
    double dt;
    double TF;
    double x[1];

    //RK4 func approx
    double k1, k2, k3, k4;

    //initial conditions
    dt = 1.0e-2;
    TF = 1.0e2;
    x[0] = 0.01; //x_0
    t = 0.0;

    ofstream fout("sim.csv");
    if (!fout)
    {
        cout << "\n error" << endl;
    }
    fout << scientific;
    fout.precision(8);

    fout << "%Time(s)" << "," << "x" << "\n";
    while (t <= TF)
    {
        //std::cout << "time = " << t << '\n';
        double t_2 = t + 0.5 * dt;
        double t_3 = t + dt;

        k1 = dt * func(x[0]);
        k2 = dt * func(x[0] + 0.5 * k1);
        k3 = dt * func(x[0] + 0.5 * k2);
        k4 = dt * func(x[0] + k3);

        x[0] = x[0] + (k1 + 2 * k2 + 2 * k3 + k4) / 6;

        //output
        fout << t << "," << x[0] << "\n";
        t = t + dt;
    }
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
