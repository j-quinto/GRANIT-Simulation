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
int s;

const char* sim_mode = "transition"; //"coefficient" or "transition profile"
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
    dt = 5e-7; //1.0e-7;
    //TF = 10.0e-2;
    t = 0.0;

    //a << 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0;
    //a << sqrt(0.02 / 2), sqrt(0.02 / 2), sqrt(0.3/2), sqrt(0.3 / 2), sqrt(0.3 / 2), sqrt(0.3 / 2), sqrt(0.3 / 2), sqrt(0.3 / 2);
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
    Bh << 2e-6, 300e-6, 2e-6; //AC

    if (strcmp(sim_mode, "coefficient") == 0) {
        v = 4.0;
        TF = L / v;
        prog_len = TF / dt;
        prog_size = prog_len / 50;
        f = 199.2; //Hz
        p = 0 * pi / 4;
        for (s = 0; s <= 1; s = s + 1) //spin states
        {
            if (s == 0) {
                cout << "\n   spin up" << endl;
                a << 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0;

                ofstream uout("coeff4-1_up_Ia2.04_Ib5.2_AC_n1234_f273.9_v4.0_ppi0_L0.5m.csv");
                if (!uout)
                {
                    cout << "\n error" << endl;
                }
                uout << scientific;
                uout.precision(8);
                uout << "%Time(s)" << "," << "a1" << "," << "a1_up" << "," << "a1_down" << "," << "a2_up" << "," << "a2_down" << "," << "a3_up" << "," << "a3_down" << "," << "a4_up" << "," << "a4_down" << "\n";
                i = 1;
                TF = L / v;
                prog_len = TF / dt;
                prog_size = prog_len / 50;
                cout << "           " << 0 << "%" << "\r" << flush;
                for (t = 0.0; t <= TF; t = t + dt)
                {
                    t_2 = t + 0.5 * dt;
                    t_3 = t + dt;

                    k1 = da_dt(t, a, B(t, v, p, f, Bh, Bo, mag_mode));
                    k2 = da_dt(t_2, a + dt * k1 / 2, B(t_2, v, p, f, Bh, Bo, mag_mode));
                    k3 = da_dt(t_2, a + dt * k2 / 2, B(t_2, v, p, f, Bh, Bo, mag_mode));
                    k4 = da_dt(t_3, a + dt * k3, B(t_3, v, p, f, Bh, Bo, mag_mode));
                    a = a + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
                    A = a.abs2();
                    //output
                    i++;
                    if (i % prog_size == 0)
                    {
                        cout << "           " << 100 * (double)i / prog_len << "%" << "\r" << flush;
                        //cout << "           " << "a1 = " << A(0) + A(1) << endl;

                    }
                    uout << t << "," << A(0) + A(1) << "," << A(0) << "," << A(1) << "," << A(2) << "," << A(3) << "," << A(4) << "," << A(5) << "," << A(6) << "," << A(7) << "\n";
                }
            }
            else if (s == 1) {
                cout << "\n   spin down" << endl;
                a << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

                ofstream dout("coeff4-1_down_Ia2.04_Ib5.2_AC_n1234_f273.9_v4.0_ppi0_L0.5m.csv");
                if (!dout)
                {
                    cout << "\n error" << endl;
                }
                dout << scientific;
                dout.precision(8);
                dout << "%Time(s)" << "," << "a1" << "," << "a1_up" << "," << "a1_down" << "," << "a2_up" << "," << "a2_down" << "," << "a3_up" << "," << "a3_down" << "," << "a4_up" << "," << "a4_down" << "\n";
                i = 1;
                TF = L / v;
                prog_len = TF / dt;
                prog_size = prog_len / 50;
                cout << "           " << 0 << "%" << "\r" << flush;
                for (t = 0.0; t <= TF; t = t + dt)
                {
                    t_2 = t + 0.5 * dt;
                    t_3 = t + dt;

                    k1 = da_dt(t, a, B(t, v, p, f, Bh, Bo, mag_mode));
                    k2 = da_dt(t_2, a + dt * k1 / 2, B(t_2, v, p, f, Bh, Bo, mag_mode));
                    k3 = da_dt(t_2, a + dt * k2 / 2, B(t_2, v, p, f, Bh, Bo, mag_mode));
                    k4 = da_dt(t_3, a + dt * k3, B(t_3, v, p, f, Bh, Bo, mag_mode));
                    a = a + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
                    A = a.abs2();
                    //output
                    i++;
                    if (i % prog_size == 0)
                    {
                        cout << "           " << 100 * (double)i / prog_len << "%" << "\r" << flush;
                        //cout << "           " << "a1 = " << A(0) + A(1) << endl;

                    }
                    dout << t << "," << A(0) + A(1) << "," << A(0) << "," << A(1) << "," << A(2) << "," << A(3) << "," << A(4) << "," << A(5) << "," << A(6) << "," << A(7) << "\n";
                }
            }
        }
               
    }
    else if (strcmp(sim_mode, "transition") == 0) 
    {
        ofstream fout("trans3-1_Ia1.41_Ib3.82_AC_n1234_f189-213_v3.0_p0.0.csv");
        if (!fout)
        {
            cout << "\n error" << endl;
        }
        fout << scientific;
        fout.precision(8);

        //fout << "driving frequency(Hz)" << "," << "phase(rad)" << "," << "v(m/s)" << "," << "a1" << "\n";
        fout << "driving frequency(Hz)" << "," << "a1" << "," << "up-a1_up" << "," << "up-a1_down" << "," << "down-a1_up" << "," << "down-a1_down" << "\n";

        df = 1;
        for (f = 189; f <= 213; f = f + df)
        {
            sum = 0;
            sum1 << 0.0, 0.0,
                0.0, 0.0;
            sum2 << 0.0, 0.0,
                0.0, 0.0;
            sum3 << 0.0, 0.0,
                0.0, 0.0;
            cout << "f = " << f << endl;
            dp = pi/3.0;
            for (p = 0; p <= 0; p = p + dp)
            {
                cout << "   " << "p = " << p << endl;
                dv = 1.5;
                for (v = 3.0; v <= 3.0; v = v + dv)
                {
                    for (s = 0; s <= 1; s = s + 1) //spin states
                    {
                        if (s == 0) {
                            cout << "   spin up" << endl;
                            a << 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0;
                        }
                        else {
                            cout << "   spin down" << endl;
                            a << 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;
                        }
                        cout << "       " << "v = " << v << endl;
                        i = 1;
                        TF = L / v;
                        prog_len = TF / dt;
                        prog_size = prog_len / 50;
                        cout << "           " << 0 << "%" << "\r" << flush;
                        for (t = 0.0; t <= TF; t = t + dt)
                        {
                            t_2 = t + 0.5 * dt;
                            t_3 = t + dt;

                            k1 = da_dt(t, a, B(t, v, p, f, Bh, Bo, mag_mode));
                            k2 = da_dt(t_2, a + dt * k1 / 2, B(t_2, v, p, f, Bh, Bo, mag_mode));
                            k3 = da_dt(t_2, a + dt * k2 / 2, B(t_2, v, p, f, Bh, Bo, mag_mode));
                            k4 = da_dt(t_3, a + dt * k3, B(t_3, v, p, f, Bh, Bo, mag_mode));
                            a = a + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
                            A = a.abs2();
                            //output
                            i++;
                            if (i % prog_size == 0)
                            {
                                cout << "           " << 100 * (double)i / prog_len << "%" << "\r" << flush;
                                //cout << "           " << "a1 = " << A(0) + A(1) << endl;

                            }
                        }
                        cout << "           " << "a1 = " << A(0) + A(1) << endl;
                        //fout << f << "," << p << "," << v << "," << A(0) + A(1) << "\n";
                        if (v == 2.5) {
                            //sum1(s) = sum1(s) + (A(0) + A(1))*0.159;
                            sum1(0,s) =  A(0);
                            sum1(1,s) =  A(1);
                        }
                        else if (v == 3.0) {
                            //sum2(s) = sum2(s) + (A(0) + A(1)) * 0.682;
                            sum2(0,s) =  A(0);
                            sum2(1,s) =  A(1);
                        }
                        else if (v == 5.5) {
                            //sum1(s) = sum1(s) + (A(0) + A(1))*0.159;
                            sum3(0,s) = A(0);
                            sum3(1,s) = A(1);
                        }
                    }
                }
            }
            sum = (sum1.sum() + sum2.sum() + sum3.sum())*0.5;
            cout << "f = " << f << ", a1 = " << sum << endl;
            fout << f << "," << sum << "," << sum2(0,0) << "," << sum2(0,1) << "," << sum2(1, 0) << "," << sum2(1, 1) << "\n";
        }
    }
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu