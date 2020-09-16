clear
%V = importdata('C:\Users\japio\Documents\GitHub\GRANIT-Simulation\a_n_Simulation\a_n_Simulation\trans4-1_Ia1.41_Ib3.82_AC_n1234_f250-300_v4.0_p0.0.csv');

%{
L = a(:,1)*4;
a1_avg = a(:,2);
plot(L,a1_avg)
peak = L(find(a1_avg == max(a1_avg)));
xline(0.16);
xline(peak,":");
legend("|a_1|^2","L_c_u_r_r_e_n_t = 0.16 m","L_m_a_X = " + string(peak) +" m");
xlabel("L")
ylabel("|a_1|^2")
title("Ground State Population (|a_1|^2) vs Transition Region Length (L): f = 199.2 Hz, v = 4 m/s, \phi = 0 rad");
%}
%{
vf = v(:,1);
vtp = v(:,2);
f = fit(vf,vtp,'gauss1');
peak = f.b1; 
plot(f,vf,vtp);
xline(peak, "--");
%xline(f.b2, ":");
xlabel('driving frequency (Hz)')
ylabel('ground state population')
%legend('Ia = 1.41, Ib = 3.82');
legend('Ia = 1.41, Ib = 3.82','fit','f_1 = '+string(peak)+' Hz');
title("Transition Probability: |4> to |1> vs driving freq (v = 4 m/s, \phi = 0 rad)");
%}

A = importdata('C:\Users\japio\Documents\GitHub\GRANIT-Simulation\a_n_Simulation\a_n_Simulation\trans3-1_Ia1.41_Ib3.82_AC_n1234_f189-213_v3.0_p0.0.csv');
a = A.data;
B = importdata('C:\Users\japio\Documents\GitHub\GRANIT-Simulation\a_n_Simulation\a_n_Simulation\trans3-1_Ia1.41_Ib3.82_AC_n1234_f180-220_v4.0_p0.0.csv');
b = B.data;
C = importdata('C:\Users\japio\Documents\GitHub\GRANIT-Simulation\a_n_Simulation\a_n_Simulation\trans4-1_Ia2.04_Ib5.2_AC_n1234_f260-290_v4.0_p0.0.csv');
c = C.data;

Af = a(:,1);
Atp = a(:,2);
Bf = b(:,1);
Btp = b(:,2);
Cf = c(:,1);
Ctp = c(:,2);

fa = fit(Af,Atp,'gauss1');
fb = fit(Bf,Btp,'gauss1');
fc = fit(Cf,Ctp,'gauss2');
peaka = fa.b1; 
peakb1 = fb.b1; 
peakc1 = fc.b1; 
plot(Af,Atp,'r*',Bf,Btp,'b*')
%xline(peaka, "--");xline(peakb1, "-");
xline((peaka+peakb1)/2);
xlabel('driving frequency (Hz)');
ylabel('ground state population');
legend('v = 3 m/s','v = 4 m/s','f^d^o^w^n_3_1 = '+string((peaka+peakb1)/2)+' Hz');
xlabel("driving frequency (Hz)");
ylabel("ground state population");
grid on
title("Transition Probability vs driving freq:(I_a, I_b) = (1.41A,3.82A),(\phi = 0 rad)");
xlim([189,213]);
%{
V1 = importdata('C:\Users\japio\Documents\GitHub\GRANIT-Simulation\a_n_Simulation\a_n_Simulation\transA1_Ia1.4_Ib3.5_AC_n1234_f141_p0_v1-7.csv');
V2 = importdata('C:\Users\japio\Documents\GitHub\GRANIT-Simulation\a_n_Simulation\a_n_Simulation\transA1_Ia1.41_Ib3.82_AC_n1234_f141_p0_v1-7.csv');
P1 = importdata('C:\Users\japio\Documents\GitHub\GRANIT-Simulation\a_n_Simulation\a_n_Simulation\transA1_Ia1.4_Ib3.5_AC_n1234_f141_p0-2pi_v4.csv');
P2 = importdata('C:\Users\japio\Documents\GitHub\GRANIT-Simulation\a_n_Simulation\a_n_Simulation\transA1_Ia1.41_Ib3.82_AC_n1234_f141_p0-2pi_v4.csv');

v1 = V1.data;
v2 = V2.data;
p1 = P1.data;
p2 = P2.data;

vel1 = v1(:,3);
vA1 = v1(:,4);
vel2 = v2(:,3);
vA2 = v2(:,4);
phi1 = p1(:,2);
pA1 = p1(:,4);
phi2 = p2(:,2);
pA2 = p2(:,4);

figure(1)
plot(vel1,vA1,"x",vel2,vA2,"o")
xlabel('velocity (m/s)')
ylabel('ground state population')
legend('Ia = 1.4, Ib = 3.5', 'Ia = 1.41, Ib = 3.82');
title("Transition Probability for |2> to |1> vs velocity (f = 141.5 Hz, phi = 0 rad)");;

figure(2)
plot(phi1,pA1,"x",phi2,pA2,"o")
xlabel('phase (rad)')
ylabel('ground state population')
legend('Ia = 1.4, Ib = 3.5', 'Ia = 1.41, Ib = 3.82');
title("Transition Probability for |2> to |1> vs phase (f = 141.5 Hz, v = 4 m/s)");
%}
%{
plot(t, a1); legend('a_1^2');
hold on
plot(to, No, "x"); legend('a_1^2: t_0 = L/v');
ylim([0,1]);
%}
%{
plot(t, a1, t, a2, t, Tot); legend('a1','a2','Total');
hold on
plot(to, No, "x"); legend('a_1^2: t_0 = L/v');
%}
%plot(ta, a1a, tb, a1b, tc, a1c, td, a1d); legend('w = 2910 rad/s','w-100 rad/s','w-200 rad/s','w-250 rad/s');
%plot(ta, a1a, tb, a1b); legend('n = 3 spin down: w-250','n = 3 spin up: w+250');
%plot(ta, a1a, tb, a1b, tc, a1c); legend('B_0 = 15G','B_0 = 50G','B_0 = 200G');
%xlabel('time (s)')
%ylabel('ground state population')
%}
