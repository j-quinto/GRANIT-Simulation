clear
A = importdata('C:\Users\japio\Documents\GitHub\GRANIT-Simulation\a_n_Simulation\a_n_Simulation\B_z_AC_n1234_f130_v5.5_ppi1.csv');

%{
A = importdata('C:\Users\japio\Documents\GitHub\GRANIT-Simulation\a_n_Simulation\a_n_Simulation\B_z_15G_n13_w-0.csv');
B = importdata('C:\Users\japio\Documents\GitHub\GRANIT-Simulation\a_n_Simulation\a_n_Simulation\B_z_15G_n13_w-100.csv');
C = importdata('C:\Users\japio\Documents\GitHub\GRANIT-Simulation\a_n_Simulation\a_n_Simulation\B_z_15G_n13_w-200.csv');
D = importdata('C:\Users\japio\Documents\GitHub\GRANIT-Simulation\a_n_Simulation\a_n_Simulation\B_z_15G_n13_w-250.csv');

a = A.data;
b = B.data;
%c = C.data;
%d = D.data;
ta = a(:,1);
a1a_up = a(:,2);
a1a_down = a(:,3);
a1a = sum(a(:,2:3),2);
tb = b(:,1);
a1b_up = b(:,2);
a1b_down = b(:,3);
a1b = sum(b(:,2:3),2);

tc = c(:,1);
a1c_up = c(:,2);
a1c_down = c(:,3);
a1c = sum(c(:,2:3),2);
td = d(:,1);
a1d_up = d(:,2);
a1d_down = d(:,3);
a1d = sum(d(:,2:3),2);
%}

d = A.data;
t = d(:,1);
a1_up = d(:,2);
a1_down = d(:,3);
%{
a3_up = d(:,4);
a3_down = d(:,5);
Tot = sum(d(:,2:5),2);
a1 = a1_up + a1_down;
a3 = a3_up + a3_down;
%}
a2_up = d(:,4);
a2_down = d(:,5);
a3_up = d(:,6);
a3_down = d(:,7);
a4_up = d(:,8);
a4_down = d(:,9);
%a5_up = d(:,10);
%a5_down = d(:,11);
Tot = sum(d(:,2:9),2);
a1 = a1_up + a1_down;
a2 = a2_up + a2_down;
a3 = a3_up + a3_down;
a4 = a4_up + a4_down;
%a5 = a5_up + a5_down;

%plot(t,a1_up,t,a1_down,t,a3_up,t,a3_down); legend('a1_u_p','a1_d_o_w_n','a3_u_p','a3_d_o_w_n','Total');
%plot(t, a1, t, a2, t, a3, t, a4, t, a5, t, Tot); legend('a1','a2','a3','a4','a5','Total');
%plot(t, a1, t, a3, t, Tot); legend('a1','a3','Total');
%plot(t, a1, t, Tot); legend('a_1','Total');

L = 0.16; 
v = 4.00;
to = L/v;
%N = interp1(t, a1, to);
N = a1(end)

fa = 60.0; Na = 0.0203;
f0 = 70.0; N0 = 0.0589; 
f1 = 85.0; N1 = 0.0326; 
f2 = 100.0; N2 = 0.1730; 
f3 = 105.0; N3 = 0.010903811900000;
f4 = 113.5; N4 = 0.4907; 
f5 = 130.0; N5 = 0.0742; 
f6 = 141.5; N6 = 0.5097; 
f7 = 160.0; N7 = 0.0321; 
%{
plot(fa,Na,'x',f0,N0,'x',f1,N1,'x',f2,N2,'x',f3,N3,'x',f4,N4,'x',f5,N5,...
    'x',f6,N6,'x',f7,N7,'x',f6,N,'o');
%plot([fa,f0,f1,f2,f3,f4,f5,f6],[Na,N0,N1,N2,N3,N4,N5,N6]);
xlim([0.0,300]);
ylim([0,0.7]);
xlabel('driving frequency (Hz)')
ylabel('ground state population')
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
