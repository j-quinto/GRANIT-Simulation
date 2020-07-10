%A = importdata('C:\Users\japio\Documents\GitHub\constantB.csv');
%A = importdata('C:\Users\japio\Documents\GitHub\B_z.csv');
A = importdata('C:\Users\japio\Documents\GitHub\GRANIT-Simulation\a_n_Simulation\a_n_Simulation\B_z.csv');
d = A.data;
t = d(:,1);
a1_up = d(:,2);
a1_down = d(:,3);
a3_up = d(:,4);
a3_down = d(:,5);
a1 = a1_up + a1_down;
a3 = a3_up + a3_down;
Tot = sum(d(:,2:5)');
%plot(t,a1_up,t,a1_down,t,a3_up,t,a3_down); legend('a1_u_p','a1_d_o_w_n','a3_u_p','a3_d_o_w_n','Total');
plot(t, a1, t, a3, t, Tot); legend('a1','a3','Total');
%plot(t, a1); legend('a1');
xlabel('time (s)')
ylabel('state population')
