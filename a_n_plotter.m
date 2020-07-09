A = importdata('C:\Users\japio\Documents\GitHub\GRANIT-Simulation\a_n_Simulation\a_n_Simulation\sim.csv');
d = A.data;
t = d(:,1);
a_up = d(:,2);
a_down = d(:,3);
Sum = a_up + a_down;
plot(t,a_up,t,a_down,t,Sum);