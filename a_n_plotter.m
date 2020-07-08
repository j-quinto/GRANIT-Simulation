A = importdata('C:\Users\japio\Documents\GitHub\GRANIT-Simulation\a_n_Simulation\a_n_Simulation\sim.csv');
d = A.data;
t = d(:,1);
a = d(:,2);

plot(t,a);