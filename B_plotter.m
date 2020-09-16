mu = 4*pi*1e-7;
f = 100; %Hz
%{
Ia = 1.4;
Ib = 3.5;
Bh = [0, 300e-6, 0];

Ia = 1.41;
Ib = 3.82;
Bh = [0, 190e-6, 0];
%}
Ia = 2*2.04;
Ib = 2*5.2;
Ia2 = 2*2.04;
Ib2 = 2*5.2;
Bh = [2e-6, 300e-6, 2e-6];
%}
L = 1e-2;
c = 1e-3;
z = 1.3e-3;
k = -mu*L/(pi^2*c^2)*exp(-2*pi*z/L)*sinh(pi*c/L)*sin(pi*c/L);
%t = 0;
%x = 0:0.00001:0.16;
t = [0:0.00001:0.04];
x = 0;
a = 2*pi*x/L;
d = 2*pi/L;
w = 2*pi*f;

BxAC = k*cos(w*t)*(Ia*((sqrt(2)-2)*cos(a)+sqrt(2)*sin(a))+Ib*((sqrt(2)+2)*sin(a)-sqrt(2)*cos(a)));
BzAC = k*cos(w*t)*(Ia*((sqrt(2)-2)*sin(a)-sqrt(2)*cos(a))-Ib*((sqrt(2)+2)*cos(a)+sqrt(2)*sin(a)));
B = sqrt((BxAC+Bh(1)).^2+Bh(2)^2+(BzAC+Bh(3)).^2);
%{
BxAC2 = k*cos(w*t)*(Ia2*((sqrt(2)-2)*cos(a)+sqrt(2)*sin(a))+Ib2*((sqrt(2)+2)*sin(a)-sqrt(2)*cos(a)));
BzAC2 = k*cos(w*t)*(Ia2*((sqrt(2)-2)*sin(a)-sqrt(2)*cos(a))-Ib2*((sqrt(2)+2)*cos(a)+sqrt(2)*sin(a)));
B2 = sqrt((BxAC2+Bh(1)).^2+Bh(2)^2+(BzAC2+Bh(3)).^2);
%}
dzB = d*((BxAC + Bh(1)).*BxAC + (BzAC+Bh(3)).*BzAC)./B;
%plot(t,abs(BxAC),t,abs(BzAC),t,abs(BxAC2),t,abs(BzAC2));
%legend("B_x","B_z","B_x'","B_z'")
%plot(t,B,t,B2,t,B2./B*1e-3);
%plot(x,B,x,B2,x,B2./B*1e-3);
legend("B_I","B_2_I","(B_2_I/B_I)*10^-^3")
xlabel("time (s)");
ylabel("Field Strength (T)")
title("Magnetic Fields for Feb Currents (I) and 2x Feb Currents (2I)");
%plot(x,dzB);
%plot(t,dzB);

dzB_fit = fit(t',dzB',"fourier2");
Beta_0 = dzB_fit.a0;

%resonance frequency
m = 1.67492749804e-27; %kg neutron mass
g = 9.81; %m/s^2
hbar = 1.054571817e-34; %J*s
gam_n = 183e6;
mu_n = -1.0*gam_n*hbar/2;
z_0 = (hbar^2/(2*m^2*g))^(1/3);
e = [2.338, 4.088, 5.521, 6.787];

f_0 = m*g*z_0/(2*pi*hbar);
f_21 = f_0*(e(2)-e(1))/2;
f_21up = f_21*(1+(Beta_0*abs(mu_n))/(m*g))^(2/3);
f_21down = f_21*(1-(Beta_0*abs(mu_n))/(m*g))^(2/3);
f_31 = f_0*(e(3)-e(1))/2;
f_31up = f_31*(1+(Beta_0*abs(mu_n))/(m*g))^(2/3);
f_31down = f_31*(1-(Beta_0*abs(mu_n))/(m*g))^(2/3);
f_41 = f_0*(e(4)-e(1))/2;
f_41up = f_41*(1+(Beta_0*abs(mu_n))/(m*g))^(2/3);
f_41down = f_41*(1-(Beta_0*abs(mu_n))/(m*g))^(2/3);