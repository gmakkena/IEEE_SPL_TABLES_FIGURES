

%%%%%%%%%%%%%% This file generates the transfer function of Morlet wavelet
%%%%%%%%%%%%%% Morlet wavelet is nothing but shifted Gaussian  (t0=3) multiplied by cosine of angular frequency 6rad/sec

%%%%%%%%%%%% In generating the transfer function of Morlet, First the
%%%%%%%%%%%% shifted Gaussian Function is approximated using one of the
%%%%%%%%%%%% proposed variants. Then the multiplication of cosine is
%%%%%%%%%%%% handled similar to the method proposed by Sandro A.P Haddad
%%%%%%%%%%%% et. al in "Analog Complex Wavelet Filters" 2005. Here the
%%%%%%%%%%%% final transfer function is given.


numa=[0.06218, -0.457072, 10.9614, -25.6303, 258.129, 2549.07, -14124.7, 106006., -208949., 416198.];
dena = [1., 5.604, 208.031, 892.3, 15757.2, 49038.3, 536445.,1.08973*10^6 , 8.07667*10^6 , 8.14444*10^6, 4.21125*10^7];
sysm=tf(numa,dena);
t=0:0.0001:20;
f1=exp(-(t - 3).^2).*cos(6.*t);
 f11 = f1/max(abs(f1));
% plot(t,f11,'b');
hold on
[ha,t]=impulse(sysm,0:0.0001:20,'b');
ha1=ha/max(abs(ha));

hold on
grid on
plot(t,f1);
plot(t,ha);

title('Morlet wavelet ')
l=legend({'Actual wavelet ','10th order $\widehat{y}-transformation$'});
set(l,'Interpreter','Latex');


