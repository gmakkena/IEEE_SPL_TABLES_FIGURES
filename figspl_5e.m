

%%%%%%%% 9th order tauhat approximation of 4th derivative of Exponential-Spline of order 8, alpha = 0.25
%%%%%%%% %%%%%%%


clear all
syms s

%%%%%%%%%%% Laplace transform of 4th derivative of Exponential-spline of order 8 with alpha=0.25 %%%%%%% 
q = s^4*(((1-exp((0.25-s)))/((s-0.25))))^8;

%%%%%% Taylor series of the laplace transform around 0 %%%%%%
 t =  taylor(q,'order',16,'ExpansionPoint',0.0);
 a= sym2poly(t);
den = fliplr(a);
%%%%%%%% constant shift parameter beta %%%%%%%%%%
beta=-3.13;
%%%%%%%% The number of Taylor coefficients to avoid, computed from Lemma 1%
p =5;

%%%%%%%%%%% Order of approximation  %%%%%%%%%%%
k1=9;
%%%%%% Computation of denominator coefficients that are independent of beta in tauhat transformation  %%%%%%%%%%
for j = 0:k1+p,
    dencoeff(1,j+1) =(-1)^(k1+p-j)* (factorial(k1+p)/(factorial(j)*factorial(k1+p-j)))/((den(1,(k1-j)+1+p)));
end


%%%%%% Computation of numerator coefficients that are independent of beta in tauhat transformation  %%%%%%%%%%
for j = 0:k1+p,
    
    for i = 0:j,
        numcoeff(j+1,i+1) = (-1)^(k1+p-i)*(factorial(k1+p)/(factorial(i)*factorial(k1+p-i)))*((den(j-i+1)))/((den(k1-i+1+p)));
    end
end

%%%%%%%%%%%%% Computation of numerator and denominator coefficients dependent on
%%%%%%%%%%%%%%%% beta of the tauhat transformation  %%%%%%%%% 
Den = zeros(1,k1+1);
for i=0:k1
  Den(i+1) = dencoeff(1,i+1)*(1+k1-i+beta+p).^(k1-1);
end
Num = zeros(k1+1,k1+1);

for j=0:k1
    for i = 0:j
  Num(j+1) = Num(j+1)+ numcoeff(j+1,i+1)*(1+k1-i+beta+p).^(k1-1);
    end
end

mydenom = [];
for j=1:k1+1
    mydenom = [mydenom,Den(k1+1-(j-1))];
end
mynum=[0];
for j=2:k1+1
    mynum = [mynum,Num(k1+1-(j-1),1)];
end

 system=tf(mynum,mydenom);
  
 [ha,t]=impulse(system,0:0.0001:20,'b');
ha1=ha/max(abs(ha));
minreal(system)
hold on
grid on 
%%%%%%%%%%% Time domain equatio of  4th derivative of Exponential-spline of order 8 with alpha=0.25 %%%%%%% 
f1=t.^3.*exp(t.*(1.0./4.0)).*(1.0./6.0)+t.^4.*exp(t.*(1.0./4.0)).*(1.0./2.4e1)+t.^5.*exp(t.*(1.0./4.0)).*(1.0./3.2e2)+t.^6.*exp(t.*(1.0./4.0)).*8.680555555555556e-5+t.^7.*exp(t.*(1.0./4.0)).*7.750496031746032e-7+heaviside(t-2.0).*exp(1.0./2.0).*(exp(t.*(1.0./4.0)-1.0./2.0).*(t-2.0).^3.*(1.0./6.0)+exp(t.*(1.0./4.0)-1.0./2.0).*(t-2.0).^4.*(1.0./2.4e1)+exp(t.*(1.0./4.0)-1.0./2.0).*(t-2.0).^5.*(1.0./3.2e2)+exp(t.*(1.0./4.0)-1.0./2.0).*(t-2.0).^6.*8.680555555555556e-5+exp(t.*(1.0./4.0)-1.0./2.0).*(t-2.0).^7.*7.750496031746032e-7).*2.8e1+heaviside(t-4.0).*exp(1.0).*(exp(t.*(1.0./4.0)-1.0).*(t-4.0).^3.*(1.0./6.0)+exp(t.*(1.0./4.0)-1.0).*(t-4.0).^4.*(1.0./2.4e1)+exp(t.*(1.0./4.0)-1.0).*(t-4.0).^5.*(1.0./3.2e2)+exp(t.*(1.0./4.0)-1.0).*(t-4.0).^6.*8.680555555555556e-5+exp(t.*(1.0./4.0)-1.0).*(t-4.0).^7.*7.750496031746032e-7).*7.0e1-heaviside(t-1.0).*exp(1.0./4.0).*(exp(t.*(1.0./4.0)-1.0./4.0).*(t-1.0).^3.*(1.0./6.0)+exp(t.*(1.0./4.0)-1.0./4.0).*(t-1.0).^4.*(1.0./2.4e1)+exp(t.*(1.0./4.0)-1.0./4.0).*(t-1.0).^5.*(1.0./3.2e2)+exp(t.*(1.0./4.0)-1.0./4.0).*(t-1.0).^6.*8.680555555555556e-5+exp(t.*(1.0./4.0)-1.0./4.0).*(t-1.0).^7.*7.750496031746032e-7).*8.0-heaviside(t-3.0).*exp(3.0./4.0).*(exp(t.*(1.0./4.0)-3.0./4.0).*(t-3.0).^3.*(1.0./6.0)+exp(t.*(1.0./4.0)-3.0./4.0).*(t-3.0).^4.*(1.0./2.4e1)+exp(t.*(1.0./4.0)-3.0./4.0).*(t-3.0).^5.*(1.0./3.2e2)+exp(t.*(1.0./4.0)-3.0./4.0).*(t-3.0).^6.*8.680555555555556e-5+exp(t.*(1.0./4.0)-3.0./4.0).*(t-3.0).^7.*7.750496031746032e-7).*5.6e1+heaviside(t-8.0).*exp(2.0).*(exp(t.*(1.0./4.0)-2.0).*(t-8.0).^3.*(1.0./6.0)+exp(t.*(1.0./4.0)-2.0).*(t-8.0).^4.*(1.0./2.4e1)+exp(t.*(1.0./4.0)-2.0).*(t-8.0).^5.*(1.0./3.2e2)+exp(t.*(1.0./4.0)-2.0).*(t-8.0).^6.*8.680555555555556e-5+exp(t.*(1.0./4.0)-2.0).*(t-8.0).^7.*7.750496031746032e-7)+heaviside(t-6.0).*exp(3.0./2.0).*(exp(t.*(1.0./4.0)-3.0./2.0).*(t-6.0).^3.*(1.0./6.0)+exp(t.*(1.0./4.0)-3.0./2.0).*(t-6.0).^4.*(1.0./2.4e1)+exp(t.*(1.0./4.0)-3.0./2.0).*(t-6.0).^5.*(1.0./3.2e2)+exp(t.*(1.0./4.0)-3.0./2.0).*(t-6.0).^6.*8.680555555555556e-5+exp(t.*(1.0./4.0)-3.0./2.0).*(t-6.0).^7.*7.750496031746032e-7).*2.8e1-heaviside(t-5.0).*exp(5.0./4.0).*(exp(t.*(1.0./4.0)-5.0./4.0).*(t-5.0).^3.*(1.0./6.0)+exp(t.*(1.0./4.0)-5.0./4.0).*(t-5.0).^4.*(1.0./2.4e1)+exp(t.*(1.0./4.0)-5.0./4.0).*(t-5.0).^5.*(1.0./3.2e2)+exp(t.*(1.0./4.0)-5.0./4.0).*(t-5.0).^6.*8.680555555555556e-5+exp(t.*(1.0./4.0)-5.0./4.0).*(t-5.0).^7.*7.750496031746032e-7).*5.6e1-heaviside(t-7.0).*exp(7.0./4.0).*(exp(t.*(1.0./4.0)-7.0./4.0).*(t-7.0).^3.*(1.0./6.0)+exp(t.*(1.0./4.0)-7.0./4.0).*(t-7.0).^4.*(1.0./2.4e1)+exp(t.*(1.0./4.0)-7.0./4.0).*(t-7.0).^5.*(1.0./3.2e2)+exp(t.*(1.0./4.0)-7.0./4.0).*(t-7.0).^6.*8.680555555555556e-5+exp(t.*(1.0./4.0)-7.0./4.0).*(t-7.0).^7.*7.750496031746032e-7).*8.0;
 plot(t,f1);
  plot(t,ha);
    title('4th derivative of Exponential-Spline of order 8, alpha = 0.25')
l=legend({'Actual wavelet ','9th order $\widehat{\tau}-transformation$'});
set(l,'Interpreter','Latex');
