



clear all

syms s
%%%%% Laplace tranform of Mexican Hat wavelet with t0=3,sigma=1 %%%%%%%%%%
q = - 2.5*10^(-15)*exp(-0.25*(s - 6.0)^2)*(2.961835298*10^11*exp(0.25*(s - 6.0)^2) + 49363921633.0*s*exp(0.25*(s - 6.0)^2) - 87495273033.0*exp(0.5*(s - 6.0)^2)*erfc(0.5*s - 3.0) + 7.089815404*s*exp(0.5*(s - 6.0)^2)*erfc(0.5*s - 3.0) - 43747636500.0*s^2*exp(0.5*(s - 6.0)^2)*erfc(0.5*s - 3.0)) - 0.0001234098041*pi^(1/2)*erfc(s/2 - 3.0)*exp((s - 6.0)^2/4);

%%%%%% Taylor series of the laplace transform around 0 %%%%%%
t =  taylor(q,'order',33,'ExpansionPoint',0.0);
a= sym2poly(t);
den = fliplr(a);

%%%%%%%% The number of Taylor coefficients to avoid, computed from Lemma 1%
p =4;

%%%%%%%% constant shift parameter beta %%%%%%%%%%
beta=-2.02;

%%%% Order of approximation %%%%%%%%%%%%%
k1=7;

%%%%%% Computation of denominator coefficients that are independent of beta in t-hat transformation  %%%%%%%%%%
for j = 0:k1+p,
    dencoeff(1,j+1) =(-1)^(k1+p-j)* (factorial(k1+p)/(factorial(j)*factorial(k1+p-j)))/((den(1,(k1-j)+1+p)));
end

%%%%%% Computation of numerator coefficients that are independent of beta in t-hat transformation  %%%%%%%%%%
for j = 0:k1+p,
    
    for i = 0:j,
        numcoeff(j+1,i+1) = (-1)^(k1+p-i)*(factorial(k1+p)/(factorial(i)*factorial(k1+p-i)))*((den(j-i+1)))/((den(k1-i+1+p)));
    end
end

%%%%%%%%%%%%% Computation of numerator and denominator coefficients dependent on
%%%%%%%%%%%%%%%% beta of the t-hat transformation  %%%%%%%%% 
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
mynum(k1+1)=0;

 system=tf(mynum,mydenom);
 minreal(system)
 [ha,t]=impulse(system,0:0.0001:20);
%%%%%% Time domain representation of Mexican Hat wavelet t0=3, sigma =1
%%%%%% %%%%
 f1=exp(-(t - 3.0).^2).*(2.0.*t - 6.0).^2 - 2.0.*exp(-(t - 3.0).^2);

 
 %%%%%%%%%% computation of MSE %%%%%%%%%%%%
 [yPSNR,yMSE,yMAXERR,yL2RAT] = measerr(f1,ha);
  disp('Mean Square Error');
 disp(yMSE);
