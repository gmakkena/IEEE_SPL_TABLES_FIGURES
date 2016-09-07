


clear all
syms s 
%%%%% Laplace tranform of gaussian wavelet with t0=2,sigma=1 %%%%%%%%%%
q = -5.0*10^(-13)*exp(-0.25*(s - 4.0)^2)*(36631277799.0*exp(0.25*(s - 4.0)^2) + 3.544907702*exp(0.5*(s - 4.0)^2)*erfc(0.5*s - 2.0) - 32463624699.0*s*exp(0.5*(s - 4.0)^2)*erfc(0.5*s - 2.0));


%%%%%% Taylor series of the laplace transform around 0 %%%%%%
t =  taylor(q,'order',33,'ExpansionPoint',0.0);
a= sym2poly(t);
den = fliplr(a);

%%%%%%%% The number of Taylor coefficients to avoid, computed from Lemma 1%
p =1;
%%%%%%%% constant shift parameter beta %%%%%%%%%%
beta=-1.09;

%%%% Order of approximation %%%%%%%%%%%%%
k1=6;

%%%%%% Computation of denominator coefficients that are independent of beta in yhat transformation  %%%%%%%%%%
for j = 0:k1+p,
    dencoeff(1,j+1) =(-1)^(k1+p-j)* (factorial(k1+p)/(factorial(j)*factorial(k1+p-j)))/((den(1,(k1-j)+1+p)));
end

%%%%%% Computation of numerator coefficients that are independent of beta in yhat transformation  %%%%%%%%%%
for j = 0:k1+p,
   
    for i = 0:j,
        numcoeff(j+1,i+1) = (-1)^(k1+p-i)*(factorial(k1+p)/(factorial(i)*factorial(k1+p-i)))*((den(j-i+1)))/((den(k1+p-i+1)));
    end
end


%%%%%%%%%%%%% Computation of numerator and denominator coefficients dependent on
%%%%%%%%%%%%%%%% beta of the yhat transformation  %%%%%%%%% 
Den = zeros(1,k1+1);
for i=0:k1
  Den(i+1) = dencoeff(1,i+1)*pochhammer((1+k1-i+beta+p),(k1-2));
end

Num = zeros(k1+1,k1+1);

for j=0:k1
    for i = 0:j
  Num(j+1) = Num(j+1)+ numcoeff(j+1,i+1)*pochhammer((1+k1-i+beta+p),(k1-2));
    end
end

 mydenom = [];
for j=1:k1+1
    mydenom = [mydenom,Den(k1+1-(j-1))];
end
mynum=[0];
for j=2:1:k1+1
    mynum = [mynum,Num(k1+1-(j-1),1)];
end

 system=tf(mynum,mydenom);

[ha,t]=impulse(system,0:0.0001:20,'r');

minreal(system)

 
%%% Time domain representation of Gaussian wavelet t0=2,sigma=1 %%%%%
f1=- exp(-(t - 2.0).^2).*(2.0.*t - 4.0);
 

%%%%%%%%%% computation of MSE %%%%%%%%%%%%
 [yPSNR,yMSE,yMAXERR,yL2RAT] = measerr(f1,ha);
  disp('Mean Square Error');
 disp(yMSE);
