


%%%%%%%%%%% 8th order approximation of Compactly supported B-Spline wavelet of order 2  %%%%%%% 

clear all
syms s

%%%%% Laplace tranform of Compactly supported B-Spline wavelet of order 2 %%%%%%%%%%
q =  0.01* ((1-exp(-s/2))/(s/2))^2/12*(1 - 6*exp(-s/2) +10*exp(-2*s/2) - 6*exp(-3*s/2) + exp(-4*s/2));

%%%%%% Taylor series of the laplace transform around 0 %%%%%%
 t =  taylor(q,'order',28,'ExpansionPoint',0.0);
 a= sym2poly(t);
den = fliplr(a);

%%%%%%%%% constant shift parameter beta %%%%%%%%%%%%%%%
beta=-7.78;

%%%%%%%% The number of Taylor coefficients to avoid, computed from Lemma 1%
p =10;

%%%%%%%%%% Order of approximation %%%%%%%%%%%%
k1=8;

%%%%%% Computation of denominator coefficients that are independent of beta in that transformation  %%%%%%%%%%
for j = 0:k1+p,
    dencoeff(1,j+1) =(-1)^(k1+p-j)* (factorial(k1+p)/(factorial(j)*factorial(k1+p-j)))/((den(1,(k1-j)+1+p)));
end

%%%%%% Computation of numerator coefficients that are independent of beta in that transformation  %%%%%%%%%%

for j = 0:k1+p,
    
    for i = 0:j,
        numcoeff(j+1,i+1) = (-1)^(k1+p-i)*(factorial(k1+p)/(factorial(i)*factorial(k1+p-i)))*((den(j-i+1)))/((den(k1-i+1+p)));
    end
end


%%%%%%%%%%%%% Computation of numerator and denominator coefficients dependent on
%%%%%%%%%%%%%%%% beta of the that transformation  %%%%%%%%% 

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


%%%%%% Since the first Taylor coefficient is zero for this function , the approximation
%%%%%% obtianed already possesses one moment. No need to make the constant
%%%%%% term of the numerator zero %%%%%%%

 system=tf(mynum,mydenom);
  
 [ha,t]=impulse(system,0:0.0001:20,'b');
ha1=ha/max(abs(ha));
minreal(system)
hold on
grid on 

%%%%%%%% Time domain representation of Comapctly supported B-Spline wavelet
%%%%%%%% of order 2 %%%%%%%%%%%%%%%
f1=t.*(1.0./3.0e2)+heaviside(t-1.0).*(t-1.0).*(2.3e1./3.0e2)+heaviside(t-2.0).*(t-2.0).*(2.3e1./3.0e2)+heaviside(t-3.0).*(t-3.0).*(1.0./3.0e2)-heaviside(t-1.0./2.0).*(t-1.0./2.0).*(2.0./7.5e1)-heaviside(t-3.0./2.0).*(t-3.0./2.0).*(8.0./7.5e1)-heaviside(t-5.0./2.0).*(t-5.0./2.0).*(2.0./7.5e1);
f11 = f1/max(abs(f1));
 f2=transpose(f11);
 plot(t,f1);
  plot(t,ha );
  title('Compactly suppoted BSW of order 2')
l=legend({'Actual wavelet ','8th order $\widehat{t}-transformation$'});
set(l,'Interpreter','Latex');