

%%%%%%%%%%% 9th order uhat approximation of 3rd derivative of B-spline of order 7 with sigma=0.5 %%%%%%% 


clear all
syms s

%%%%%%%%%%% Laplace transform of 3rd derivative of B-spline of order 7 with sigma=0.5 %%%%%%% 
q =(s/2)^3*((1-exp(-s/2))/(s/2))^7;

%%%%%% Taylor series of the laplace transform around 0 %%%%%%
t =  taylor(q,'order',33,'ExpansionPoint',0.0);
a= sym2poly(t);
den = fliplr(a);
%%%%%%%% constant shift parameter beta %%%%%%%%%%
beta=-2.32;
%%%%%%%% The number of Taylor coefficients to avoid, computed from Lemma 1%
p =4;

%%%%%%%%%%%% Order of approximation %%%%%%%%% 
k1=9;

%%%%%% Computation of denominator coefficients that are independent of beta in uhat transformation  %%%%%%%%%%

for j = 0:k1+p,
    dencoeff(1,j+1) =(-1)^(k1+p-j)* (factorial(k1+p)/(factorial(j)*factorial(k1+p-j)))/((den(1,(k1-j)+1+p)));
end
%%%%%% Computation of numerator coefficients that are independent of beta in uhat transformation  %%%%%%%%%%
for j = 0:k1+p,
    
    for i = 0:j,
        numcoeff(j+1,i+1) = (-1)^(k1+p-i)*(factorial(k1+p)/(factorial(i)*factorial(k1+p-i)))*((den(j-i+1)))/((den(k1-i+1+p)));
    end
end
%%%%%%%%%%%%% Computation of numerator and denominator coefficients dependent on
%%%%%%%%%%%%%%%% beta of the uhat transformation  %%%%%%%%% 
Den = zeros(1,k1+1);
for i=0:k1
  Den(i+1) = dencoeff(1,i+1)*(1+k1-i+beta+p).^(k1-2);
end
Num = zeros(k1+1,k1+1);

for j=0:k1
    for i = 0:j
  Num(j+1) = Num(j+1)+ numcoeff(j+1,i+1)*(1+k1-i+beta+p).^(k1-2);
    end
end

mydenom = [];
for j=1:k1+1
    mydenom = [mydenom,Den(k1+1-(j-1))];
end
mynum=[0];
for j=5:k1+1
    mynum = [mynum,Num(k1+1-(j-1),1)];
end
 system=tf(mynum,mydenom);
  [ha,t]=impulse(system,0:0.0001:20,'b');
minreal(system)
hold on
grid on 
%%%%%%%%%%% Time domain  equation of 3rd derivative of B-spline of order 7 with sigma=0.5 %%%%%%% 
f1=heaviside(t-1.0).*(t-1.0).^3.*5.6e1+heaviside(t-2.0).*(t-2.0).^3.*(2.8e2./3.0)+heaviside(t-3.0).*(t-3.0).^3.*(5.6e1./3.0)-heaviside(t-1.0./2.0).*(t-1.0./2.0).^3.*(5.6e1./3.0)-heaviside(t-3.0./2.0).*(t-3.0./2.0).^3.*(2.8e2./3.0)-heaviside(t-5.0./2.0).*(t-5.0./2.0).^3.*5.6e1-heaviside(t-7.0./2.0).*(t-7.0./2.0).^3.*(8.0./3.0)+t.^3.*(8.0./3.0);
 plot(t,f1);
  plot(t,ha);
  title('3rd derivative of B-spline of order 7 with sigma=0.5')
l=legend({'Actual wavelet ','9th order $\widehat{u}-transformation$'});
set(l,'Interpreter','Latex');
 