
%%%%%%%%%%% 7th order approximation of Mexican Hat wavelet with shift 3 and
%%%%%%%%%%% scaling parameter 1 obtained by proposed uhat transformation 



clear all
syms s

%%%%% Laplace tranform of Mexican Hat wavelet with t0=3,sigma=1 %%%%%%%%%%
q = - 2.5*10^(-15)*exp(-0.25*(s - 6.0)^2)*(2.961835298*10^11*exp(0.25*(s - 6.0)^2) + 49363921633.0*s*exp(0.25*(s - 6.0)^2) - 87495273033.0*exp(0.5*(s - 6.0)^2)*erfc(0.5*s - 3.0) + 7.089815404*s*exp(0.5*(s - 6.0)^2)*erfc(0.5*s - 3.0) - 43747636500.0*s^2*exp(0.5*(s - 6.0)^2)*erfc(0.5*s - 3.0)) - 0.0001234098041*pi^(1/2)*erfc(s/2 - 3.0)*exp((s - 6.0)^2/4);

%%%%%% Taylor series of the laplace transform around 0 %%%%%%
t =  taylor(q,'order',33,'ExpansionPoint',0.0);
a= sym2poly(t);
den = fliplr(a);

%%%%%%%% constant shift parameter beta %%%%%%%%%%
beta=-1.16;

%%%%%%%% The number of Taylor coefficients to avoid, computed from Lemma 1%
p =2;

%%%% Order of approximation %%%%%%%%%%%%%
k1=7;
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
for j=2:k1+1
    mynum = [mynum,Num(k1+1-(j-1),1)];
end

%%%%%%%%%%% Making the constant term of the numerator polynomial zero to
%%%%%%%%%%% ensure the admissibility property %%%%%%%%%%%%%%
 mynum(k1+1)=0;
 
system=tf(mynum,mydenom);
[ha,t]=impulse(system,0:0.0001:20,'b');
ha1=ha/max(abs(ha));
minreal(system)
hold on
grid on 

%%%%%% Time domain representation of Mexican Hat wavelet t0=3, sigma =1
%%%%%% %%%%
f1=exp(-(t - 3.0).^2).*(2.0.*t - 6.0).^2 - 2.0.*exp(-(t - 3.0).^2);
 plot(t,f1);
  plot(t,ha);
  title('Mexican hat wavelet t0=3, sigma=1');
l=legend({'Actual wavelet ','7th order $\widehat{u}-transformation$'});
set(l,'Interpreter','Latex');