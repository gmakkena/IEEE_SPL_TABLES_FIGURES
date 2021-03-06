
clear all
syms s 

%%%%% Laplace tranform of gaussian wavelet with t0=2,sigma=1 %%%%%%%%%%
q = -5.0*10^(-13)*exp(-0.25*(s - 4.0)^2)*(36631277799.0*exp(0.25*(s - 4.0)^2) + 3.544907702*exp(0.5*(s - 4.0)^2)*erfc(0.5*s - 2.0) - 32463624699.0*s*exp(0.5*(s - 4.0)^2)*erfc(0.5*s - 2.0));

%%%%%% Taylor series of the laplace transform around 0 %%%%%%
t =  taylor(q,'order',33,'ExpansionPoint',0);

%%%%%%%%%% Conversion of Taylor series in symbolic form to Polynomial form
 a= sym2poly(t);
den = fliplr(a);

beta=0;

%%%%%%%%% The Value of p, selected from Lemma 1 %%%%%%%
p =1;

%%%%% Order of approximation %%%%%%%%
k1=5;

%%%%%% Computation of denominator coefficients that are common for proposed variants  %%%%%%%%%%
for j = 0:k1+p,
    dencoeff(1,j+1) =(-1)^(k1+p-j)* (factorial(k1+p)/(factorial(j)*factorial(k1+p-j)))/((den(1,(k1-j)+1+p)));
end


%%%%%% Computation of numerator  coefficients that are common for proposed variants  %%%%%%%%%%
for j = 0:k1+p,
    
    for i = 0:j,
        numcoeff(j+1,i+1) = (-1)^(k1+p-i)*(factorial(k1+p)/(factorial(i)*factorial(k1+p-i)))*((den(j-i+1)))/((den(k1-i+1+p)));
    end
end


%%%%%%%%%%%%% Computation of numerator and denominator coefficients dependent on
%%%%%%%%%%%%%%%% beta of the uhat transformation as beta is varied  %%%%%%%%% 

w=0;
for beta=-1.9:0.01:0.5
    w=w+1;
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

system=tf(mynum,mydenom);
 [ha,t]=impulse(system,0:0.0001:20,'b');

[ze,po,kg] = zpkdata(system);
poles=po{1,1};
    for m = 1:k1
         pqf(m) = -abs(poles(m))/(2*real(poles(m)));
    end
    pol(w)=max(pqf);
end

l=-1.9:0.01:0.5;
semilogy(l,pol,'r');
hold on


%%%%%%%%%%%%% Computation of numerator and denominator coefficients dependent on
%%%%%%%%%%%%%%%% beta of the that transformation as beta is varied  %%%%%%%%% 

w=0;
for beta=-1.9:0.01:0.5
    w=w+1;
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




[ze,po,kg] = zpkdata(system);
poles=po{1,1};
    for m = 1:k1
         pqf(m) = -abs(poles(m))/(2*real(poles(m)));
    end
    pol(w)=max(pqf);
end

l=-1.9:0.01:0.5;
semilogy(l,pol,'b');

%%%%%%%%%%%%% Computation of numerator and denominator coefficients dependent on
%%%%%%%%%%%%%%%% beta of the yhat transformation as beta is varied  %%%%%%%%% 

w=0;
for beta=-1.9:0.01:0.5
    w=w+1;
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
for j=1:1:k1+1
    mynum = [mynum,Num(k1+1-(j-1),1)];
end

 system=tf(mynum,mydenom);
 [ha,t]=impulse(system,0:0.0001:20,'b');


[ze,po,kg] = zpkdata(system);
poles=po{1,1};
    for m = 1:k1
         pqf(m) = -abs(poles(m))/(2*real(poles(m)));
    end
    pol(w)=max(pqf);
end
l=-1.9:0.01:0.5;
semilogy(l,pol,'g');

%%%%%%%%%%%%% Computation of numerator and denominator coefficients dependent on
%%%%%%%%%%%%%%%% beta of the tauhat transformation as beta is varied  %%%%%%%%% 

w=0;
for beta=-1.9:0.01:0.5
    w=w+1;
    Den = zeros(1,k1+1);
for i=0:k1
  Den(i+1) = dencoeff(1,i+1)*pochhammer((1+k1-i+beta+p),(k1-1));
end

Num = zeros(k1+1,k1+1);

for j=0:k1
    for i = 0:j
  Num(j+1) = Num(j+1)+ numcoeff(j+1,i+1)*pochhammer((1+k1-i+beta+p),(k1-1));
    end
end

 mydenom = [];
for j=1:k1+1
    mydenom = [mydenom,Den(k1+1-(j-1))];
end
mynum=[0];
for j=1:1:k1+1
    mynum = [mynum,Num(k1+1-(j-1),1)];
end

 system=tf(mynum,mydenom);
 [ha,t]=impulse(system,0:0.0001:20,'b');

[ze,po,kg] = zpkdata(system);
poles=po{1,1};
    for m = 1:k1
         pqf(m) = -abs(poles(m))/(2*real(poles(m)));
    end
    pol(w)=max(pqf);
end
l=-1.9:0.01:0.5;
semilogy(l,pol,'y');

title('Variation of Highest pole quality factor with \beta');
xlabel('\beta') % x-axis label
ylabel('Highest pole quality factor') % y-axis label
l= legend({'$\widehat{u}-transformation$','$\widehat{t}-transformation$','$\widehat{y}-transformation$','$\widehat{\tau}-transformation$'});
set(l,'Interpreter','Latex');




