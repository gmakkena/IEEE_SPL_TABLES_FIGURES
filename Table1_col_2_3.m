clear all
syms s

%%%%% Laplace tranform of Gaussian wavelet with t0=2,sigma=1 %%%%%%%%%%
q = -5.0*10^(-13)*exp(-0.25*(s - 4.0)^2)*(36631277799.0*exp(0.25*(s - 4.0)^2) + 3.544907702*exp(0.5*(s - 4.0)^2)*erfc(0.5*s - 2.0) - 32463624699.0*s*exp(0.5*(s - 4.0)^2)*erfc(0.5*s - 2.0));
%%%%%% Taylor series of the laplace transform around 0 %%%%%%
t =  taylor(q,'order',33,'ExpansionPoint',0);
 a= sym2poly(t);
den = fliplr(a);


%%%%% Order of the approximation
k1=5;
%%%%%% constant shift parameter %%%%
beta=0;

%%%%%%%%%% u-transformation %%%%%%%%%%%

%%%%%% Computation of denominator coefficients that are independent of beta %%%%%%%%%%
for j = 0:k1,
    dencoeff(1,j+1) =(-1)^(k1-j)* (factorial(k1)/(factorial(j)*factorial(k1-j)))/((den(1,(k1-j)+1)));
end

%%%%%% Computation of numerator  coefficients that are common for proposed variants  %%%%%%%%%%
for j = 0:k1,
    
    for i = 0:j,
        numcoeff(j+1,i+1) = (-1)^(k1-i)*(factorial(k1)/(factorial(i)*factorial(k1-i)))*((den(j-i+1)))/((den(k1-i+1)));
    end
end

%%%%%%%%%%%%% Computation of numerator and denominator coefficients
%%%%%%%%%%%%% dependent on beta %%%%%%%%%%%%%%%%%%%%%%%%
Den = zeros(1,k1+1);
for i=0:k1
  Den(i+1) = dencoeff(1,i+1)*(1+k1-i+beta).^(k1-2);
end
Num = zeros(k1+1,k1+1);

for j=0:k1
    for i = 0:j
  Num(j+1) = Num(j+1)+ numcoeff(j+1,i+1)*(1+k1-i+beta).^(k1-2);
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
[z,p,k] = zpkdata(system);
disp('zeros of u-transformation')
disp(z{1})
disp('poles of u-transformation')
disp(p{1})
%%%%%%%%%% t-transformation %%%%%%%%%%%

%%%%%% Computation of denominator coefficients that are independent of beta %%%%%%%%%%
for j = 0:k1,
    dencoeff(1,j+1) =(-1)^(k1-j)* (factorial(k1)/(factorial(j)*factorial(k1-j)))/((den(1,(k1-j)+1)));
end

%%%%%% Computation of numerator coefficients that are independent of beta %%%%%%%%%%
for j = 0:k1,
    
    for i = 0:j,
        numcoeff(j+1,i+1) = (-1)^(k1-i)*(factorial(k1)/(factorial(i)*factorial(k1-i)))*((den(j-i+1)))/((den(k1-i+1)));
    end
end

%%%%%%%%%%%%% Computation of numerator and denominator coefficients
%%%%%%%%%%%%% dependent on beta %%%%%%%%%%%%%%%%%%%%%%%%
Den = zeros(1,k1+1);
for i=0:k1
  Den(i+1) = dencoeff(1,i+1)*(1+k1-i+beta).^(k1-1);
end
Num = zeros(k1+1,k1+1);

for j=0:k1
    for i = 0:j
  Num(j+1) = Num(j+1)+ numcoeff(j+1,i+1)*(1+k1-i+beta).^(k1-1);
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
[z,p,k] = zpkdata(system);
disp('zeros of t-transformation')
disp(z{1})
disp('poles of t-transformation')
disp(p{1})
%%%%%%%%%%%%%%%%%% y-transformation %%%%%%%%%%%%

%%%%%% Computation of denominator coefficients that are independent of beta %%%%%%%%%%
for j = 0:k1,
    dencoeff(1,j+1) =(-1)^(k1-j)* (factorial(k1)/(factorial(j)*factorial(k1-j)))/((den(1,(k1-j)+1)));
end

%%%%%% Computation of numerator  coefficients that are independent of beta %%%%%%%%%%
for j = 0:k1,
   
    for i = 0:j,
        numcoeff(j+1,i+1) = (-1)^(k1-i)*(factorial(k1)/(factorial(i)*factorial(k1-i)))*((den(j-i+1)))/((den(k1-i+1)));
    end
end


%%%%%%%%%%%%% Computation of numerator and denominator coefficients
%%%%%%%%%%%%% dependent on beta %%%%%%%%%%%%%%%%%%%%%%%%
Den = zeros(1,k1+1);
for i=0:k1
  Den(i+1) = dencoeff(1,i+1)*pochhammer((1+k1-i+beta),(k1-2));
end

Num = zeros(k1+1,k1+1);

for j=0:k1
    for i = 0:j
  Num(j+1) = Num(j+1)+ numcoeff(j+1,i+1)*pochhammer((1+k1-i+beta),(k1-2));
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
 [z,p,k] = zpkdata(system);
disp('zeros of y-transformation')
disp(z{1})
disp('poles of y-transformation')
disp(p{1})

   %%%%%%%%%%%% tau-transformation %%%%%%%%%%%%%%%%
   
   %%%%%% Computation of denominator coefficients that are independent of beta %%%%%%%%%%
for j = 0:k1,
    dencoeff(1,j+1) =(-1)^(k1-j)* (factorial(k1)/(factorial(j)*factorial(k1-j)))/((den(1,(k1-j)+1)));
end

%%%%%% Computation of numerator coefficients that are independent of beta %%%%%%%%%%
for j = 0:k1,
   
    for i = 0:j,
        numcoeff(j+1,i+1) = (-1)^(k1-i)*(factorial(k1)/(factorial(i)*factorial(k1-i)))*((den(j-i+1)))/((den(k1-i+1)));
    end
end


%%%%%%%%%%%%% Computation of numerator and denominator coefficients
%%%%%%%%%%%%% dependent on beta %%%%%%%%%%%%%%%%%%%%%%%%
Den = zeros(1,k1+1);
for i=0:k1
  Den(i+1) = dencoeff(1,i+1)*pochhammer((1+k1-i+beta),(k1-1));
end

Num = zeros(k1+1,k1+1);

for j=0:k1
    for i = 0:j
  Num(j+1) = Num(j+1)+ numcoeff(j+1,i+1)*pochhammer((1+k1-i+beta),(k1-1));
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
 [z,p,k] = zpkdata(system);
disp('zeros of tau-transformation')
disp(z{1})
disp('poles of tau-transformation')
disp(p{1})
