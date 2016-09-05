clear all
syms s

%%%%% Laplace tranform of MExican Hat wavelet with t0=2,sigma=1 %%%%%%%%%%
q=2.5*10^(-13)*exp(-0.25*(s - 4.0)^2)*(1.298544987*10^11*exp(0.5*(s - 4.0)^2)*erfc(0.5*s - 2.0) - 73262555566.0*s*exp(0.25*(s - 4.0)^2) - 2.930502222*10^11*exp(0.25*(s - 4.0)^2) + 7.089815404*s*exp(0.5*(s - 4.0)^2)*erfc(0.5*s - 2.0) + 64927249377.0*s^2*exp(0.5*(s - 4.0)^2)*erfc(0.5*s - 2.0)) - 0.01831563889*pi^(1/2)*erfc(s/2 - 2.0)*exp((s - 4.0)^2/4);

%%%%%% Taylor series of the laplace transform around 0 %%%%%%
t =  taylor(q,'order',33,'ExpansionPoint',0);
 a= sym2poly(t);
den = fliplr(a);


%%%%% Order of the approximation
k1=5;
%%%%%% constant shift parameter %%%%
beta=0;
%%%%% u-transformation compuation of denominator polynomial%%%%%%%%%
for j = 0:k1,
    dencoeff(1,j+1) =(-1)^(k1-j)* (factorial(k1)/(factorial(j)*factorial(k1-j)))/((den(1,(k1-j)+1)));
end

%%%%% u-transformation compuation of numerator polynomial%%%%%%%%%
for j = 0:k1,
    
    for i = 0:j,
        numcoeff(j+1,i+1) = (-1)^(k1-i)*(factorial(k1)/(factorial(i)*factorial(k1-i)))*((den(j-i+1)))/((den(k1-i+1)));
    end
end

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
%%%%% transfer function obtained from u-transformation %%%%%%%%%%
 system=tf(mynum,mydenom);
[z,p,k] = zpkdata(system);
disp('zeros of u-transformation')
disp(z{1})
disp('poles of u-transformation')
disp(p{1})

for j = 0:k1,
    dencoeff(1,j+1) =(-1)^(k1-j)* (factorial(k1)/(factorial(j)*factorial(k1-j)))/((den(1,(k1-j)+1)));
end


for j = 0:k1,
    
    for i = 0:j,
        numcoeff(j+1,i+1) = (-1)^(k1-i)*(factorial(k1)/(factorial(i)*factorial(k1-i)))*((den(j-i+1)))/((den(k1-i+1)));
    end
end

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

for j = 0:k1,
    dencoeff(1,j+1) =(-1)^(k1-j)* (factorial(k1)/(factorial(j)*factorial(k1-j)))/((den(1,(k1-j)+1)));
end

for j = 0:k1,
   
    for i = 0:j,
        numcoeff(j+1,i+1) = (-1)^(k1-i)*(factorial(k1)/(factorial(i)*factorial(k1-i)))*((den(j-i+1)))/((den(k1-i+1)));
    end
end

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


for j = 0:k1,
    dencoeff(1,j+1) =(-1)^(k1-j)* (factorial(k1)/(factorial(j)*factorial(k1-j)))/((den(1,(k1-j)+1)));
end

for j = 0:k1,
   
    for i = 0:j,
        numcoeff(j+1,i+1) = (-1)^(k1-i)*(factorial(k1)/(factorial(i)*factorial(k1-i)))*((den(j-i+1)))/((den(k1-i+1)));
    end
end

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
