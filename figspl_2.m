
%%%%%% Location of Poles on the complex plane as n is varied for a 5th order
%%%%%% Gaussian wavelet obtained with u-transformation with a translation of 2, sigma 1%%%%%%%%


clear all
syms s

%%%%% Laplace tranform of gaussian wavelet with t0=2,sigma=1 %%%%%%%%%%
q = -5.0*10^(-13)*exp(-0.25*(s - 4.0)^2)*(36631277799.0*exp(0.25*(s - 4.0)^2) + 3.544907702*exp(0.5*(s - 4.0)^2)*erfc(0.5*s - 2.0) - 32463624699.0*s*exp(0.5*(s - 4.0)^2)*erfc(0.5*s - 2.0));

%%%%%% Taylor series of the laplace transform around 0 %%%%%%
t =  taylor(q,'order',33,'ExpansionPoint',0);
 a= sym2poly(t);
den = fliplr(a);


%%%%% Order of the approximation
k1=5;
%%%%%% constant shift parameter %%%%
beta=0;

%%%%%%%%%%%%%% Denoninator polynomial with n=1 %%%%%%%%%%%%%%%%%%
for j = 0:k1,
    dencoeff(1,j+1) =(-1)^(k1-j)* (factorial(k1)/(factorial(j)*factorial(k1-j)))/((den(1,(k1-j)+1)));
end



Den = zeros(1,k1+1);
for i=0:k1
  Den(i+1) = dencoeff(1,i+1)*(1+k1-i+beta).^(k1-2);
end

mydenom = [];
for j=1:k1+1
    mydenom = [mydenom,Den(k1+1-(j-1))];
end

 system=tf(1,mydenom);
 %%%%% n=1 %%%%%%
pzplot(system);
hold on

%%%%%%%%%%%%%% Denoninator polynomial with n=2 %%%%%%%%%%%%%%%%%%
for j = 0:k1,
    dencoeff(1,j+1) =(-1)^(k1-j)* (factorial(k1)/(factorial(j)*factorial(k1-j)))/((den(1,(k1-j)+2)));
end



Den = zeros(1,k1+1);
for i=0:k1
  Den(i+1) = dencoeff(1,i+1)*(2+k1-i+beta).^(k1-2);
end

mydenom = [];
for j=1:k1+1
    mydenom = [mydenom,Den(k1+1-(j-1))];
end

 system=tf(1,mydenom);
pzplot(system)

%%%%%%%%%%%%%% Denoninator polynomial with n=3 %%%%%%%%%%%%%%%%%%

for j = 0:k1,
    dencoeff(1,j+1) =(-1)^(k1-j)* (factorial(k1)/(factorial(j)*factorial(k1-j)))/((den(1,(k1-j)+3)));
end



Den = zeros(1,k1+1);
for i=0:k1
  Den(i+1) = dencoeff(1,i+1)*(3+k1-i+beta).^(k1-2);
end

mydenom = [];
for j=1:k1+1
    mydenom = [mydenom,Den(k1+1-(j-1))];
end

 system=tf(1,mydenom);
pzplot(system)


legend('n=1','n=2','n=3')
% opt1=[];
% % opt1.LineStyle = {'-'};
% opt1.FontSize=11;
% opt1.Colors = [ 
%     0.15, 0.15, 0.15;
%     0.15,0.15,0.15;
%     0.15, 0.15, 0.15;
%     0.15, 0.15, 0.15;
%     0.15, 0.15, 0.15;
%     0.15,0.15,0.15;
%     0.15, 0.15, 0.15;
%     0.15, 0.15, 0.15;
%     0.15, 0.15, 0.15;
%     0.15,0.15,0.15;
%     0.15, 0.15, 0.15;
%     0.15, 0.15, 0.15;
%    0.15,0.15,0.15;
%     0.15, 0.15, 0.15;   
%     0.15, 0.15, 0.15;
%     ];
% opt1.Markers = {'s','s','s','s','s','+','+','+','+','+','o','o','o','o','o'};
% opt1.XLabel='Real Axis';
% opt1.YLabel='Imaginary Axis'
% 
% 
% opt1.BoxDim=[7 3.0];
% opt1.XGrid = 'off';
% opt1.YGrid = 'off';
% opt1.FontName = 'Times New Roman';
% opt1.FontSize = [12];
% opt1.XMinorTick = 'on';
% opt1.YMinorTick = 'on';
% opt1.AxisLineWidth = [1];
% opt1.ShowBox = 'on';
% opt1.LegendBox = 'on';
% opt1.LegendLoc ='North';
% setPlotProp(opt1);





