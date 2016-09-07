
clear all 

syms s t b

Coeff = [];
%%%%% computation of Taylor coefficients of Gaussian Wavelet (First diff of
%%%%% Gaussian function) as the translation is varied %%%%%%%%%%%%%
for b = 0:0.25:3.25
    
    f = exp(-(t - b).^2);
    g = diff(f);
    h = laplace(g,t,s);
    tay =  taylor(h,'order',6,'ExpansionPoint',0);
    a= sym2poly(tay);
    coeff = fliplr(a);
    Coeff = vertcat(Coeff,coeff);     
end

%%%%%%%%%% Computationf of absolute percentage change of a0 - a5 %%%%%
pdiff=[]; final = []; 
for j = 1:6
         for i = 1:13
                pd= abs((Coeff(i,j)-Coeff(i+1,j))/Coeff(i+1,j))*100;
                pdiff=horzcat(pdiff,pd);
          end

final = vertcat(final,pdiff);
pdiff = [];
end 
%%%%%%%%%% Plotting the Absolute percentage change vs translation %%%%%%
t1 = 0.25:0.25:3.25;
for k = 1:6
    
    plot(t1,final(k,:));
     hold on
end

%%%%plot properties %%%%%%%%%%%
    legend('a0','a1','a2','a3','a4','a5','Location','northwest');
     axis([0 3 0 300]);
     xlabel('Translation  (s)') % x-axis label
     ylabel('Absolute Percentage Change') % y-axis label
     title('Absolute % change of Taylor coefficents a0 - a5 of Gaussian wavelet as translation is increased')
    