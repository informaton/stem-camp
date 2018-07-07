function y = f(p)
close all;
%initial conditions
T = 10; %10 samples;
po = 1;
p = 20; %number of protein molecules
m = 10;  %number of mRNA molecules

p  = 1;
k = 2.5;
po = 1;
n = 2;
% y = k./(1+p.^n./po.^n);

iters = 10000
y = zeros(iters,1);

dt = 1e-2;

%decay constants
b = 1;
c = 1;
a = 1; %production rate of new protein
for i = 1:iters
    
    dp = a*m-b*p;    
    dm = k./(1+p.^n./po.^n) - c*m;
    m = m+dm*dt;
    p = p+dp*dt;
    y(i) = k./(1+p.^n./po.^n);

    
end

plot(1:iters,y);
