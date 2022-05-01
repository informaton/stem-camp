% my first animation
% author: Hyatt 5
% date : 2/19/2022
N=100;
h=50;

for n=0:N
    
    h=50; y =rand (1,h);x=rand(1,h);  plot(x,2*y,'r-', x, y,'bx','markersize', 2);
    a=gca;
    blue=(1-n/N);
    a.Color = [0 1 blue];
    title(num2str(n))
    pause(.05);  
    
end
