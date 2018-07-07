%% Molecular Dynamics Simulator
% <md.m> Mark D. Shattuck 7/22/2010

% revision history:
% 7/22/2010 Mark D. Shattuck <mds> md.m
%           MD Demo for HandsOn 2010
%           000 mds Initial Conditions and visualization

%% Simulation Parmeters

N=10;  % number of particles
D=2;   % diameter of particles
M=3;   % mass of particles

Lx=10*D;  % size of box
Ly=10*D; 

%% Initial Conditions

[x y]=ndgrid(D/2:D:Lx-D/2,D/2:D:Ly-D/2);  
[junk ii]=sort(rand(1,numel(x)));
x=x(ii(1:N));
y=y(ii(1:N));

%% Setup Plotting
clf;
h=zeros(1,N);
for np=1:N
  h(np)=rectangle('Position',[x(np)-.5*D y(np)-.5*D D D],'Curvature',[1 1],'edgecolor','b');
end
axis('equal');
axis([0 Lx 0 Ly]);

%% Main Loop

for t=1:1000
  for np=1:N
    set(h(np),'Position',[x(np)-.5*D y(np)-.5*D D D]);
  end
  drawnow;
end

