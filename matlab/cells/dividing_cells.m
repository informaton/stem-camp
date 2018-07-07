function dividing_cells(num_cells, diameter, speed)
% function moving_cells(num_cells, diameter,speed)
% usage moving_cells() - produce an animation of 10 blue cells that move at
% a constant speed of ..
% optional arguments num_cells sets the number of cells to be used and the
% optional parameter speed sets the speed for these cells
%This function was written in attempt to replicate Ingmar's first
%assignment which was an animation of several blue circles (outlines) that
%moved about an axes in random directions at a constant speed
close all; 

global MAX_SIZE
global cell_h;
MAX_SIZE = 2;
delete(timerfindall);
if(nargin<3)
    speed_limit = 0.2;
end
if(nargin<1)
    num_cells = 10;
end
if(nargin<2)
    min_size = 0.2;
    max_size = 0.8;
     diameter = min(rand(num_cells,1)+min_size,max_size); %between 1 and 2;
%     diameter = ones(num_cells,1);
end

half_width = 10;
half_height = 10;
border_width = 20;
border_height = 20;

f = figure('renderer','opengl');%painters may be preferred to set double buffering

a = axes('parent',f,'xlim',[-half_width, half_width],'ylim',[-half_height, half_height],'box','on','xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[]); 

pos = half_width*(2*rand([num_cells,2])-1); %falls in the range of the x and y limits of the axes

% give them all a constant (though random) velocity for now
% velocity = (1-rand([num_cells,2])/2)*speedlimit;
velocity = max(min(randn(num_cells,2),speed_limit),-speed_limit);

%[-1 1]*speed_limit
% from Mathworks.com: linkdata buffers updates to data and dispatches them to plots at roughly
% half-second intervals. This makes data linking not suitable for smoothly
% animating changes in data values unless they are updated in loops that
% are forced to execute two times per second or less.
cell_h = zeros(num_cells,1);
for k=1:num_cells
    cell_h(k) = rectangle('position',[pos(k,:)-diameter(k)/2,diameter(k), diameter(k)],'curvature',[1,1],'parent',a,'edgecolor','blue','linewidth',1,'userdata',velocity(k,:));
end



fps = 20;
t = timer('timerfcn',{@update_cells,border_height, border_width},'period',1/fps,'executionmode','fixeddelay');

set(f,'userdata',t,'closerequestfcn','t = get(gcf,''userdata'');stop(t);delete(t);delete(gcf);'); %clean up the figure and remove timer object from memory
tic; %setup timer in seconds...
start(t);

%alternatively could use this, but don't want to get the wrong one!
%all_timers = timerfind;

function update_cells(hObject,eventdata, border_height, border_width)
%timer update function that adjusts the position of the cells based on the
%speed and whether they are bouncing around or not

%would need a check here in the case that only one cell handle (and the
%cell2mat function would break;
global MAX_SIZE;
global cell_h;
try
      dt = toc;
%      dt = 1e-2;

    pos = cell2mat(get(cell_h,'position'));
    
    velocity = cell2mat(get(cell_h,'userdata'));
    diameter = pos(:,3);
    x = pos(:,1)+diameter/2;
    y = pos(:,2)+diameter/2;
    
    NumCells = numel(cell_h);
    diameter_rep = repmat(diameter,1,NumCells)/2+repmat(diameter',NumCells,1)/2;

    K=1000; % spring constant for harmonic force law
    M=3; %diameter; %min(diameter,3);   % mass of particles

%     pos(:,1:2) = pos(:,1:2)+velocity;
    
    
%     for k=1:numel(cell_h)
%         %handle the border constraints/boundaries
%         if(pos(k,1)+cell_diameter(k)>border_width)
%             pos(k,1)=border_width-cell_diameter(k);
%             velocity(k,1)=-velocity(k,1);
%         elseif(pos(k,1)<-border_width)
%             pos(k,1)=-border_width;
%             velocity(k,1)=-velocity(k,1);
%         end;
%         if(pos(k,2)+cell_diameter(k)>border_height)
%             pos(k,2)=border_height-cell_diameter(k);
%             velocity(k,2)=-velocity(k,2);
%         elseif(pos(k,2)<-border_height)
%             pos(k,2)=-border_height;
%             velocity(k,2)=-velocity(k,2);
%         end;
%         
%         otherCell = hitOtherCell(k,pos,cell_diameter);
%         if(~isempty(otherCell))
%             this_theta = atan2(velocity(k,2),velocity(k,1));
%             other_theta = atan2(velocity(otherCell,2),velocity(otherCell,1));
%             incident_angle = abs(diff(this_theta,other_theta));
%             rotate_by = incident_angle/2;
%             new_velocity = velocity
%         end
%         
%         set(cell_h(k),'position',pos(k,:),'userdata',velocity(k,:));
%     end

    
% Interaction detector and Force Law
  dx=repmat(x,1,NumCells);
  dx=dx-dx';
  dy=repmat(y,1,NumCells);
  dy=dy-dy';
  
  dnm=sqrt(dx.^2+dy.^2)/2+diameter_rep/2; % Distance between all particles
  dnm(1:NumCells+1:end)=diameter;     % set diagonal to diameter

  Fx=K*sum(-(1-diameter_rep./dnm).*dx.*(dnm<diameter_rep),2);  % particle-particle Force Law
  Fy=K*sum(-(1-diameter_rep./dnm).*dy.*(dnm<diameter_rep),2);
  
  Fx=Fx-K*(x-(-border_width/2+diameter/2)).*(x<-border_width/2+diameter/2);  % Left wall
  Fy=Fy-K*(y-(-border_height/2+diameter/2)).*(y<-border_height/2+diameter/2);  % Bottom wall

  Fx=Fx-K*(x-(border_width/2-diameter/2)).*(x>border_width/2-diameter/2);  % Right wall
  Fy=Fy-K*(y-(border_height/2-diameter/2)).*(y>border_height/2-diameter/2);  % Top wall

  ax=Fx./M;
  ay=Fy./M;
  
  % integrator
  velocity(:,1)=velocity(:,1)+ax*dt;
  velocity(:,2)=velocity(:,2)+ay*dt;
  
  x=x+velocity(:,1)*dt;
  y=y+velocity(:,2)*dt;

  growth_rate = .2;
  diameter = diameter+growth_rate*dt;
  

  for k=1:NumCells
     if(diameter(k)>MAX_SIZE)
        %tack on another one....
        x(k) = x(k)-diameter(k)/2;
        y(k) = y(k)-diameter(k)/2;
        diameter(k) = diameter(k)/4;

%         velocity(k) = velocity(k)*0.5;
        cell_h(end+1) = rectangle('position',[[x(k),y(k)]+diameter(k)/2,diameter(k), diameter(k)],'curvature',[1,1],'parent',gca,'edgecolor','blue','linewidth',1,'userdata',velocity(k,:)*-1);

     end
      set(cell_h(k),'position',[[x(k),y(k)]-diameter(k)/2, diameter(k), diameter(k)],'userdata',velocity(k,:));
  end
  
catch ME
   disp(ME); 
end
tic

function otherCell = hitOtherCell(cell_ind,positions,diameters)
%ind is the position of the cell in question and positions is an NumCellsX2
%Matrix whose columns represent the x,y coordinates of each row/cell body
%diameters is the matrix of diameters for each cell

try
x = positions(:,1)+diameters/2; %get to the center of the cell
y = positions(:,2)+diameters/2; %get to the center of the cell


all_ind = (1:size(positions,1))';
cell_x = x(cell_ind);
cell_y = x(cell_ind);

threshold_distance = (diameters+diameters(cell_ind))/2; %figure out the minimum distance between each group

distance_x = abs(cell_x-x);
distance_y = abs(cell_y-y);

hits = distance_x<threshold_distance&distance_y<threshold_distance;

%% alternatively 
% geometric_dist = sqrt((cell_x-x).^2+(cell_y-y).^2);
% hits = geometric_dist<threshold_distance;


otherCell = find(hits&(all_ind~=cell_ind),1);
catch ME
   disp(ME); 
end
