function moving_cells(num_cells, diameter, speed)
% function moving_cells(num_cells, diameter,speed)
% usage moving_cells() - produce an animation of 10 blue cells that move at
% a constant speed of ..
% optional arguments num_cells sets the number of cells to be used and the
% optional parameter speed sets the speed for these cells
%This function was written in attempt to replicate Ingmar's first
%assignment which was an animation of several blue circles (outlines) that
%moved about an axes in random directions at a constant speed
close all; 
if(nargin<3)
    speed = 0.1;
end
if(nargin<2)
    diameter = 1;
end
if(nargin<1)
    num_cells = 10;
end

half_width = 10;
half_height = 10;
f = figure('renderer','opengl');%painters may be preferred to set double buffering

a = axes('parent',f,'xlim',[-half_width, half_width],'ylim',[-half_height, half_height],'box','on','xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[]); 

pos = half_width*(2*rand([num_cells,2])-1); %falls in the range of the x and y limits of the axes

% give them all a constant (though random) velocity for now
velocity = (0.5-rand([num_cells,2]))*speed;  

% from Mathworks.com: linkdata buffers updates to data and dispatches them to plots at roughly
% half-second intervals. This makes data linking not suitable for smoothly
% animating changes in data values unless they are updated in loops that
% are forced to execute two times per second or less.
h = zeros(num_cells,1);
for k=1:num_cells
    h(k) = rectangle('position',[pos(k,:),diameter, diameter],'curvature',[1,1],'parent',a,'edgecolor','blue','linewidth',1,'userdata',velocity(k,:));    
end


fps = 20;
t = timer('timerfcn',{@update_cells,h,diameter,half_height, half_width},'period',1/fps,'executionmode','fixeddelay');

set(f,'userdata',t,'closerequestfcn','t = get(gcf,''userdata'');stop(t);delete(t);delete(gcf);'); %clean up the figure and remove timer object from memory
start(t);

%alternatively could use this, but don't want to get the wrong one!
%all_timers = timerfind;

function update_cells(hObject,eventdata,cell_h, cell_diameter, border_height, border_width)
%timer update function that adjusts the position of the cells based on the
%speed and whether they are bouncing around or not

%would need a check here in the case that only one cell handle (and the
%cell2mat function would break;
pos = cell2mat(get(cell_h,'position'));
velocity = cell2mat(get(cell_h,'userdata'));

pos(:,1:2) = pos(:,1:2)+velocity;

for k=1:numel(cell_h)
    %handle the border constraints/boundaries
    if(pos(k,1)+cell_diameter>border_width)
        pos(k,1)=border_width-cell_diameter;
        velocity(k,1)=-velocity(k,1);
    elseif(pos(k,1)<-border_width)
        pos(k,1)=-border_width;
        velocity(k,1)=-velocity(k,1);
    end;
    if(pos(k,2)+cell_diameter>border_height)
        pos(k,2)=border_height-cell_diameter;
        velocity(k,2)=-velocity(k,2);
    elseif(pos(k,2)<-border_height)
        pos(k,2)=-border_height;
        velocity(k,2)=-velocity(k,2);
    end;
    
    set(cell_h(k),'position',pos(k,:),'userdata',velocity(k,:));
end
