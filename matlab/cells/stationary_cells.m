function stationary_cells(num_cells, diameter, speed)
% function stationary_cells(num_cells, diameter,speed)
% 
% moving_cells() - produce an animation of 10 blue cells that move at
% a constant speed of ..
% optional arguments num_cells sets the number of cells to be used and the
% optional parameter speed sets the speed for these cells
%This function was written in attempt to replicate Ingmar's first
%assignment which was an animation of several blue circles (outlines) that
%moved about an axes in random directions at a constant speed
close all; 
if(nargin<3)
    speed = 0; %good velocities can be found between +/- 0.1
end
if(nargin<2)
    diameter = 1;
end
if(nargin<1)
    num_cells = 10;
end

%% initial graphics setup
half_width = 10;
half_height = 10;
fig_h = figure('renderer','opengl');%painters may be preferred to set double buffering
axes_h = axes('parent',fig_h,'xlim',[-half_width, half_width],'ylim',[-half_height, half_height],'box','on','xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[]); 

pos = half_width*(2*rand([num_cells,2])-1); %falls in the range of the x and y limits of the axes


% Programming note on why we do not use Matlab's linkdata functionality
% from Mathworks.com: linkdata buffers updates to data and dispatches them to plots at roughly
% half-second intervals. This makes data linking not suitable for smoothly
% animating changes in data values unless they are updated in loops that
% are forced to execute two times per second or less.

%h is a an array of handles that keep track of the graphic objects created
%here
h = zeros(num_cells,1);
for k=1:num_cells
    userdata = [];
    h(k) = rectangle('position',[pos(k,:),diameter, diameter],'curvature',[1,1],...
        'parent',axes_h,'edgecolor','blue','linewidth',1,'userdata',userdata);
end

%frames per second
fps = 20;
t = timer('period',1/fps,'executionmode','fixeddelay');

% This portion may need updating to pass more parameters
% through to the update_cells function (hObject and eventdata are 
% automatically added and do not require declaration here
set(t,'timerfcn',{@update_cells,h});

%clean up the figure and remove timer object from memory
set(fig_h,'userdata',t,'closerequestfcn','t = get(gcf,''userdata'');stop(t);delete(t);delete(gcf);'); 
%alternatively could do this:
%   all_timers = timerfind; stop(all_timers);

start(t);



function update_cells(hObject,eventdata,cell_h)
%timer update function that is caled at a fixed delay rate of 1/fps seconds.
% hObject is a handle to the parent callback object (i.e. the timer)
% eventdata is not used in this case, but is required by Matlab as a place
% holder for the second argument.
%
%  This function needs to  be updated so that the cells' positions are 
%  updated on each call.
%

pos = get(cell_h,'position');
if(iscell(pos))
    pos = cell2mat(pos);
end

%loop through the cells - currently they remain stationary
for k=1:size(pos,1)
    userdata = get(cell_h(k),'userdata');
   set(cell_h(k),'position',pos(k,:)); 
end

