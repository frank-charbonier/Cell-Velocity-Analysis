function plot_cell_trajectories(traj_name, nstart, nend, savename, invisible)
    arguments
        % Name of data containing cell trajectories
        traj_name = 'cell_trajectories_tstart_end.mat';
        % Choose time points to use in the analysis. Select timepoints as a
        % fraction as the total number of time points available, where 0
        % corresponds to the first time point and 1 corresponds to the last time
        % point.
        nstart = 0;
        nend = 1;
        % Name to save plot 
        savename = "Trajectories";
        % Set to [] to make figure visible. Set to 1 to make figure invisible.
        invisible = 1;
    end
% 
% Plot cell trajectories. Run this script after running
% compute_cell_trajectories.m.
% 
% If running as a batch, uncomment the statement function at the top and
% comment the clear command
%
% Written by Jacob Notbohm, University of Wisconsin-Madison, 2021

% clear;
close all;
% clc;


%% --- USER INPUTS ---

% Name of data containing cell trajectories
% traj_name = 'cell_trajectories_tstart_end.mat';
% Option to downsample number of trajectories when plotting them. Typically
% only do this if the trajectories are too dense to see. This value must be
% a positive integer. Be carefule here, because if the script
% compute_cell_trajectories.m is run correctly, then the number of data
% points approximately equals the number of cells, meaning that
% downsampling will plot fewer trajectories than there are cells, which may
% be misleading. To avoid downsampling, enter a value of 1.
fd = 1;
% Set axis limits. To set the limits automatically, set to empty array [].
% Format [xmin xmax ymin ymax]. Units: um
xy_lim = []; 


%% --- PLOT TRAJECTORIES ---

load(traj_name);
% Units: um

hf = make_fig([0.5 0.5 1 1]);
% Set axes to nearly fill window
set(hf,'DefaultAxesPosition',[0.1 0.1 .86 .86]);
if invisible
        set(hf,'visible','off');
end
hold on

% --- Get time points of interest ---
K = size(traj_x,2);
idx = round(nstart*K)+1 : round(nend*K);
traj_x = traj_x(:,idx);
traj_y = traj_y(:,idx);

% Downsample for plotting
traj_x = downsample(traj_x,fd);
traj_y = downsample(traj_y,fd);

% % Plot trajectories with the same pre-selected color (not typically used)
% plot(traj_x',traj_y','b')
% Or, more common is...
% Let Matlab automatically choose color. Note that colors are randomized
% and repeated. The different colors are only to help identify the
% different trajectories.
plot(traj_x',traj_y')
xlabel('\mum','fontsize',11);
ylabel('\mum','fontsize',11);
axis equal;
set(gca,'box','on','fontsize',11);
% Axis limits
if ~isempty(xy_lim)
    axis(xy_lim);
else
    % Otherwise, get the automatically
    x1 = min(traj_x(:));
    x2 = max(traj_x(:));
    y1 = min(traj_y(:));
    y2 = max(traj_y(:));
    axis([x1 x2 y1 y2])
end

% Option to set ticks
% set(gca,'xtick',[-500:250:500],'ytick',[-500:250:500]);

set(hf,'Paperpositionmode','auto');
full_savename = savename + "_" + nstart + "-" + nend + ".png";
print('-dpng','-r300',full_savename);

% Option: Save as eps
% print('-depsc',savename);

% % Option: Save plot in a different directory--useful when comparing
% % trajectories from multiple positions. Uncomment the two lines below and
% % make sure the path exists
% [~,curdir,~] = fileparts(pwd);
%  print('-dpng','-r300',['../Trajectories/',curdir]);





