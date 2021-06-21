% Create figures for for use in the presentation
% This mainly adapts figures from the paper to the PowerPoint formatting
% Before:
% * Run case_study/ParRob_nullspace_trajectory.m
%   * Run with option "debug_plot" or manually run the section to create
%     the image with slices of the performance map.
% * Run case_study/ParRob_nullspace_static_pose.m

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2021-06
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear
close all

thispath = fileparts(which('traj_figures_presentation.m'));

%% Modify performance map for trajectory evaluation
% (Fig. 5,a in the paper)
% Load figure of the performance map
uiopen(fullfile(thispath, '..', 'case_study', ['nullspace_traj','.fig']),1);
fh = gcf();
fch = get(fh, 'children');
% Format the plot. Remove lines.
axhdl = gca();
leghdl = legend();
% Delete labels and legend. Will be inserted in PowerPoint
ach = get(axhdl, 'children'); % update variable
delete(ach(strcmp(get(ach,'Type'),'text')));
xlabel(''); ylabel('');
delete(fch(strcmp(get(fch, 'Type'),'legend')));
cb = colorbar();
ylabel(cb,'');
set_font_fontsize(fh,'Times',12)
% Export as image for PowerPoint
set_size_plot_subplot(fh, ...
  11,11,axhdl,...
  0.10,0.16,0.01,0.06,... %l r u d
  0,0) % x y
% Save with all lines
exportgraphics(fh, fullfile(thispath, ['nullspace_traj_alllines','.png']), ...
  'Resolution','300');
ach = get(axhdl, 'children'); % update variable
for k = 1:length(ach)
  % Remove lines other than 3T3R
  if strcmp(get(ach(k), 'Type'), 'line')
    if ~strcmp(get(ach(k), 'DisplayName'), 'ref.')
    	delete(ach(k));
    end
  end
end
% Save with only 3T3R line
exportgraphics(fh, fullfile(thispath, ['nullspace_traj_only3T3Rline','.png']), ...
  'Resolution','300');

% Remove last line and save again with no lines
ach = get(axhdl, 'children'); % update variable
delete(ach(strcmp(get(ach,'Type'),'line')));
exportgraphics(fh, fullfile(thispath, ['nullspace_traj_noline','.png']), ...
  'Resolution','300');

%% Modify Figure for Time Evolution of Performance Criterion in Traj.
% (Fig. 5,b in the paper)
% Load figure of the performance map
uiopen(fullfile(thispath, '..', 'case_study', ['nullspace_traj_condition','.fig']),1);
fh = gcf();
fch = get(fh, 'children');
% Format the plot. Remove lines.
axhdl = gca();
% Delete labels and legend. Will be inserted in PowerPoint
ach = get(axhdl, 'children'); % update variable
delete(ach(strcmp(get(ach,'Type'),'text')));
xlabel(''); ylabel('');
delete(fch(strcmp(get(fch, 'Type'),'legend')));
set_font_fontsize(fh,'Times',12)
% Export as image for PowerPoint
set_size_plot_subplot(fh, ...
  9,9,axhdl,...
  0.12,0.03,0.02,0.09,... %l r u d
  0,0) % x y
exportgraphics(fh, fullfile(thispath, ['nullspace_traj_condition','.png']), ...
  'Resolution','300');

%% Create image with slices of the performance map
% (corresponding to Fig. 5 in the paper, additional material)
uiopen(fullfile(thispath, '..', 'case_study', ['performance_map_slices','.fig']),1);
fh = gcf(); axhdl = gca();
figure_format_publication(fh);
set_font_fontsize(fh,'Times',8)
xlabel('');ylabel('');title('');
set_size_plot_subplot(fh, ...
  5,5,axhdl,...
  0.10,0.02,0.03,0.09,... %l r u d
  0,0) % x y
exportgraphics(fh, fullfile(thispath, ['performance_map_slices','.png']), ...
  'Resolution','300');

%% Create figure for exiting a type-II singularity
% (modified version of Fig. 4c in the paper)
close all;
uiopen(fullfile(thispath, '..', 'case_study', ['pkm_nullspace_case2_overview','.fig']),1);
fh = gcf(); axhdl = gca();
fch = get(fh, 'children');
% Delete everything except subplot for performance criterion
for i = 1:length(fch)
  if strcmp(get(fch(i),'type'),'axes')
    axhdl = fch(i);
  else
    delete(fch(i));
    continue
  end
  if ~strcmp(get(axhdl, 'yscale'), 'log')
    % the last subplot has logscale y axis
    axch = get(axhdl, 'children');
    delete(axch);
    delete(axhdl);
  end
end
axch = get(axhdl, 'children');
delete(axch(strcmp(get(axch,'type'), 'text')));
xlim([-180, 180])
ylim([30, 1e4])
delete(findall(fh,'type','annotation'))
xlabel(''); ylabel('');
set(axhdl, 'xtick', -180:45:180)
set_font_fontsize(fh,'Times',12)
set_size_plot_subplot(fh, ...
  8,8,axhdl,...
  0.10,0.04,0.03,0.12,... %l r u d
  0,0) % x y
exportgraphics(fh, fullfile(thispath, ['pkm_nullspace_case2_criterion','.png']), ...
  'Resolution','300');