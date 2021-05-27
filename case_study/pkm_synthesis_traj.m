% Perform a dimensional synthesis of the hexapod parallel robot
% Determine the kinematics parameters that give good conditioning for the
% trajectory.

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2021-05
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear

% Task DoF
DoF = [1 1 1 1 1 0];
Set = cds_settings_defaults(struct('DoF', DoF));
%% Initialize Trajectory
Set.task.Ts = 5e-2;
Set.task.Tv = 1e-1;
Set.task.amax = 2;
% Define intermediate poses to interpolate
% Rechteck-Trajektorie
d1=0.15; d2=0.2; % dimensions of the rectangle
ta = pi/6; % tilting angle, point from outside on edge of rectangle
X0 = [ [0.00;0.00;0.6]; [0;0;0]*pi/180 ];
% First corner without tilting
k=1;   XE(k,:) = X0' + [-d1/4,+d2/4,+0.1,0,0,0];
% Edge 1: Move in x direction
k=k+1; XE(k,:) = XE(1,:) + [ 0,0,0, r2eulxyz(rotx(ta))'];
k=k+1; XE(k,:) = XE(1,:) + [ d1,0,0, r2eulxyz(rotx(ta))'];
% Edge 2: Move in -y direction
k=k+1; XE(k,:) = XE(1,:) + [ d1, 0,0, r2eulxyz(rotz(-pi/2)*rotx(ta))'];
k=k+1; XE(k,:) = XE(1,:) + [ d1,-d2,0, r2eulxyz(rotz(-pi/2)*rotx(ta))'];
% Edge 3: Move in -x direction
k=k+1; XE(k,:) = XE(1,:) + [ d1, -d2,0, r2eulxyz(rotz(pi)*rotx(ta))'];
k=k+1; XE(k,:) = XE(1,:) + [ 0,  -d2,0, r2eulxyz(rotz(pi)*rotx(ta))'];
% Edge 4: Move in +y direction
k=k+1; XE(k,:) = XE(1,:) + [ 0,  -d2,0, r2eulxyz(rotz(pi/2)*rotx(ta))'];
k=k+1; XE(k,:) = XE(1,:) + [ 0,  0,0, r2eulxyz(rotz(pi/2)*rotx(ta))'];

[X_ges,XD_ges,XDD_ges,T_ges] = traj_trapez2_multipoint(XE, ...
  Set.task.vmax, Set.task.vmax/Set.task.amax, Set.task.Tv, Set.task.Ts, 0);
Traj = struct('X', X_ges, 'XD', XD_ges, 'XDD', XDD_ges, 't', T_ges, 'XE', XE);

%% Settings for Robot Synthesis
Set.optimization.objective = {'condition'};
Set.optimization.optname = 'hexapod_icinco_20210506_v1';
Set.optimization.NumIndividuals = 50;
Set.optimization.MaxIter = 20;
Set.optimization.obj_limit = 1e3; % Bei erstem Erfolg aufhören
Set.general.verbosity = 2;
Set.optimization.constraint_collisions = false;
Set.general.eval_figures = {};
Set.general.animation_styles = {};
Set.general.plot_details_in_fitness = 0;%4e9; % Debug-Plots für Selbstkollision
% Dimensions of the parallel robot. Fix to 600mm/200mm
Set.optimization.base_size = false;
Set.optimization.base_size_limits = [0.6, 0.6];
Set.optimization.platform_size = false;
Set.optimization.platform_size_limits = [0.2, 0.2];
Set.optimization.base_morphology = true;
Set.optimization.platform_morphology = true;
Set.optimization.movebase = true;
Set.optimization.basepos_limits = repmat([-0.1, 0.1],3,1);
Set.optimization.rotate_base = false;
Set.optimization.ee_rotation = false;
Set.optimization.ee_translation = false;
% Set.optimization.max_range_prismatic = 1.2; % according to the dimensions in the paper example
% Perform combined structural and dimensional synthesis. Parallel robots
% with 3T3R platform DoF and UPS joint structure.
Set.structures.use_serial = false;
Set.structures.min_task_redundancy = 1;
Set.structures.max_task_redundancy = 1;
Set.structures.joint_filter = 'RRPRRR';
Set.structures.num_tech_joints = 3;
Set.structures.mounting_parallel = 'floor';
Set.structures.whitelist = {'P6RRPRRR14V3G1P4A1'}; % Debug. Only one robot
cds_start