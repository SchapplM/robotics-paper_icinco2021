% Perform nullspace motion of a hexapod robot. This creates the results of
% Sec. 6.2 of the paper (added after the initial submission).
% Figures are saved into the directory paper/figures and update latex doc.

% This file is adapted from the example ParRob_benchmark_3T3R_task_redundancy.m
% located at examples_tests/ParRob in this Git repository:
% https://github.com/SchapplM/robotics-toolbox

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-07 - 2021-05
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clear
clc

if isempty(which('serroblib_path_init.m'))
  error('The serial robot database is not initialized in Matlab.');
end
if isempty(which('parroblib_path_init.m'))
  error('The parallel robot database is not initialized in Matlab.');
end

%% User Settings and Initialization
%#ok<*UNRCH>
use_mex_functions = true; % use mex functions (much faster)
usr_recreate_mex = false; % recreate and recompile mex functions from templates
usr_short_traj = false; % Trajektorie stark abkürzen, um prinzipielle Funktionalität zu zeigen
usr_discretization_type = 'position-level'; % alternative 'trajectory'
usr_create_anim = false; % create an animation video of the robot motion
usr_anim_realtime = false; % save real-time animation (video as long as trajectory in seconds)
usr_highres_distrfig = true; % high resolution of the paper figure for performance criterion map
debug_plot = false;% Create debug plots
usr_save_figures = true; % save figures to disk
usr_load_discretization = true; % load a previously computed performance map (if trajectory stays the same)
usr_load_traj = false; % load a previously computed joint space trajectory
respath = fileparts(which('ParRob_nullspace_trajectory.m'));
paperfig_path = fullfile(respath, '..', 'paper', 'figures');
usr_only_test_keypoints = false; % if true: Abort script after key point computation
usr_only_global_discretization = false; % if true: Abort script after performance map computation
usr_test_timing = true; % multiple calls of ik functions to determine runtime for evaluation in Sec. 6.1.3.

% Formatierung der Linien der verschiedenen Methoden
format_mlines = { 'r', 'v', '-', 8; ...
                  'g', 'd', '-', 5; ...
                  'b', 's', '--', 7; ...
                  'k', 'x', '--', 9; ...
                  'm', 'o', ':', 6; ...
                  'c', '^', '-', 3; ...
                  'r', '+', ':', 6};
optimcrit_limits_hyp_deact = 0.9;  % Hyperbolische Funktion nur nahe an Grenzen
%% Klasse für PKM erstellen (basierend auf serieller Beinkette)
% Robotermodell aus PKM-Bibliothek laden.
robnr = 2; % from old test script
%% Klasse für PKM erstellen (basierend auf PKM-Bibliothek)
RP = parroblib_create_robot_class('P6RRPRRR14V3G1P4A1', 0.6, [0.150;0.100]);
T_W_0 = trotx(pi)*transl([0;0;1]); % PKM is on ceiling pointing down
RP.update_base(T_W_0(1:3,4), r2eulxyz(T_W_0(1:3,1:3)));
% Insert design parameters (for Plot)
for i = 1:RP.NLEG
  RP.Leg(i).DesPar.joint_type(RP.I_qa((RP.I1J_LEG(i):RP.I2J_LEG(i)))) = 5;
  % hollow cylinder with diameter 50mm
  RP.Leg(i).DesPar.seg_par=repmat([5e-3,50e-3],RP.Leg(i).NL,1);
end
RP.DesPar.platform_par(end) = 5e-3;
if usr_recreate_mex
  serroblib_create_template_functions({RP.Leg(1).mdlname}, false, false);
  parroblib_create_template_functions({RP.mdlname(1:end-2)}, false, false);
  matlabfcn2mex({[RP.Leg(1).mdlname,'_invkin_eulangresidual']});
  matlabfcn2mex({[RP.Leg(1).mdlname,'_invkin_traj']});
  matlabfcn2mex({[RP.mdlname(1:end-6), '_invkin']});
  matlabfcn2mex({[RP.mdlname(1:end-6), '_invkin3']});
  matlabfcn2mex({[RP.mdlname(1:end-6), '_invkin_traj']});
end
RP.fill_fcn_handles(use_mex_functions,true);
% Definition der Freiheitsgrade (vollständig und reduziert)
I_EE_full = RP.I_EE;
I_EE_red = logical([1 1 1 1 1 0]);
I_EE_full_str = sprintf('%dT%dR', sum(I_EE_full(1:3)), sum(I_EE_full(4:6)));
I_EE_red_str = sprintf('%dT%dR', sum(I_EE_red(1:3)), sum(I_EE_red(4:6)));
for i = 1:RP.NLEG
  % Begrenze die Winkel der Kugel- und Kardangelenke auf +/- 360°
  RP.Leg(i).qlim = repmat([-2*pi, 2*pi], RP.Leg(i).NQJ, 1);
  % Begrenze die Länge der Schubgelenke
  % Sieht auf Bild komisch aus, wenn die Maximallänge sehr groß ist. Dann
  % muss der Hubzylinder entsprechend lang sein. Je enger die Grenzen
  % sind, desto weniger darf ein Überschwingen stattfinden.
  qpris_minmax = [0.6, 1.2]; % TODO: Wert reduzieren und bessere Grenzüberwachung
  RP.Leg(i).qlim(RP.Leg(i).MDH.sigma==1,:) = ...
    repmat(qpris_minmax, sum(RP.Leg(i).MDH.sigma==1), 1);
  % Setze Geschwindigkeit und Beschleunigung auf moderate Werte, damit
  % die Schrittweite in der Traj.-IK nicht zu groß wird und die Klein-
  % winkelnäherung nicht verletzt wird.
  RP.Leg(i).qDlim = repmat([-45, 45]*pi/180, RP.Leg(i).NQJ, 1); % 45deg/s is moderately high
  RP.Leg(i).qDlim(RP.Leg(i).MDH.sigma==1,:) = ...
    repmat([-2, 2], sum(RP.Leg(i).MDH.sigma==1), 1); % 2m/s
  % Set moderately high acceleration limits which could be feasible for a
  % high-performance parallel robot.
  RP.Leg(i).qDDlim = repmat([-20, 20], RP.Leg(i).NQJ, 1);
  RP.Leg(i).qDDlim(RP.Leg(i).MDH.sigma==1,:) = ...
    repmat([-20, 20], sum(RP.Leg(i).MDH.sigma==1), 1);
end
%% Startpose bestimmen
% Mittelstellung im Arbeitsraum
X0 = [ [0.00;0.00;0.5]; [0;0;0]*pi/180 ];
q0 = 0.5+rand(RP.NJ,1); % Startwerte für numerische IK (zwischen 0.5 und 1.5 rad)
q0(RP.MDH.sigma==1) = 0.5; % mit Schubaktor größer Null anfangen (damit Konfiguration nicht umklappt)
[qs, Phis] = RP.invkin_ser(X0, rand(RP.NJ,1), struct('retry_on_limitviol',true));
if any(abs(Phis) > 1e-6)
  error('Inverse Kinematik (für jedes Bein einzeln) konnte in Startpose nicht berechnet werden');
end
if any(qs(RP.MDH.sigma==1) < 0)
  error('Start-Konfiguration ist umgeklappt mit Methode Seriell.');
end
% TODO: Gelenkgrenzen nicht anpassen. Ist wesentlich aufwändiger in IK
% mit strengen Grenzen. Nehme die Gelenkkonfig. in der Mittelstellung des
% Roboters.
for i = 1:RP.NLEG
  q0_i = qs(RP.I1J_LEG(i):RP.I2J_LEG(i));
  % TODO: Begrenze die Winkel der Kugel- und Kardangelenke auf +/- 30°
  % TODO: Geht erst, wenn Grenzen richtig eingehalten werden können.
  % Ignoriere Grenzen zunächst. Für Paper-Ergebnisse nicht kritisch
  RP.Leg(i).qlim(1:2,:) = q0_i(1:2) + 10*repmat([-360, 360]*pi/180, 2, 1);
  RP.Leg(i).qlim(4:6,:) = q0_i(4:6) + 10*repmat([-360, 360]*pi/180, 3, 1);
end
% qlim für gesamte PKM festlegen
qlim   = cat(1, RP.Leg.qlim);
qDlim  = cat(1, RP.Leg.qDlim);
qDDlim = cat(1, RP.Leg.qDDlim);
%% Initialisierung Teil 2
% Roboter auf 3T2R einstellen
RP.update_EE_FG(I_EE_full, I_EE_red);

%% Eckpunkte für Beispiel-Trajektorie bestimmen und IK prüfen
% IK-Grundeinstellungen
s = struct('Phit_tol', 1e-12, 'Phir_tol', 1e-12, ... % sehr genau rechnen
  'maxstep_ns', 1e-5, ... % Schrittweite für Nullraum-Inkremente gering halten
  'wn', [0;1], ... % keine Vorgabe von K oder Kn (Standard-Werte)
  'scale_lim', 0.7, ...
  'retry_limit', 0);
% Rechteck-Trajektorie
d1=0.25; d2=0.2; % dimensions of the rectangle
ta = pi/4; % tilting angle, point from outside on edge of rectangle
% First corner without tilting
X1 = X0' + [-d1/5,+d2/5,+0.2,0,0,0];
% Edge 1: Move in x direction
k=1;   XL(k,:) = X1 +      [ 0,0,0, r2eulxyz(rotx(ta))']; % first point
k=k+1; XL(k,:) = X1 + [ d1,0,0, r2eulxyz(rotx(ta))'];
% Edge 2: Move in -y direction
k=k+1; XL(k,:) = X1 + [ d1, 0,0, r2eulxyz(rotz(-pi/2)*rotx(ta))'];
k=k+1; XL(k,:) = X1 + [ d1,-d2,0, r2eulxyz(rotz(-pi/2)*rotx(ta))'];
% Edge 3: Move in -x direction
k=k+1; XL(k,:) = X1 + [ d1, -d2,0, r2eulxyz(rotz(pi)*rotx(ta))'];
k=k+1; XL(k,:) = X1 + [ 0,  -d2,0, r2eulxyz(rotz(pi)*rotx(ta))'];
% Edge 4: Move in +y direction
k=k+1; XL(k,:) = X1 + [ 0,  -d2,0, r2eulxyz(rotz(pi/2)*rotx(ta))'];
k=k+1; XL(k,:) = X1 + [ 0,  0,0, r2eulxyz(rotz(pi/2)*rotx(ta))'];

% Determine tilting angles
TA = NaN(size(XL,1), 1);
for i = 1:size(XL,1)
  R_i = eulxyz2r(XL(i,4:6)');
  TA(i) = acos(R_i(:,3)'*[0;0;1]);
end
% Debug:
% 180/pi*TA

% Zeitverlauf der Trajektorie generieren
X_t = traj_trapez2_multipoint(XL, 3, 0.05, 0.01, 2e-3, 0);

% Debug: Zeichne Roboter mit Trajektorie
if debug_plot
  change_current_figure(100);clf;
  set(100,'Name','Rob','NumberTitle','off');
  title(sprintf('Robot in Rest Pose'));
  hold on;grid on;
  xlabel('x in m');ylabel('y in m');zlabel('z in m');
  view(3);
  s_plot = struct(  'ks_legs', [], 'straight', 0, 'mode', 4);
  RP.plot( qs, X0, s_plot );
  Traj_0 = struct('X', XL);
  Traj_W = RP.transform_traj(Traj_0);
  plot3(Traj_W.X(:,1), Traj_W.X(:,2), Traj_W.X(:,3), 'k-');
  Xii_lastplot = inf(1,6);
  for ii = 1:size(X_t,1)
    T_ii = RP.x2t(X_t(ii,:)');
    % Reduce density of the coordinate systems in the plot
    if ii > 1 && norm(X_t(ii,1:3)-Xii_lastplot(1:3))<0.02, continue; end
    Xii_lastplot = X_t(ii,:);
    trplot(T_W_0*T_ii, 'frame', '', 'rgb', 'length', 0.02, 'labels', '   ')
  end
  xlabel('x in mm');ylabel('y in mm');zlabel('z in mm');
  axis equal
  if usr_save_figures
    saveas(100, fullfile(respath, 'robot_with_traj.fig'));
    saveas(100, fullfile(respath, 'robot_with_traj.png'));
  end
end

%% Berechne IK zu den einzelnen Eckpunkten (zum Testen)
% Berechne IK zu den einzelnen Eckpunkten. Wenn das nicht geht, bringt die
% Trajektorie sowieso nichts. Benutze die IK mit Aufgabenredundanz 
s_ep = s; % Einstellungen für die Eckpunkte-IK
s_ep.wn = [1;0;0]; % Nehme Zielfunktion 1. Damit Überschreitung der Ränder in Zwischenständen möglich
s_ep.n_max = 5000; % Mehr Versuche (Abstände zwischen Punkten größer als bei Traj.-IK)
% s_ep.maxrelstep_ns = 0.05; % Große Werte, Größere Nullraumbewegung pro Zeitschritt
s_ep.retry_on_limitviol = true;
s_ep.retry_limit = 100; % Neuversuche erlauben (bei Einzelpunkt i.O.)
s_ep.normalize = false;
s_ep.finish_in_limits = true;
s_ep.scale_lim = 0;
QL = NaN(size(XL,1),RP.NJ);
h1 = NaN(size(XL,1),1); h2 = h1; % Zielfunktion für Aufgabenredundanz
hcond = NaN(size(XL,1),1);
t0 = tic();
for i = 1:size(XL,1)
  t1 = tic();
  % Berechne IK mit reduziertem FG
  RP.update_EE_FG(I_EE_full, I_EE_red);
  % Berechnung mit aus Vorlagendatei generierter Funktion
  [q_i, Phi_i, ~, Stats] = RP.invkin4(XL(i,:)', qs, s_ep);
  assert(all(size(Phi_i)==[RP.NLEG*sum(I_EE_full)-1 1]), ...
    'ZB Phi_i aus ParRob/invkin4 hat die falsche Dimension');
  QL(i,:) = q_i;
  % Prüfe nochmals die Richtigkeit mit anderer Modellierung
  x_i = RP.fkineEE_traj(q_i')'; % tatsächliche EE-Drehung (Freiheitsgrad der Aufg.Red.)
  [~,Phi_i_voll] = RP.constr1(q_i, x_i(:));
  assert(all(abs(Phi_i_voll)<1e-8), ...
    sprintf('Ergebnis der %s-IK (Parallel) ist falsch', I_EE_red_str));
  % Berechne die IK mit Seriell-Methode (zum Testen).
  [q_i_test, Phi_i_test] = RP.invkin_ser(XL(i,:)', qs, s_ep); % rmfield(s_ep,{'maxstep_ns','maxrelstep_ns'}));
  assert(all(size(Phi_i_test)==[RP.NLEG*sum(I_EE_full)-1 1]), ...
    'ZB Phi_i aus ParRob/invkin_ser hat die falsche Dimension');
  % Prüfe auch hiervon die Richtigkeit
  x_i_test = RP.fkineEE_traj(q_i_test');
  [~,Phi_i_voll] = RP.constr1(q_i_test, x_i_test(:));
  assert(all(abs(Phi_i_voll)<1e-9), ...
    sprintf('Ergebnis der %s-IK (Seriell) ist falsch', I_EE_red_str));

  fprintf('Eckpunkt %d/%d berechnet. Dauer %1.1fs (tpl-Funktion). Bis hier %1.1fs.\n', ...
    i, size(XL,1), toc(t1), toc(t0));
  if max(abs(Phi_i)) > 1e-6
    error('Eckpunkt %d geht nicht', i);
  end
  h1(i) = invkin_optimcrit_limits1(q_i, qlim);
  h2(i) = invkin_optimcrit_limits2(q_i, qlim);
  if any(isinf(h2(i)))
    I_limviol = find(q_i < qlim(:,1) | q_i > qlim(:,2));
    warning('Grenzverletzung bei Eckpunkt %d (kann an Einstellungen liegen)', i);
  end
  [~,Phi_q_voll] = RP.constr4grad_q(q_i);
  [~,Phi_x_voll] = RP.constr4grad_x(x_i);
  Jinv_voll = -Phi_q_voll\Phi_x_voll;
  hcond(i) = cond(Jinv_voll(RP.I_qa,RP.I_EE));
end

QL_norm = (QL-repmat(qlim(:,1)', size(QL, 1), 1)) ./ ...
          repmat(qlim(:,2)'-qlim(:,1)', size(QL, 1), 1);
fprintf(['Range of normalized joint positions in key poses:\n', ...
  'Mean: [%s]%%\n\t(Ideal: 50%%, in mid range)\n', ...
  'Range Left: [%s]%%\n\tIdeal: 50%%, highest reserve possible.\n'], ...
  disp_array(100*mean(QL_norm),'%1.0f'), ...
  disp_array(min(100*[max(1-QL_norm); -max(-QL_norm)]),'%1.0f'));
fprintf('Condition numbers for key poses: min %1.1f, max %1.1f\n', ...
  min(hcond), max(hcond));
if usr_only_test_keypoints
  return
end

% Give table of support points of the trajectory
fprintf('Table of support points:\n');
for i = 1:size(XL,1)
  fprintf('%d & %1.0f & %1.0f & %1.0f & %1.0f & %1.0f\\\\\n', i-1, ...
    1e3*XL(i,1), 1e3*XL(i,2), 1e3*XL(i,3), 180/pi*XL(i,4), 180/pi*XL(i,5));
end
%% Global Discretization of the Redundant Coordinate
% Create equidistant trajectory samples. Very fine sampling
[X_tref,~,~,~,IL_ref] = traj_trapez2_multipoint(XL, 3, 0.01, 0.01, 1e-4, 0);
% Get normalized path coordinate
s_ref = NaN(size(X_tref,1),1);
% Get progress related to key points of the trajectory
IL_ref = [IL_ref; size(X_tref,1)];
for i = 1:size(XL,1)-1
  I1 = IL_ref(i);
  I2 = IL_ref(i+1);
  p_all = ( X_tref(I1:I2,:)-repmat(XL(i,:),I2-I1+1,1) ) ./ ...
           repmat(XL(i+1,:)-XL(i,:),I2-I1+1,1);
  p = mean(p_all,2, 'omitnan');
  s_ref(I1:I2) = (i-1)+p;
end
% Remove detailed data and limit to resolution 1mm/1deg
last_x = inf(1,6);
last_s = inf;

if usr_highres_distrfig % settings for high resulution of performance map
  % three values for horizontal resolution of the performance map
  mapres_thresh_eepos = 1e-3; % 1mm
  mapres_thresh_eerot = 3*pi/180; % 1deg
  mapres_thresh_pathcoordres = 0.01; % min 100 values between two trajectory key points
  % vertical resolution of the performance map. Resolution for the
  % redundant coordinate phi_z
  mapres_redcoord_dist_deg = 0.5; % deg
else % settings for low resulution of map
  mapres_thresh_eepos = 3e-3; % 3mm
  mapres_thresh_eerot = 3*pi/180; % 3deg
  mapres_thresh_pathcoordres = 0.05;% min 20 values between key points
  mapres_redcoord_dist_deg = 5; % deg
end
for i = 1:size(X_tref,1)
  if norm(X_tref(i,1:3)-last_x(1:3)) < mapres_thresh_eepos && ...
     norm(X_tref(i,4:6)-last_x(4:6)) < mapres_thresh_eerot && ...
     abs(s_ref(i)-last_s) < mapres_thresh_pathcoordres
    X_tref(i,:) = NaN; % diesen Punkt löschen
  else
    last_x = X_tref(i,:);
    last_s = s_ref(i);
  end
end
I_remove = isnan(X_tref);
X_tref = X_tref(any(~I_remove,2),:);
s_ref = s_ref(any(~I_remove,2));
% Debug: For fast check of validity
% X_tref = X_tref(1:5:30,:); % Debug
% s_ref = s_ref(1:5:30,:); % Debug
% Create range of values for the redundant coordinate. Only consider one
% sign here and consider negative signs below
phiz_range_calc = (0:mapres_redcoord_dist_deg:180)*pi/180;
phiz_range = NaN(1,2*length(phiz_range_calc)-1); % for both signs
filename_pre = sprintf('pkm_traj_%dtraj_%dphiz_discretization_%s', ...
  size(X_tref,1), length(phiz_range_calc), usr_discretization_type);
filename_discr = fullfile(respath, [filename_pre, '_data.mat']);
data_loaded_offline = false;
if usr_load_discretization && ~exist(filename_discr, 'file')
  warning('Unable to load discretization results from %s', filename_discr);
elseif usr_load_discretization && exist(filename_discr, 'file')
  d = load(filename_discr);
  if all(size(XL)==size(d.XL))
    test_XL = XL - d.XL;
  else
    test_XL = 1; % raises error
  end
  if all(size(phiz_range)==size(d.phiz_range))
    test_phiz = phiz_range - d.phiz_range;
  else
    test_phiz = 1;
  end
  if any(abs(test_XL(:)) > 1e-8) || any(abs(test_phiz(:)) > 1e-8)
    warning('Data from file does not match the settings');
  elseif length(intersect(fields(d), {'H_all', 'XL', 'phiz_range', 'Q_all'}))~=4
    warning('Data from file misses fields');
  else
    data_loaded_offline = true;
    % Load results
    H_all = d.H_all;
    Q_all = d.Q_all;
    phiz_range = d.phiz_range;
  end
end
fprintf(['The performance map contains %d trajectory samples and %d ', ...
  'values for the redundant coordinates. %d evaluations in total\n'], ...
  size(X_tref,1), 2*length(phiz_range_calc), size(X_tref,1)*2*length(phiz_range_calc));

if ~data_loaded_offline
fprintf('Start discretization of the inverse kinematics\n');
H_all = NaN(size(X_tref,1), 2*length(phiz_range_calc)-1, 4);
Q_all = NaN(2*length(phiz_range_calc)-1, RP.NJ, size(X_tref,1));
% Einstellungen für Dummy-Berechnung ohne Änderung der Gelenkwinkel.
s_ep_dummy = s_ep;
s_ep_dummy.retry_limit = 0;
s_ep_dummy.wn = ones(4,1); % hierdurch werden die Kriterien berechnet
s_ep_dummy.K = zeros(RP.NJ,1); % hierdurch keine Bewegung und damit ...
s_ep_dummy.Kn = zeros(RP.NJ,1); % ... sofortiger Abbruch
s_ep_dummy.optimcrit_limits_hyp_deact = optimcrit_limits_hyp_deact;
% Einstellung für IK bei globaler Diskretisierung über Trajektorie
s_ep_glbdscr = s_ep; % leave as is (temporarily permit limit violations; no scaling until limit)
s_ep_glbdscr.retry_limit = 10;
s_ep_glbdscr.normalize = false; % no normalization (due to joint limits)
s_ep_glbdscr = rmfield(s_ep_glbdscr, 'finish_in_limits'); % does not work without redundancy
s_traj_glbdscr = struct('simplify_acc', true);

for ii_sign = 0:1 % move redundant coordinate in positive and negative direction
  fprintf('Start computation for sign %+d\n', (-1)^ii_sign);
  H_all_ii = NaN(size(X_tref,1), length(phiz_range_calc), 4);
  Q_all_ii = NaN(length(phiz_range_calc), RP.NJ, size(X_tref,1));
  t_lastmessage = tic();
  t1 = tic();
  phiz_range_ii = (-1)^ii_sign * phiz_range_calc;
  for i = 1:size(X_tref,1) % loop trajectory samples
    t1_i = tic();
    % Get joint configurations using a virtual trajectory
    if i > 1
      q0_j = Q_all_ii(1, :, i-1)';
    else
      q0_j = q0;
    end
    if strcmp(usr_discretization_type, 'trajectory')
      RP.update_EE_FG(I_EE_full,I_EE_full); % Roboter auf 3T3R einstellen
      X_i_traj = [repmat(X_tref(i,1:5), length(phiz_range_ii), 1), phiz_range_ii'];
      XD_i_traj = [zeros(length(phiz_range_ii), 5), [0;diff(phiz_range_ii(:))]];
      XDD_i_traj = zeros(length(phiz_range_ii), 6);
      T_i_traj = (0:1:length(phiz_range_ii)-1)';
      [Q_i,~,~,~,~,~,~,Stats] = RP.invkin2_traj(X_i_traj, XD_i_traj, XDD_i_traj, T_i_traj, q0_j, s_traj_glbdscr); 
      RP.update_EE_FG(I_EE_full,I_EE_red); % Roboter auf 3T2R einstellen
      for j = 1:length(phiz_range_ii)
        x_j = [X_tref(i,1:5), phiz_range_ii(j)]';
        q_j = Q_i(j,:)';
        % IK benutzen, um Zielfunktionswerte zu bestimmen (ohne Neuberechnung)
        [q_dummy, Phi,~,Stats_dummy] = RP.invkin4(x_j, q_j, s_ep_dummy);
        if any(abs(q_j - q_dummy) > 1e-8)
          error('IK-Ergebnis hat sich bei Test verändert');
        end
        Q_all_ii(j, :, i) = q_j;
        H_all_ii(i,j,:) = Stats_dummy.h(Stats_dummy.iter+1,2:end);
      end
    elseif strcmp(usr_discretization_type, 'position-level')
    % Get joint configurations using position-level IK (less efficient)
    % Here enforcing the joint limits is easier, but not necessary to
    % determine the condition numbers
      for j = 1:length(phiz_range_ii)
        if j == 1
          % assure joint limits by retrying, only they match
          s_ep_glbdscr.retry_on_limitviol = true;
        else
          % Not able to retry on position limits violation. Otherwise,
          % platform angles greater than 180° can not be checked. For this
          % the initial value has to be the next point
          s_ep_glbdscr.retry_on_limitviol = false;
        end
        RP.update_EE_FG(I_EE_full,I_EE_full); % Roboter auf 3T3R einstellen
        % IK benutzen, um Zielfunktionswerte zu bestimmen (ohne Neuberechnung)
        t1_j = tic();
        x_j = [X_tref(i,1:5), phiz_range_ii(j)]';
        if j > 1
          delta_x = [zeros(5,1);phiz_range_ii(j)-phiz_range_ii(j-1)];
          [~,Phi_q] = RP.constr4grad_q(Q_all_ii(j-1, :, i)');
          [~,Phi_x] = RP.constr4grad_x([X_tref(i,1:5),phiz_range_ii(j-1)]');
          Jtilde_inv_x = -Phi_q\Phi_x; % Full coordinate Jacobian (equ. 17 in paper)
          q0_j = Q_all_ii(j-1, :, i)' + Jtilde_inv_x*delta_x;
        else
          if i > 1
            q0_j = []; % will be selected below
          else
            q0_j = q0;
          end
        end
        if any(isnan(q0_j))
          % The previous pose was not computed successfully. Take the
          % second last and so on
          q0_j = q0; % overwrite this later
          for jjj = j-1:-1:1 % look back until the first platform rotation
            if ~any(isnan(Q_all_ii(jjj, :, i)))
              q0_j = Q_all_ii(jjj, :, i)';
              break;
            end
          end
        end
        % create list of initial values for the IK and compute IK
        q0_list = q0_j';
        if j > 1 && ~any(isnan(Q_all_ii(j-1, :, i)))
          q0_list = [q0_list; Q_all_ii(j-1, :, i)];  %#ok<AGROW>
        end
        if i > 1 && j > 1 && ~any(isnan(Q_all_ii(j-1, :, i-1)))
          q0_list = [q0_list; Q_all_ii(j-1, :, i-1)];  %#ok<AGROW>
        end
        if i > 1 && ~any(isnan(Q_all_ii(j, :, i-1)))
          q0_list = [q0_list; Q_all_ii(j, :, i-1)];  %#ok<AGROW>
        end
        if i == 1 && j == 1
          % Add random poses to be able to start with the best pose
          q0_list = [q0_list; repmat(qlim(:,1)',200,1)+rand(200, RP.NJ).* ...
            repmat(qlim(:,2)' - qlim(:,1)',200,1)];  %#ok<AGROW>
        end
        % In case of second run overwrite everything and directly take
        % results for phi=0. Otherwise their may be a discontinuity
        if ii_sign == 1 && j == 1 % first value corresponds to 0
          q0_list = Q_all(length(phiz_range_ii),:,i);
        end
        if any(isnan(q0_list(:))), error('An Initial value is NaN'); end % this makes a random new initial seed, which is not desired here
        q_j_list = NaN(size(q0_list));
        for k = 1:size(q0_list,1)
          [q_k, Phi, ~, Stats] = RP.invkin2(x_j, q0_list(k,:)', s_ep_glbdscr);
          if any(abs(Phi) > 1e-8)
            % Try other IK method
            [q_k, Phi, ~, Stats] = RP.invkin4(x_j, q0_list(k,:)', s_ep_glbdscr);
            % Mark this as not working, if other method fails as well
            if any(abs(Phi) > 1e-8)
              q_j_list(k,:) = NaN;
              continue
            end
          end
          % normalize joint angles (to stay near the limits' center)
          % Do not do this to see what happens after a full rotation
          % q_k(RP.MDH.sigma==0) = normalizeAngle(q_k(RP.MDH.sigma==0), ...
          %   mean(qlim(RP.MDH.sigma==0,:),2)); % normalize to center of limits
          q_j_list(k,:) = q_k;
        end
        if all(isnan(q_j_list(:)))
          warning('IK did not find a solution for phi_z=%1.1fdeg', 180/pi*x_j(6));
          continue % this IK configuration did not work.
        end
        % Compute performance criteria
        RP.update_EE_FG(I_EE_full,I_EE_red); % Roboter auf 3T2R einstellen
        h_list = NaN(size(q_j_list,1),4);
        q_dist = NaN(size(q_j_list,1),4);
        for k = 1:size(q_j_list,1)
          % IK benutzen, um Zielfunktionswerte zu bestimmen (ohne Neuberechnung)
          [q_dummy, Phi,~,Stats_dummy] = RP.invkin4(x_j, q_j_list(k,:)', s_ep_dummy);
          if any(abs(q_j_list(k,:)' - q_dummy) > 1e-8)
            error('IK-Ergebnis hat sich bei Test verändert');
          end
          h_list(k,:) = Stats_dummy.h(Stats_dummy.iter+1,2:end);
        end
        % select the joint angles that are nearest to the previous pose
        if ~(i==1 && j == 1) % not possible for first sample
          [~,Ibest_k] = min(sum((q_j_list-repmat(q0_list(1,:),size(q_j_list,1),1)).^2,2));
        else
          % Alternative: Pick the best joint configuration (away from limits)
          % This may require reconfigurations within the plane phiz-s.
          % Therefore only use for the first sample.
          [~, Ibest_k] = min(h_list(:,2)+h_list(:,1)); % squared and hyperbolic limit. Use best
        end
        if any(isnan(q_j_list(Ibest_k,:))), error('Unexpected NaN'); end
        Q_all_ii(j, :, i) = q_j_list(Ibest_k,:);
        H_all_ii(i,j,:) = h_list(Ibest_k,:);
      end
    else
      error('Mode not defined');
    end
    if toc(t_lastmessage)-toc(t1_i) > 20 || i == 1
      fprintf(['Duration for trajectory sample number %d/%d: %1.1fs. Remaining %d ', ...
        'samples (estimated %1.1fmin)\n'], i, size(X_tref,1), toc(t1_i), ...
        (size(X_tref,1)-i), toc(t1_i)*(size(X_tref,1)-i)/60);
      t_lastmessage = tic();
    end
  end
  % Save data for the sign value ii
  if ii_sign == 0 % positive sign
    % take all the result data
    I_phi = 1:length(phiz_range_ii);
    % store in second half of result variables
    I_all = length(phiz_range_ii):(2*length(phiz_range_ii)-1);
  else
    % Beginning of the interval. The computation starts with 0, but the
    % lowest value has to be saved as the first one. Therefore `flipör`.
    % Remove the 0 entry because it is already included in the positive
    % values.
    I_phi = fliplr(2:length(phiz_range_ii));
    I_all = 1:length(phiz_range_ii)-1;
  end
  phiz_range(I_all) = phiz_range_ii(I_phi);
  for k = 1:size(Q_all_ii,3)
    Q_all(I_all,:,k) = Q_all_ii(I_phi,:,k);
  end
  for k = 1:size(H_all_ii,3)
    H_all(:,I_all,k) = H_all_ii(:,I_phi,k);
  end
end
fprintf('Finished discretization of the trajectory. Total %1.1fmin)\n', toc(t1)/60);
save(filename_discr, 'H_all', 'Q_all', 'phiz_range', 'XL');
end
if length(unique(phiz_range))~=length(phiz_range)
  error('Something went wrong when assembling phiz_range');
end
%% Plot global distribution
if debug_plot % Distribution of PKM Jacobian condition number
  change_current_figure(3);clf; set(3, 'Name', 'Distr_PKMJacCrit', 'NumberTitle', 'off');
  [X,Y] = meshgrid(s_ref,180/pi*phiz_range);
  Z = zeros(size(X));
  % Create color code from Jacobian condition number
  CC = H_all(:,:,4)';
  % high condition numers all get the same dark color to be able to
  % distinguish the good values better.
  CC(CC>1e6) = 1e6; % limit to control colors.
  condsat_limit = 2e2;
  I_exc = CC > condsat_limit;
  CC(I_exc) = condsat_limit+10*log10(CC(I_exc)/condsat_limit); % last term gives 0...40 for condsat_limit=1e2
  % Remove values with invalid joint configurations
  I_limviol = isinf(H_all(:,:,2)');
  % CC(I_limviol) = NaN;
  % Create color plot
  surf(X,Y,Z,CC, 'EdgeColor', 'none');
  xlabel('Trajectory coordinate (norm. per point)');
  ylabel('redundant coordinate in deg');
  title(sprintf('Jacobian Condition Number Criterion (%1.0f=singular)', condsat_limit));
  colorbar
  view([0,90])
  % Set Colormap: low condition numbers white, high/singularity dark red.
  colormap(flipud(hot(1024)));
  if usr_save_figures
    saveas(3, fullfile(respath, [filename_pre,'_pkmjac_condition.fig']));
    saveas(3, fullfile(respath, [filename_pre,'_pkmjac_condition.png']));
  end
end
if debug_plot % Distribution of inverse kinematics Jacobian condition number
  change_current_figure(6);clf; set(6, 'Name', 'Distr_IKJacCrit', 'NumberTitle', 'off');
  [X,Y] = meshgrid(s_ref,180/pi*phiz_range);
  Z = zeros(size(X));
  % Create color code from inverse kinematics Jacobian condition number
  CC = H_all(:,:,3)';
  % high condition numers all get the same dark color to be able to
  % distinguish the good values better.
  CC(CC>1e6) = 1e6; % limit to control colors.
  condsat_limit = 2e2;
  I_exc = CC > condsat_limit;
  CC(I_exc) = condsat_limit+10*log10(CC(I_exc)/condsat_limit); % last term gives 0...40 for condsat_limit=1e2
  % Remove values with invalid joint configurations
  I_limviol = isinf(H_all(:,:,2)');
  % CC(I_limviol) = NaN;
  % Create color plot
  surf(X,Y,Z,CC, 'EdgeColor', 'none');
  xlabel('Trajectory coordinate (norm. per point)');
  ylabel('redundant coordinate in deg');
  title(sprintf('IK Jacobian Condition Number Criterion (%1.0f=singular)', condsat_limit));
  colorbar
  view([0,90])
  % Set Colormap: low condition numbers white, high/singularity dark red.
  colormap(flipud(hot(1024)));
  if usr_save_figures
    saveas(6, fullfile(respath, [filename_pre,'_ikjac_condition.fig']));
    saveas(6, fullfile(respath, [filename_pre,'_ikjac_condition.png']));
  end
end
if debug_plot % Distribution of joint limits criterion
  change_current_figure(5);clf; set(5, 'Name', 'Distr_JointLimCrit', 'NumberTitle', 'off');
  % Create color code from joint limits violation criterion
  CC = H_all(:,:,2)';
%   CC(isinf(CC)) = NaN;
  CC(CC>1e6) = 1e6; % limit to control colors.
  limcrit_satlimit = 1e3;
  I_exc = CC > limcrit_satlimit;
  CC(I_exc) = limcrit_satlimit+10*log10(CC(I_exc)/limcrit_satlimit); % last term gives 0...40 for condsat_limit=1e2
  surf(X,Y,Z,CC, 'EdgeColor', 'none');
  xlabel('Trajectory coordinate (norm. per point)');
  ylabel('redundant coordinate in deg');
  title('Joint Limit Criterion (1000=invalid)');
  colorbar
  view([0,90])
  % Set Colormap:
  colormap(flipud(winter(1024)));
  if usr_save_figures
    saveas(5, fullfile(respath, [filename_pre,'_limitviol_crit.fig']));
    saveas(5, fullfile(respath, [filename_pre,'_limitviol_crit.png']));
  end
end
if debug_plot % Distribution of combined criterion (condition+limits)
  change_current_figure(8);clf; set(8, 'Name', 'Distr_CombCrit', 'NumberTitle', 'off');
  % Create color code from combined criterion
  CC = H_all(:,:,2)' + H_all(:,:,4)';
  CC(isinf(CC)) = 1e6;
  CC(CC>1e6) = 1e6; % limit to control colors.
  critsat_limit = 2e2;
  I_exc = CC > critsat_limit;
  CC(I_exc) = critsat_limit+10*log10(CC(I_exc)/critsat_limit);
  surf(X,Y,Z,CC, 'EdgeColor', 'none');
  xlabel('Trajectory coordinate (norm. per point)');
  ylabel('redundant coordinate in deg');
  title('Combined Criterion (>200=invalid/singular)');
  colorbar
  view([0,90])
  % Set Colormap:
  colormap(flipud(hot(1024)));
  if usr_save_figures
    saveas(8, fullfile(respath, [filename_pre,'_combined_crit.fig']));
    saveas(8, fullfile(respath, [filename_pre,'_combined_crit.png']));
  end
  change_current_figure(9);clf; set(9, 'Name', 'Distr_CombCrit2', 'NumberTitle', 'off');
  % Create color code from combined criterion
  CC = H_all(:,:,2)' + H_all(:,:,3)' + H_all(:,:,4)';
  CC(isinf(CC)) = 1e6;
  CC(CC>1e6) = 1e6; % limit to control colors.
  critsat_limit = 2e2;
  I_exc = CC > critsat_limit;
  CC(I_exc) = critsat_limit+10*log10(CC(I_exc)/critsat_limit);
  surf(X,Y,Z,CC, 'EdgeColor', 'none');
  xlabel('Trajectory coordinate (norm. per point)');
  ylabel('redundant coordinate in deg');
  title('Combined Criterion 2 (>200=invalid/singular)');
  colorbar
  view([0,90])
  % Set Colormap:
  colormap(flipud(hot(1024)));
  if usr_save_figures
    saveas(9, fullfile(respath, [filename_pre,'_combined_crit2.fig']));
    saveas(9, fullfile(respath, [filename_pre,'_combined_crit2.png']));
  end
end

if debug_plot
  change_current_figure(4); clf;
  for i = 1:6
    if i < 4
      scale = 1e3;
      ylt = sprintf('%s in mm', char(119+i));
    else
      scale = 180/pi;
      ylt = sprintf('phi %s in deg', char(119+i-3));
    end
    subplot(2,3,i);
    plot(s_ref, scale*X_tref(:,i));
    xlabel('Trajectory coordinate (norm. per point)');
    grid on;
    ylabel(ylt);
  end
end

if false && debug_plot % For Debugging the global discretization
  i_traj = 1;
  [~,i_phi] = min(abs(phiz_range-0));
  % Plot joint angles for one platform rotation
  change_current_figure(901);clf;
  ii = 0;
  for i = 1:RP.NLEG
    for j = 1:RP.Leg(i).NJ
      ii = ii + 1;
      subplot(RP.NLEG,RP.Leg(1).NJ,sprc2no(RP.NLEG, RP.Leg(1).NJ, i, j));
      hold on; grid on;
      stairs(180/pi*phiz_range, Q_all(:,ii,i_traj));
      plot(180/pi*phiz_range([1,end]), repmat(RP.Leg(i).qlim(j,:),2,1), 'r-');
      if j == 1, ylabel(sprintf('BK %d', i)); end
      if i == 1, title(sprintf('q %d', j)); end
    end
  end
  sgtitle(sprintf('Joint configurations for phi_z at s=%1.2f (sample %d)', ...
    s_ref(i_traj), i_traj));
  linkxaxes
end
%% Create paper figure with robot in initial pose of the trajectory
% Plot arbitrary configuration to check if the platform rotation worked
i_traj = 1; % beginning of the trajectory
[~,i_phi] = min(abs(phiz_range-0)); % value for phiz=0
q_test = Q_all(i_phi,:,i_traj)';
x_test = [X_tref(i_traj,1:5), phiz_range(i_phi)]';
% Test the loaded configuration
[~,Phi_test] = RP.constr1(q_test, x_test);
assert(all(abs(Phi_test)<1e-10), 'Kinematic constraints in loaded config. do not match');
change_current_figure(900);clf;
set(900,'Name','Rob','NumberTitle','off');
% title(sprintf('Robot in Selected Pose: traj idx %d, phiz=%1.0fdeg', ...
%   i_traj, 180/pi*x_test(6)));
hold on;grid on;
view(3);
s_plot = struct(  'ks_legs', [], 'straight', 0, 'mode', 4);
% Modify the robot pose to avoid overlapping of the platform and the
% trajectory. Only minor change for visibility.
x_paperfig = x_test;
x_paperfig(1:3) = x_paperfig(1:3) + eulxyz2r(x_paperfig(4:6))*[0;0;-20e-3];
[q_paperfig, Phi_paperfig] = RP.invkin_ser(x_paperfig, q_test);
assert(all(abs(Phi_paperfig)<1e-6), 'IK for robot paper figure not successful');
RP.plot( q_paperfig, x_paperfig, s_plot );
Traj_0 = struct('X', XL);
Traj_W = RP.transform_traj(Traj_0);
plot3(Traj_W.X(:,1), Traj_W.X(:,2), Traj_W.X(:,3), 'k-');
Xii_lastplot = inf(1,6);
for ii = 1:size(X_t,1)
  T_ii = RP.x2t(X_t(ii,:)');
  % Reduce density of the coordinate systems in the plot
  if ii > 1 && norm(X_t(ii,1:3)-Xii_lastplot(1:3))<0.02, continue; end
  Xii_lastplot = X_t(ii,:);
  trplot(T_W_0*T_ii, 'frame', '', 'rgb', 'length', 0.02, 'labels', '   ')
end
% xlabel('x in m');ylabel('y in m');zlabel('z in m');
axis equal
% Format plot
figure_format_publication(gca);
set(gca, 'Box', 'off');
set(gca,'XTICKLABEL',{});set(gca,'YTICKLABEL', {});set(gca,'ZTICKLABEL',{});
set(gca,'xtick',[],'ytick',[],'ztick',[]);
set(get(gca, 'XAxis'), 'visible', 'off');
set(get(gca, 'YAxis'), 'visible', 'off');
set(get(gca, 'ZAxis'), 'visible', 'off');
ch = get(gca, 'Children');
for jj = 1:length(ch)
  % KS-Texte entfernen
  if strcmp(ch(jj).Type, 'hgtransform')
    chjjch = ch(jj).Children;
    if strcmp(chjjch(1).Type, 'text') && ~strcmp(chjjch(1).String, ' ')
      % this handle is a coordinate system with text annotation (i.e. 0/E)
      % Delete this frame (not necessary in this figure)
      delete(ch(jj));
      continue;
    end
%     % Entferne den KS-Text wieder
%     delete(chjjch(1:4));
  end
end
% Remove all transparency. Produces problems in png image export
ch = get(gca, 'Children');
for jj = 1:length(ch)
  if strcmp(ch(jj).Type, 'surface') || strcmp(ch(jj).Type, 'patch')
    % Remove transparency (set everything to non-transparent)
    set(ch(jj), 'FaceAlpha', 1.0, 'EdgeAlpha', 1.0);
    % set edge color to face color (for robot link elements)
    if ~isempty(regexp(get(ch(jj), 'DisplayName'), '^Leg_([\d]+)_Link_([\d]+)$', 'match'))
      set(ch(jj), 'FaceColor', get(ch(jj), 'EdgeColor'));
    end
    if ~isempty(regexp(get(ch(jj), 'DisplayName'), '^Leg_([\d]+)_Joint_([\d]+)$', 'match'))
%       set(ch(jj), 'EdgeColor', get(ch(jj), 'FaceColor'));
%       chch = get(ch(jj), 'children')
%       for k = 1:length(chch)
%         set(chch(k), 'EdgeColor', get(chch(k), 'FaceColor'));
%       end
    end
    % remove the cuboid representation of prismatic joints
    if ~isempty(regexp(get(ch(jj), 'DisplayName'), '^Leg_([\d]+)_Joint_3$', 'match'))
      delete(ch(jj));
    end
  end
end
% Set different camera angle to avoid obfuscation of the trajectory
view([140, 20]);

set(900, 'windowstyle', 'normal');
set_size_plot_subplot(900, ...
  10,10,gca,...
  0,0,0,0,0,0)
drawnow();
% Bild speichern
if usr_save_figures
  filename_savepre = sprintf('pkm_traj_start_pose_phi%1.0f', 180/pi*phiz_range(i_phi));
  exportgraphics(gcf, fullfile(paperfig_path, [filename_savepre,'.png']),'Resolution',600)
  saveas(900, fullfile(respath, [filename_savepre,'.fig']));
end
if usr_only_global_discretization
  return
end

%% Generate Trajectory
% Use moderate end effector velocity and acceleration to have reserves for
% nullspace motion
vmax = 50e-3; % 50mm/s
phiDmax = 10*pi/180; % 10°/s

[X_t, XD_t, XDD_t, t, IL] = traj_trapez2_multipoint(XL, [vmax*[1;1;1];phiDmax*[1;1;1]]', ...
  0.01, 0.01, 1e-3, 0);
if usr_short_traj
  n = 2000;
else
  n = length(t);
end
IL = IL(IL<length(t));
II = 1:n;
t = t(II);
X_t = X_t(II,:);
XD_t = XD_t(II,:);
XDD_t = XDD_t(II,:);
fprintf('Length of the trajectory: %1.1fmm, %1.1fdeg\n', ...
  1e3*sum(sum(abs(diff(X_t(:,1:3))))), 180/pi*sum(sum(abs(diff(X_t(:,4:5))))))
fprintf('Duration of the trajectory: %1.1fs\n', t(end));
%% Inverse Kinematik zum Startpunkt der Trajektorie
% Inverse Kinematik berechnen;  Lösung der IK von oben als Startwert
% Wird nicht für die Trajektorien-IK benutzt, da die optimale Startkon-
% figuration von den benutzten Nebenbedingungen abhängt.
tic();
s_start = s;
% Toleranz maximal stark setzen, damit es keinen Sprung im Verlauf gibt
% (durch die vielen Nullraumiterationen ist die Aufgabentoleranz später
% sowieso ungefähr Null.
s_start.Phit_tol = 1e-12;
s_start.Phir_tol = 1e-12;
s_start.normalize = false;
s_start.finish_in_limits = true;
s_start.scale_lim = 0;
s_start.maxstep_ns = 1e-10; % Nullraumbewegung sollte schon zum Optimum konvergiert sein
% Berechne IK mit 3T2R
RP.update_EE_FG(I_EE_full, I_EE_red);
warning on
% Berechne Ersten Punkt der Trajektorie mit Aufgabenredundanz.
% Dadurch bestmögliche Startkonfiguration
[q1, Psi_num1, ~, Stats1] = RP.invkin3(X_tref(1,:)', qs, s_start);
x1 = RP.fkineEE_traj(q1')';
if any(abs(Psi_num1) > 1e-4)
  error('IK konvergiert nicht für Startpunkt der Trajektorie');
end
% % Normiere die Start-Gelenkwinkel auf Bereich 0 bis 1
qsnorm = (qs-qlim(:,1)) ./ (qlim(:,2) - qlim(:,1)); % Muss 0.5 sein per Definition
if any(abs(qsnorm(~RP.I_qa)-0.5)>1e-6)
  error('Die erste Startkonfiguration für die IK ist nicht in der Mitte des Gelenkraums');
end
% Prüfe den Anfangswert für die IK (Optimal im Startpunkt)
q1norm = (q1-qlim(:,1)) ./ (qlim(:,2) - qlim(:,1));
if any(q1norm > 1) || any(q1norm < 0) % Winkel mit +- 2*pi berücksichtigen
  error('Anfangs-Konfiguration für Trajektorie verletzt bereits die Grenzen');
end
fprintf('Best platform orientation in starting pose of trajectory: %1.0fdeg\n', 180/pi*x1(6));

%% IK für Trajektorie berechnen (Vorbereitung)
RP.update_EE_FG(I_EE_full, I_EE_red);
% Berechne Beispiel-Trajektorie für 3T2R
% Einstellungen für Trajektorie: Kein Normalisieren, damit Traj. nicht
% springt. Muss doch normalisieren, damit Gelenkwinkelgrenzen korrekt
% erkannt werden.
s_Traj = struct('Phit_tol', 1e-12, 'Phir_tol', 1e-12);
% Abspeichern der Gelenkwinkel für verschiedene Varianten der Berechnung
Namen_Methoden = cell(1,5);
Namen_Methoden_Leg_Paper = Namen_Methoden;
Q_t_all = NaN(length(t), RP.NJ, length(Namen_Methoden)); QD_t_all = Q_t_all; QDD_t_all = Q_t_all;
Q_t_norm_all = Q_t_all;
H1_all = NaN(length(t), 1+RP.NLEG, length(Namen_Methoden)); H2_all = H1_all;
H1D_all = H1_all; H2D_all = H1_all;
Hcond_all = NaN(length(t), 2, length(Namen_Methoden));
XE_all = NaN(length(t), 6, length(Namen_Methoden));
XDE_all = XE_all;

%% IK für Trajektorie berechnen (verschiedene Methoden)
qlim_backup = cat(1, RP.Leg.qlim);
qDlim_backup = cat(1, RP.Leg.qDlim); % damit überschriebene Werte wieder ...
qDDlim_backup = cat(1, RP.Leg.qDDlim); % ... hergestellt werden können

filename_traj = fullfile(respath, [filename_pre, '_traj.mat']);
if usr_load_traj && ~exist(filename_traj, 'file')
  warning('Unable to load trajectory results from %s', filename_traj);
  usr_load_traj = false;
end
if usr_load_traj
  d = load(filename_traj);
  if size(d.Q_t_all,1) ~= length(t)
    warning('Saved trajectory has wrong length');
    usr_load_traj = false;
  end
end
if usr_load_traj
  Q_t_all = d.Q_t_all;
  QD_t_all = d.QD_t_all;
  QDD_t_all = d.QDD_t_all;
  XE_all = d.XE_all;
  XDE_all = d.XDE_all;
  Q_t_norm_all = d.Q_t_norm_all;
  H_all = d.H_all;
  H1_all = d.H1_all;
  H1D_all = d.H1D_all;
  H2_all = d.H2_all;
  H2D_all = d.H2D_all;
  Hcond_all = d.Hcond_all;
  Namen_Methoden = d.Namen_Methoden;
  Namen_Methoden_Leg_Paper = d.Namen_Methoden_Leg_Paper;
  fprintf('Loaded results for trajectory without calculation\n');
  if usr_test_timing
    warning('Timinig evaluation can not be performed');
  end
end
for kk = 1:length(Namen_Methoden)*(~usr_load_traj)
  set_kk = s_Traj;
  set_kk.debug = true; % check for errors in trajectory IK function
  wn_traj_default = zeros(10,1);
  wn_traj_default(2) = 1; % K_P (hyperb. limit)
  wn_traj_default(3) = 0.03; % K_v
  wn_traj_default(6) = 0.05; % K_P (cond)
  wn_traj_default(8) = 0.6; % K_D (limit)
  wn_traj_default(10) = 0.01; % K_D (cond)
  set_kk.thresh_ns_qa = 1; % always use full joint projector
  for j = 1:RP.NLEG
    RP.Leg(j).qlim = qlim_backup(RP.I1J_LEG(j):RP.I2J_LEG(j),:);
    RP.Leg(j).qDlim = qDlim_backup(RP.I1J_LEG(j):RP.I2J_LEG(j),:);
    RP.Leg(j).qDDlim = qDDlim_backup(RP.I1J_LEG(j):RP.I2J_LEG(j),:);
  end
  phiz_0 = NaN; %#ok<NASGU>
  phiz_const = 0; % is ignored in 3T2R IK
  wn_traj = wn_traj_default;
  switch kk
    case 1
      name_method=sprintf('%s-IK', I_EE_full_str);
      name_method_leg = 'reference';
      I_EE_Task_kk = I_EE_full;
      set_kk.wn(:) = 0; % Keine zusätzliche Optimierung möglich
      phiz_0 = 30*pi/180; % take approximate average best value from global map
      % phiz_0 = X_t(1,6); % take initial value from reference traj.
      phiz_const = phiz_0; % leave platform rotation constant
    case 2
      name_method='phi-30';
      name_method_leg = '$\varphi_{z,0} = -30^\circ$';
      I_EE_Task_kk = I_EE_red;
      phiz_0 = -30*pi/180;
    case 3
      name_method='phi45';
      name_method_leg = '$\varphi_{z,0} = 45^\circ$';
      I_EE_Task_kk = I_EE_red;
%       for j = 1:RP.NLEG
%         RP.Leg(j).qDlim = 3*qDlim_backup(RP.I1J_LEG(j):RP.I2J_LEG(j),:);
%         RP.Leg(j).qDDlim = 5*qDDlim_backup(RP.I1J_LEG(j):RP.I2J_LEG(j),:);
%       end
      phiz_0 = 45*pi/180;
      % Select lower gains to avoid overshoot
      wn_traj = wn_traj_default;
      wn_traj(3) = 0.1; % K_v
      wn_traj(6) = wn_traj(6)/5; % K_P (cond)
      wn_traj(10) = wn_traj(10)/5; % K_D (cond)
    case 4
      name_method='phi90';
      name_method_leg = '$\varphi_{z,0} = 90^\circ$';
      I_EE_Task_kk = I_EE_red;
      phiz_0 = 90*pi/180;
      % Select lower gains to avoid overshoot
      wn_traj = wn_traj_default;
      wn_traj(3) = 0.25; % K_v
      wn_traj(6) = wn_traj(6)/10; % K_P (cond)
      wn_traj(10) = wn_traj(10)/10; % K_D (cond)
    case 5
      name_method='phi45_lowgain';
      name_method_leg = '$\varphi_{z,0} = 45^\circ$ (low gains)';
      I_EE_Task_kk = I_EE_red;
      phiz_0 = 45*pi/180;
      % Select lower gains to avoid overshoot
      wn_traj = wn_traj_default;
      wn_traj(3) = 0.5; % K_v
      wn_traj(6) = wn_traj(6)/20; % K_P (cond)
      wn_traj(10) = wn_traj(10)/20; % K_D (cond)
%     case 3
%       name_method='noposlimit';
%       I_EE_Task_kk = I_EE_red;
%       for j = 1:RP.NLEG
%         RP.Leg(j).qlim = repmat([-inf, inf], RP.Leg(j).NJ, 1);
%         RP.Leg(j).qDlim = 3*qDlim_backup(RP.I1J_LEG(j):RP.I2J_LEG(j),:);
%         RP.Leg(j).qDDlim = 5*qDDlim_backup(RP.I1J_LEG(j):RP.I2J_LEG(j),:);
%       end
    otherwise
      error('Fall %d noch nicht definiert', kk);
  end
  if isnan(phiz_0)
    error('phiz_0 is not set');
  end
  set_kk.wn = wn_traj;
  % Positions-IK zum Startpunkt der Trajektorie mit genau den gleichen
  % Optimierungs-Nebenbedingungen wie in der Trajektorie. Dadurch keine
  % Nullraumbewegung am Anfang (Start in lokalem Optimum)
  s_pik_kk = struct('scale_lim',0.5); % Grenzen dürfen nicht überschritten werden
  % s_pik_kk.wn = s_kk.wn([1 2 5]); % Positions-Grenzen und Kondition
  % Wähle immer die gleichen Nebenbedingungen, damit alle mit gleicher
  % Konfiguration starten (besser vergleichbar)
  RP.update_EE_FG(I_EE_full, I_EE_full);
  X_t_kk = X_t; XD_t_kk = XD_t; XDD_t_kk = XDD_t;
  X_t_kk(:,6) = phiz_const; XD_t_kk(:,6) = 0; XDD_t_kk(:,6) = 0;
  [qs_kk, Phi_s, ~, Stats_s] = RP.invkin2([X_t(1,1:5),phiz_0]', q1, ...
    struct('retry_on_limitviol', true, 'Phit_tol', 1e-12, 'Phir_tol', 1e-12));
  if any(abs(Phi_s)>1e-6)
    error('IK des Startpunkts fehlgeschlagen');
  end
  % Normalize angle. For a wide range of joint limits the IK solution can
  % be far from the center which later leads to oscillations in the traj-
  % ectory.
  qs_kk(RP.MDH.sigma==0) = normalizeAngle(qs_kk(RP.MDH.sigma==0), mean(qlim(RP.MDH.sigma==0,:),2));
  assert(all(abs(RP.constr3(qs_kk,[X_t(1,1:5),phiz_0]'))<1e-10), 'angle normalization failed');
%   [qs_kk, Phi_s, ~, Stats_s] = RP.invkin3(X_t(1,:)', qs, s_pik_kk);
  qskknorm = (qs_kk-qlim(:,1)) ./ (qlim(:,2) - qlim(:,1));
  if any(qs_kk<qlim(:,1) | qs_kk>qlim(:,2))
    I_viol = find(qs_kk<qlim(:,1) | qs_kk>qlim(:,2));
    error('Startpunkt liegt außerhalb der Grenzen. Verletzt: [%s]', disp_array(I_viol','%d'));
  end
  xs_kk_test = RP.fkineEE_traj(qs_kk')';
  assert(abs(xs_kk_test(6)-phiz_0)< 1e-10, 'redundant coordinate in starting pose is not as requested');
  RP.update_EE_FG(I_EE_full, I_EE_Task_kk);
  Namen_Methoden{kk} = name_method;
  Namen_Methoden_Leg_Paper{kk} = name_method_leg;
  % Debug: Check the computation time in profiler
  t1 = tic();
  RP.invkin2_traj(X_t_kk(1:10,:), XD_t_kk(1:10,:), XDD_t_kk(1:10,:), t(1:10,:), qs_kk, set_kk);
  fprintf(['Tested Trajectory execution for 10 samples. Avg %1.1fms per ', ...
    'sample. Estimated %1.1fmin for full trajectory of %d samples.\n'], ...
    1e3*toc(t1)/10, toc(t1)/10*length(t)/60, length(t));
  assert(all(~isnan(X_t_kk(:))), 'Trajectory contains NaN. Syntax error');
  % Compute the whole trajectory
  comptimes_traj = NaN(11,1);
  for iit = 1:(1+10*usr_test_timing) % repeat computation for timing measurements
  t1=tic();
  [Q_t_kk, QD_t_kk, QDD_t_kk, Phi_t_kk,~,~,~,Stats_kk] = RP.invkin2_traj( ...
    X_t_kk, XD_t_kk, XDD_t_kk, t, qs_kk, set_kk);
  comptimes_traj(iit) = toc(t1);
  end
  I_err = abs(Phi_t_kk) > max(set_kk.Phit_tol,set_kk.Phir_tol) | isnan(Phi_t_kk);
  % Prüfe, ob die Trajektorie vorzeitig abbricht. Das ist kein Fehler, da
  % bei ungünstiger Parametrierung des IK-Algorithmus keine Lösung
  % gefunden werden kann. Bei guter Parametrierung ist dann eine Lösung
  % möglich.
  n_iO = n; % Anzahl der i.O. Trajektorienpunkte
  if any(I_err(:))
    I1 = find(sum(I_err,2),1);
    n_iO = I1-1; % Bis einen Punkt vor dem Fehler war es noch i.O.
    warning(['Fehler in Trajektorie zu groß. Zuerst bei Zeitschritt %d/%d ', ...
      '(t=%1.3fs). Traj.-IK nicht vollständig berechenbar'], I1, size(Q_t_kk,1), t(I1));
  end
  fprintf(['Traj.-IK Fall %d (%s) berechnet. Dauer: %1.1fs für %d/%d ', ...
    'Bahnpunkte. %1.1fms pro i.O.-Punkt\n']', kk, name_method, comptimes_traj(1), n_iO, ...
    n, 1e3*toc(t1)/n_iO);
  if usr_test_timing
    t_ps = comptimes_traj / n_iO;
    fprintf(['Timing Evaluation of trajectory inverse kinematics:', ...
          '\n\tin total: avg %1.2fs, std %1.3fs, n=%d, data [%s]s (ignore first);', ...
          '\n\tper sample: avg %1.2fms, std %1.6fms, data [%s]ms (ignore first); %d samples\n'], ...
          mean(comptimes_traj(2:end)), std(comptimes_traj(2:end)), length(comptimes_traj(2:end)), disp_array(comptimes_traj(:)', '%1.2f'), ...
          1e3*mean(t_ps(2:end)),       1e3*std(t_ps(2:end)), disp_array(1e3*t_ps(:)', '%1.2f'), n_iO);
  end
  Q_t_all(:,:,kk) = Q_t_kk;
  QD_t_all(:,:,kk) = QD_t_kk;
  QDD_t_all(:,:,kk) = QDD_t_kk;
  % Wiederherstellung der Grenzwerte, die temporär zurückgesetzt wurden
  for j = 1:RP.NLEG
    RP.Leg(j).qlim = qlim_backup(RP.I1J_LEG(j):RP.I2J_LEG(j),:);
    RP.Leg(j).qDlim = qDlim_backup(RP.I1J_LEG(j):RP.I2J_LEG(j),:);
    RP.Leg(j).qDDlim = qDDlim_backup(RP.I1J_LEG(j):RP.I2J_LEG(j),:);
  end
  % Actual platform trajectory
  [X_ist, XD_ist] = RP.fkineEE2_traj(Q_t_kk, QD_t_kk);
  % Get platform pose from integration to avoid restriction to +/- pi
  X_ist_int = repmat(X_ist(1,:), length(t), 1) + cumtrapz(t, XD_ist);
  % Normalize platform angles from direct calculation using angles from
  % integration as center. Gives exact solution without limitation to +/-pi
  X_ist(:,4:6) = normalizeAngle(X_ist(:,4:6), X_ist_int(:,4:6));
  % Check if numeric integration matches the direct calculation
%   X_ist1norm = X_ist; X_ist(:,4:6) = normalizeAngle(X_ist(:,4:6));
%   test_X = [X_ist(:,1:3) - X_ist1norm(:,1:3), angleDiff(X_ist(:,4:6),X_ist1norm(:,4:6))];
%   assert(all(abs(test_X(:))<1e-3), ['Numeric integration of platform ', ...
%     'velocity to pose does not match direct calculation']);
%   max(abs(test_X(:)))
  % Plattform-Pose aus direkter Kinematik abspeichern
  XE_all(:,:,kk) = X_ist(:,1:6); % Nehme Plattform-Pose von erster Beinkette gesehen.
  XDE_all(:,:,kk) = XD_ist(:,1:6); % Nehme Plattform-Pose von erster Beinkette gesehen.
  % Zielfunktionen für Position
  h1 = NaN(size(H1_all(:,:,kk))); h2=h1;
  for ii = 1:length(t)
    h1(ii,1) = invkin_optimcrit_limits1(Q_t_kk(ii,:)', qlim);
    h2(ii,1) = invkin_optimcrit_limits2(Q_t_kk(ii,:)', qlim);
    for kkk = 1:RP.NLEG % Kriterien für jedes Bein einzeln berechnen
      h1(ii,1+kkk) = invkin_optimcrit_limits1(Q_t_kk(ii,RP.I1J_LEG(kkk):RP.I2J_LEG(kkk))', RP.Leg(kkk).qlim);
      h2(ii,1+kkk) = invkin_optimcrit_limits2(Q_t_kk(ii,RP.I1J_LEG(kkk):RP.I2J_LEG(kkk))', RP.Leg(kkk).qlim);
    end
  end
  H1_all(:,:,kk) = h1;
  H2_all(:,:,kk) = h2;
  % Zielfunktionen für Geschwindigkeit
  h1D = NaN(size(H1D_all(:,:,kk))); h2D=h1;
  for ii = 1:length(t)
    h1D(ii,1) = invkin_optimcrit_limits1(QD_t_kk(ii,:)', cat(1, RP.Leg.qDlim));
    h2D(ii,1) = invkin_optimcrit_limits2(QD_t_kk(ii,:)', cat(1, RP.Leg.qDlim));
    for kkk = 1:RP.NLEG % Kriterien für jedes Bein einzeln berechnen
      h1D(ii,1+kkk) = invkin_optimcrit_limits1(QD_t_kk(ii,RP.I1J_LEG(kkk):RP.I2J_LEG(kkk))', RP.Leg(kkk).qDlim);
      h2D(ii,1+kkk) = invkin_optimcrit_limits2(QD_t_kk(ii,RP.I1J_LEG(kkk):RP.I2J_LEG(kkk))', RP.Leg(kkk).qDlim);
    end
  end
  H1D_all(:,:,kk) = h1D;
  H2D_all(:,:,kk) = h2D;
  % Zielfunktion für Konditionszahl
  for jj = 1:n_iO
    % Konditionszahl der IK-Jacobi-Matrix der PKM
    Phi_q = RP.constr3grad_q(Q_t_kk(jj,:)', X_t(jj,:)');
    Hcond_all(jj,1, kk) = cond(Phi_q);
%       % Konditionszahl der Jacobi jeder einzelnen Beinkette
%       for kkk = 1:RP.NLEG % Kriterien für jedes Bein einzeln berechnen
%         Hcond_all(jj,2+kkk, kk) = cond(RP.Leg(kkk).jacobig(Q_t_kk(jj,RP.I1J_LEG(kkk):RP.I2J_LEG(kkk))'));
%       end
    % PKM-Jacobi-Matrix
    [~,Phi_q_voll] = RP.constr4grad_q(Q_t_kk(jj,:)');
    [~,Phi_x_voll] = RP.constr4grad_x(X_ist(jj,1:6)');
    Jinv_voll = -Phi_q_voll\Phi_x_voll;
    Hcond_all(jj,2, kk) = cond(Jinv_voll(RP.I_qa,RP.I_EE));
  end

  % Gelenkkoordinaten normieren
  Q_t_norm_all(:,:,kk) = (Q_t_all(:,:,kk) - repmat(qlim(:,1)',n,1)) ... % untere Grenze abziehen
                          ./ repmat(qlim(:,2)'-qlim(:,1)',n,1); % geteilt durch Spannweite
end

if ~usr_load_traj
  save(fullfile(respath, [filename_pre, '_traj.mat']), 'Q_t_all', 'QD_t_all', ...
    'QDD_t_all', 'XE_all', 'XDE_all', 'Q_t_norm_all', 'H_all', 'H1_all', ...
    'H1D_all', 'H2_all', 'H2D_all', 'Hcond_all', 'Namen_Methoden', ...
    'Namen_Methoden_Leg_Paper', 'comptimes_traj');
end
% For Debugging
save(fullfile(respath, [filename_pre, '_all_data.mat']));
%% Paper-Bild: Konditionszahl-Karte über redundante Koordinate und Traj

% Repeat distribution of criteria for pi-periodicity
% Update. Do not do this. The condition number is periodic, but not the
% joint limit criterion!
phiz_range_ext = phiz_range; %[phiz_range-2*pi, phiz_range, phiz_range+2*pi];
H_all_ext = H_all; % NaN(size(H_all,1),3*size(H_all,2),size(H_all,3));
% for k = 1:size(H_all,3)
%   H_all_ext(:,:,k) = repmat(H_all(:,:,k),1,3);
% end
fighdl = change_current_figure(2400);clf;hold on;
set(fighdl, 'Name', 'PaperFigure_PerfDistr', 'NumberTitle', 'off');
% T_tref = (1:size(X_tref,1))/size(X_tref,1);
[X_ext,Y_ext] = meshgrid(s_ref,180/pi*phiz_range_ext);
Z_ext = zeros(size(X_ext));
% Create color code from Jacobian condition number
CC_ext = H_all_ext(:,:,4)';
CC_ext = H_all_ext(:,:,4)' + H_all_ext(:,:,2)';

% high condition numers all get the same dark (or magenta) color to be able
% to better distinguish the good values.
colorlimit = 1e4;
% Information for the text in Sec. 6.2
fprintf('Condition numbers above %1.1e are excluded from the color code (magenta)\n', colorlimit);
CC_ext(CC_ext>colorlimit) = colorlimit;%1e6; % limit to control colors.
maxcond_alltraj = max(max(Hcond_all(:,2,:)));
mincond_alltraj = min(min(Hcond_all(:,2,:)));
% Saturate all values above 400 to have more colors in the range of low
% condition numbers
condsat_limit = 400;
% Information for the text in Sec. 6.2
fprintf('Condition numbers above %1.1e are saturated in the color code\n', condsat_limit);
I_exc = CC_ext > condsat_limit;
CC_ext(I_exc) = condsat_limit+10*log10(CC_ext(I_exc)/condsat_limit); % last term gives 0...30 for condsat_limit=1e3
% Remove values with invalid joint configurations
I_limviol = isinf(H_all_ext(:,:,2)');
% CC(I_limviol) = NaN;
% Create color plot
surf(X_ext,Y_ext,Z_ext,CC_ext, 'EdgeColor', 'none');
xlabel('Normalized trajectory progress $s$ (per point of support)', 'interpreter', 'latex');
ylabel('Redundant coordinate $\varphi_z$ in deg', 'interpreter', 'latex');
colorbar
view([0,90])
% Set Colormap: low condition numbers white, high/singularity dark red.
colors_map = flipud(hot(1024));
numcolors_sat = 20;  % add magenta for worst values. This determines the size of the magenta color range
colors_map = [colors_map; repmat([255 0 255]/255, numcolors_sat,1)];
cm = colormap(colors_map);
% Debug: Other colormaps
% colormap(jet(1024));
% colormap(summer(1024));
% colormap(winter(1024));
% colormap(flipud(parula(1024)))
% Use Linewidth 1 and MarkerSize 6
plot_lw = 1.5;
plot_ms = 6; % default
format = {[0 255 255]/255,  '', '-', 12, plot_lw, plot_ms; ... % cyan
          [0 80 155 ]/255, '^', '-', 16, plot_lw, plot_ms; ... %imesblau
          [231 123 41 ]/255, 'o', '-', 5, plot_lw, plot_ms; ... %imesorange
          'r', 'x', '--', 7, plot_lw, plot_ms; ...
          'k', 's',  '--', 9, plot_lw, plot_ms; ...
          'g', 'v', '-', 34, plot_lw, plot_ms};

linhdl = NaN(length(Namen_Methoden), 1);

for kk = 1:length(Namen_Methoden)
  X_ist = XE_all(:,1:6,kk);
  
  s_kk = NaN(size(X_ist,1),1);
  % Get progress related to key points of the trajectory
  for i = 1:length(IL)
    I1 = IL(i);
    if i < length(IL), I2 = IL(i+1);
    else,              I2 = size(X_ist,1);
    end
    p_all = ( X_ist(I1:I2,1:5)-repmat(XL(i,1:5),I2-I1+1,1) ) ./ ...
             repmat(XL(i+1,1:5)-XL(i,1:5),I2-I1+1,1);
    p_all(isinf(p_all))=NaN;
    p = mean(p_all,2, 'omitnan');
    s_kk(I1:I2) = (i-1)+p;
  end
  % insert trajectory into plot
  change_current_figure(2400);
  linhdl(kk) = plot(s_kk, 180/pi*X_ist(:,6), 'g-'); % style will be overwritten
end
xlim([0, ceil(max(s_ref))]);
ylim([-120, 160]);
linhdl_leg = line_format_publication(linhdl, format);
cb = colorbar();
cbyh = ylabel(cb,'Performance criterion (condition number, joint limits)', ...
  'Rotation',90, 'interpreter', 'latex');
% Namen_Methoden_Leg_Paper = {'ref.','$\varphi_{z,0}{=}{-}30^\circ$', ...
%   '$\varphi_{z,0}{=}45^\circ$','$\varphi_{z,0}{=}90^\circ$','$\varphi_{z,0}{=}45^\circ$'};
Namen_Methoden_Leg_Paper = {'ref.','set 1', 'set 2','set 3','set 4'};
lh = legend(linhdl_leg, Namen_Methoden_Leg_Paper, 'interpreter', 'latex');
% text with subfig number (a)
thdl = text(0,0,'(a)');
[x_off, x_slope] = get_relative_position_in_axes(gca, 'x');
[y_off, y_slope] = get_relative_position_in_axes(gca, 'y');
set(thdl, 'Position', [x_off+x_slope*(-1.2), y_off+y_slope*(-1.2), 0]);
set(thdl, 'FontWeight', 'bold');
% adjust axes labels text
ylh = get(gca, 'ylabel');
set(ylh, 'Position', [x_off+x_slope*(-1.15), y_off+y_slope*0, 0]);
% color bar text may look different in pdf than in figure (if not using
% latex interpreter above). Move label to the left to avoid cutting it off.
set(cbyh, 'HorizontalAlignment', 'left');
set(cbyh, 'Position', [2.3, 0, 0]);
figure_format_publication(fighdl);
set(gca, 'Box', 'off');
set(fighdl, 'windowstyle', 'normal');
set_size_plot_subplot(fighdl, ...
  11,7,gca,...
  0.09,0.18,0.1,0.12,... %l r u d
  0,0) % x y
drawnow();
% Update Legend size after changing the figure size. Otherwise the legend
% size sometimes becomes unexpectedly large
set(lh, 'orientation', 'horizontal', 'position', [0.1,0.92,0.8,0.05]);
drawnow();
if usr_save_figures
t1 = tic();
% For this to work, the Java heap memory has to be high enough. For
% high-resolution image 2GB not enough, 4GB worked.
% https://de.mathworks.com/matlabcentral/answers/255083-matlab-and-robotic-system-toolbox-java-lang-outofmemoryerror-gc-overhead-limit-exceeded#answer_318119
exportgraphics(fighdl, fullfile(paperfig_path, ['nullspace_traj','.pdf']),'ContentType','vector') % ,'Resolution','100'
fprintf('Exported performance map as vector graphics. Duration: %1.1fs\n', toc(t1));
% exportgraphics(fighdl, fullfile(paperfig_path, ['nullspace_traj','.png']),'Resolution','600')  
% export_fig(fullfile(paperfig_path, ['nullspace_traj_exportfig','.pdf']), fighdl)
end

%% Paper-Bild: Verlauf der Konditionszahl für verschiedene Trajektorien
fighdl = change_current_figure(2500);clf;
Hcond_PKMJac_all = reshape(Hcond_all(:,2,:), n, length(Namen_Methoden));
hdl=plot(s_kk, Hcond_PKMJac_all);
leghdl=line_format_publication(hdl, format);
grid on;
xlabel('Norm. trajectory progress $s$', 'interpreter', 'latex');
% ylabel('Jacobian condition number $\mathrm{cond}(\boldmath{J})$', 'interpreter', 'latex');
set(gca, 'xtick', 0:7);
set(gca, 'ytick', [50, 75, 100, 150, 200:100:700]);
set(gca, 'YScale', 'log')
xlim([0, max(s_kk)+1e-3])
ylim([50, 600])
% text with subfig number (b)
thdl = text(0,0,'(b)');
[x_off, x_slope] = get_relative_position_in_axes(gca, 'x');
[y_off, y_slope] = get_relative_position_in_axes(gca, 'y');
set(thdl, 'Position', [x_off+x_slope*(-1.33), y_off+y_slope*(-1.04), 0]);
set(thdl, 'FontWeight', 'bold');
figure_format_publication(fighdl);
set_size_plot_subplot(fighdl, ...
  4.6,7,gca,...
  0.14,0.02, 0.05,0.12,... %l r u d
  0,0) % x y
drawnow();
if usr_save_figures
%   saveas(2500, fullfile(respath, 'condition_number_traj.fig'));
%   saveas(2500, fullfile(respath, 'condition_number_traj.png'));
  exportgraphics(fighdl, fullfile(paperfig_path, ['nullspace_traj_condition','.pdf']),'ContentType','vector')
%   exportgraphics(gcf, fullfile(paperfig_path, ['nullspace_traj_condition','.png']),'Resolution','600')
end
% legend(leghdl, Namen_Methoden);

%% Bild: Schnitte verschiedener Bahnkoordinaten s durch die Karte
if debug_plot
  change_current_figure(2600);clf; hold all;
  s_select = 5.2:0.05:5.5;
  for i = 1:length(s_select)
    [~,II_i] = min(abs(s_ref-s_select(i)));
  %   s_ref(II_i)
    plot(180/pi*phiz_range, H_all(II_i,:,4));
  end
  xlim([-120, 160]);
  set(gca, 'YScale', 'log')
  ylim([50, 1000]);
  grid on;
  xlabel('platform rotation (redundant coordinate)');
  ylabel('condition number');
  title('Horizontal slices of the performance map');
  legtext = cell(1,length(s_select));
  for i = 1:length(s_select)
    legtext{i} = sprintf('$s{=}%1.2f$', s_select(i));
  end
  legend(legtext, 'interpreter', 'latex');
  
  change_current_figure(2601);clf; hold all;
  phiz_select = sort([60,30,0,(-15:-5:-30)]*pi/180);
  for i = 1:length(phiz_select)
    [~,II_i] = min(abs(phiz_range-phiz_select(i)));
  %   s_ref(II_i)
    plot(s_ref, H_all(:,II_i,4));
  end
  xlim([0, 7]);
  set(gca, 'YScale', 'log')
  xlabel('normalized trajectory progress');
  ylim([50, 1000]);
  title('Vertical slices of the performance map');
  grid on;
  ylabel('condition number');
  legtext = cell(1,length(phiz_select));
  for i = 1:length(phiz_select)
    legtext{i} = sprintf('$\\varphi_z{=}%1.0f^\\circ$', 180/pi*phiz_select(i));
  end
  legend(legtext, 'interpreter', 'latex');
end

%% Konsistenz von Position, Geschwindigkeit und Beschleunigung testen
for kk = 1:length(Namen_Methoden)
  Q = Q_t_all(:,:,kk);
  QD = QD_t_all(:,:,kk);
  QDD = QDD_t_all(:,:,kk);
  % Integriere die Beschleunigung und Geschwindigkeit
  QD_int = cumtrapz(t, QDD) + repmat(QD(1,:),n,1);
  % Integriere die Geschwindigkeit zur Position
  Q_int = cumtrapz(t, QD) + repmat(Q(1,:),n,1);
  % Differenziere die Geschwindigkeit zur Beschleunigung
  QDD_diff = [diff(QD)./repmat(diff(t),1,size(QDD,2));zeros(1,size(QDD,2))];
  % Rechnerischer Vergleich der Verläufe
  corrQD = diag(corr(QD_int, QD));
  corrQ = diag(corr(Q_int, Q));
  corrQDD = diag(corr(QDD_diff, QDD));
  conserr = false;
  if any(corrQDD < 0.95) || any(corrQD < 0.95) || any(corrQ < 0.98)
    conserr = true;
  end
  if conserr
    % Vergleiche die Verläufe graphisch: Position
    change_current_figure(100*robnr+40+kk);clf;
    set(100*robnr+40+kk, 'Name', sprintf('Rob%d_M%d_Kons_q_qD', robnr, kk), 'NumberTitle', 'off');
    sgtitle(sprintf('Konsistenz q-int(qD) M%d (%s)', kk, Namen_Methoden{kk}), 'interpreter', 'none');
    ii = 0;
    for i = 1:RP.NLEG
      for j = 1:RP.Leg(i).NJ
        ii = ii + 1;
        subplot(RP.NLEG,RP.Leg(1).NJ,sprc2no(RP.NLEG, RP.Leg(1).NJ, i, j));
        hold on; grid on;
        hdl1=stairs(t, Q(:,ii));
        hdl2=stairs(t, Q_int(:,ii), '--');
        ylim([-0.1, +0.1]+minmax2([Q(:,ii);Q_int(:,ii)]')); % Gelenkgrenzen nicht für Plot-Grenzen berücksichtigen
        plot(t([1,end]), repmat(RP.Leg(i).qlim(j,:),2,1), 'r-');
        if j == 1, ylabel(sprintf('BK %d', i)); end
        if i == 1, title(sprintf('q %d', j)); end
        if ii == RP.NJ, legend([hdl1;hdl2],{'direkt', 'integral'}); end
      end
    end
    linkxaxes
    if usr_save_figures
      saveas(100*robnr+40+kk, fullfile(respath,sprintf('Rob%d_M%d_Konsistenz_q', robnr, kk)));
    end
    % Vergleiche die Verläufe graphisch: Geschwindigkeit
    change_current_figure(100*robnr+50+kk);clf;
    set(100*robnr+50+kk, 'Name', sprintf('Rob%d_M%d_Kons_qD_qDD', robnr, kk), 'NumberTitle', 'off')
    sgtitle(sprintf('Konsistenz qD-int(qDD) M%d (%s)', kk, Namen_Methoden{kk}), 'interpreter', 'none');
    ii = 0;
    for i = 1:RP.NLEG
      for j = 1:RP.Leg(i).NJ
        ii = ii + 1;
        subplot(RP.NLEG,RP.Leg(1).NJ,sprc2no(RP.NLEG, RP.Leg(1).NJ, i, j));
        hold on; grid on;
        hdl1=stairs(t, QD(:,ii));
        hdl2=stairs(t, QD_int(:,ii), '--');
        plot(t([1,end]), repmat(RP.Leg(i).qDlim(j,:),2,1), 'r-');
        ylim([-0.1, +0.1]+minmax2([QD(:,ii);QD_int(:,ii)]')); % Gelenkgrenzen nicht für Plot-Grenzen berücksichtigen
        if j == 1, ylabel(sprintf('BK %d', i)); end
        if i == 1, title(sprintf('qD %d', j)); end
        if ii == RP.NJ, legend([hdl1;hdl2],{'direkt', 'integral'}); end
      end
    end
    linkxaxes
    if usr_save_figures
      saveas(100*robnr+50+kk, fullfile(respath,sprintf('Rob%d_M%d_Konsistenz_qD', robnr, kk)));
    end
    % Vergleiche die Verläufe graphisch: Beschleunigung
    change_current_figure(100*robnr+60+kk);clf;
    set(100*robnr+60+kk, 'Name', sprintf('Rob%d_M%d_Kons_qDD', robnr, kk), 'NumberTitle', 'off')
    sgtitle(sprintf('Konsistenz diff(qD)-qDD M%d (%s)', kk, Namen_Methoden{kk}), 'interpreter', 'none');
    ii = 0;
    for i = 1:RP.NLEG
      for j = 1:RP.Leg(i).NJ
        ii = ii + 1;
        subplot(RP.NLEG,RP.Leg(1).NJ,sprc2no(RP.NLEG, RP.Leg(1).NJ, i, j));
        hold on; grid on;
        hdl1=stairs(t, QDD(:,ii));
        hdl2=stairs(t, QDD_diff(:,ii), '--');
        plot(t([1,end]), repmat(RP.Leg(i).qDDlim(j,:),2,1), 'r-');
        if j == 1, ylabel(sprintf('BK %d', i)); end
        if i == 1, title(sprintf('qDD %d', j)); end
        if ii == RP.NJ, legend([hdl1;hdl2],{'direkt', 'differenz'}); end
      end
    end
    linkxaxes
    if usr_save_figures
      saveas(100*robnr+50+kk, fullfile(respath,sprintf('Rob%d_M%d_Konsistenz_qDD', robnr, kk)));
    end
  end
end
%% Trajektorie: Bild für einzelne Beine
if debug_plot
  fprintf('Zeichne Bilder für Zusammenfassung von einzelnen Beinketten (Trajektorien-IK)\n');
  for kk = 1:length(Namen_Methoden)
    % Definitionen für die Ergebnisse dieser Methode laden
    Q_t = Q_t_all(:,:,kk);
    QD_t = QD_t_all(:,:,kk);
    QDD_t = QDD_t_all(:,:,kk);
    Q_t_norm = Q_t_norm_all(:,:,kk);
    Name = Namen_Methoden{kk};
    H1_t = H1_all(:,:,kk);
    H2_t = H2_all(:,:,kk);
    % Zeichnen
    change_current_figure(100*robnr+10+kk);clf;
    set(100*robnr+10+kk, 'Name', sprintf('Rob%d_M%d_Beine', robnr, kk), 'NumberTitle', 'off')
    sgtitle(sprintf('Trajektorie Q (Beine); %s', Name));
    axhdl = NaN(5,RP.NLEG);
    for i = 1:RP.NLEG
      % Gelenkkoordinaten
      axhdl(1,i)=subplot(5,RP.NLEG,sprc2no(5,RP.NLEG,1,i));
      plot(t, Q_t(:,RP.I1J_LEG(i):RP.I2J_LEG(i)));
      ylabel(sprintf('q_%d', i)); grid on;
      if i == RP.NLEG % letzte Spalte der Subplots
        l = {};
        for j = 1:RP.Leg(1).NJ
          l = [l, {sprintf('q_%d', j)}];
        end
        legend(l);
      end
      title(sprintf('Beinkette %d', i));
      % Gelenkgeschwindigkeiten
      axhdl(2,i)=subplot(5,RP.NLEG,sprc2no(5,RP.NLEG,2,i));
      plot(t, QD_t(:,RP.I1J_LEG(i):RP.I2J_LEG(i)));
      ylabel(sprintf('qD_%d', i)); grid on;
      % Gelenkbeschleunigungen
      axhdl(3,i)=subplot(5,RP.NLEG,sprc2no(5,RP.NLEG,3,i));
      plot(t, QDD_t(:,RP.I1J_LEG(i):RP.I2J_LEG(i)));
      ylabel(sprintf('qDD_%d', i)); grid on;
      % Debug: Normierte Gelenk-Koordinaten.
      % plot(t, Q_t_norm(:,RP.I1J_LEG(i):RP.I2J_LEG(i)));
      % ylabel(sprintf('q_%d (norm)', i)); grid on;
      % Zielfunktion 1
      axhdl(4,i)=subplot(5,RP.NLEG,sprc2no(5,RP.NLEG,4,i));
      plot(t, H1_t(:,1+i)); grid on;
      ylabel(sprintf('ZF 1 BK %d', i));
      % Zielfunktion 2
      axhdl(5,i)=subplot(5,RP.NLEG,sprc2no(5,RP.NLEG,5,i));
      plot(t, H2_t(:,1+i)); grid on;
      ylabel(sprintf('ZF 2 BK %d', i));
    end
    linkxaxes
    remove_inner_labels(axhdl,'x')
    if usr_save_figures
      saveas(100*robnr+10+kk, fullfile(respath,sprintf( ...
      'Rob%d_M%d_%s_Trajektorie_Beine_%s.fig',robnr, kk, RP.mdlname, Name)));
    end
  end
end
%% Trajektorie: Bild für Gesamt-Zielfunktionen
if debug_plot
  fprintf('Zeichne Bilder für Zielfunktionen (Trajektorien-IK)\n');
  % Zielfunktionen bezogen auf Positions-Grenzen
  change_current_figure(100*robnr+20); clf;
  set(100*robnr+20, 'Name', sprintf('Rob%d_Zielf_q', robnr), 'NumberTitle', 'off');
  sgtitle('Zielfunktionen (Position)');
  H1_PKM_all = reshape(H1_all(:,1,:), n, length(Namen_Methoden));
  H2_PKM_all = reshape(H2_all(:,1,:), n, length(Namen_Methoden));
  subplot(2,2,1); hold on;
  hdl={};
  hdl{1}=plot(t, H1_PKM_all);
  title('Zielfunktion 1 (nicht opt.)');
  ylabel('Zielfkt 1'); grid on;
  subplot(2,2,3); hold on;
  hdl{2}=plot(t, log10(H1_PKM_all));
  ylabel('Log.-Zielfkt 1 (n.o.)'); grid on;
  subplot(2,2,2); hold on;
  hdl{3}=plot(t, H2_PKM_all);
  title('Optimierungs-Zielfunktion 2');
  ylabel('Zielfkt 2'); grid on;
  subplot(2,2,4); hold on;
  hdl{4}=plot(t, log10(H2_PKM_all));
  ylabel('Log.-Zielfkt 2'); grid on;
  if usr_save_figures
    saveas(100*robnr+20, fullfile(respath,sprintf('Rob%d_Zielf_Pos', robnr)));
  end
  % Linien nachträglich neu formatieren (bessere Lesbarkeit)
  for k = 1:4, leghdl=line_format_publication(hdl{k}, format_mlines); end
  legend(leghdl, Namen_Methoden);
  linkxaxes

  % Zielfunktionen bezogen auf Geschwindigkeits-Grenzen
  change_current_figure(100*robnr+21); clf;
  set(100*robnr+21, 'Name', sprintf('Rob%d_Zielf_qD', robnr), 'NumberTitle', 'off');
  sgtitle('Zielfunktionen (Geschwindigkeit)');
  H1D_PKM_all = reshape(H1D_all(:,1,:), n, length(Namen_Methoden));
  H2D_PKM_all = reshape(H2D_all(:,1,:), n, length(Namen_Methoden));
  subplot(2,2,1); hold on;
  hdl={};
  hdl{1}=plot(t, H1D_PKM_all);
  title('Zielfunktion 1 (nicht opt.)');
  ylabel('Zielfkt 1'); grid on;
  subplot(2,2,3); hold on;
  hdl{2}=plot(t, log10(H1D_PKM_all));
  ylabel('Log.-Zielfkt 1 (n.o.)'); grid on;
  subplot(2,2,2); hold on;
  hdl{3}=plot(t, H2D_PKM_all);
  title('Optimierungs-Zielfunktion 2');
  ylabel('Zielfkt 2'); grid on;
  subplot(2,2,4); hold on;
  hdl{4}=plot(t, log10(H2D_PKM_all));
  ylabel('Log.-Zielfkt 2'); grid on;
  if usr_save_figures
    saveas(100*robnr+21, fullfile(respath,sprintf('Rob%d_Zielf_Geschw', robnr)));
  end
  % Linien nachträglich neu formatieren (bessere Lesbarkeit)
  for k = 1:4, leghdl=line_format_publication(hdl{k}, format_mlines); end
  legend(leghdl, Namen_Methoden);
  linkxaxes

  % Zielfunktionen bezogen auf Konditionszahl
  change_current_figure(100*robnr+22); clf;
  set(100*robnr+22, 'Name', sprintf('Rob%d_Zielf_cond', robnr), 'NumberTitle', 'off');
  subplot(2,1,1); hold on;
  Hcond_IKJac_all = reshape(Hcond_all(:,1,:), n, length(Namen_Methoden));
  hdl=plot(t, Hcond_IKJac_all);
  leghdl=line_format_publication(hdl, format_mlines);
  legend(leghdl, Namen_Methoden);
  grid on; ylabel('cond(Phi_q)', 'interpreter', 'none');
  subplot(2,1,2); hold on;
  Hcond_PKMJac_all = reshape(Hcond_all(:,2,:), n, length(Namen_Methoden));
  hdl=plot(t, Hcond_PKMJac_all);
  leghdl=line_format_publication(hdl, format_mlines);
  legend(leghdl, Namen_Methoden);
  grid on; ylabel('cond(J)', 'interpreter', 'none');
  sgtitle('Konditionszahl (mögliche Opt. Zielf.)');
  linkxaxes
  if usr_save_figures
    saveas(100*robnr+22, fullfile(respath,sprintf('Rob%d_Zielf_cond', robnr)));
  end
end
%% Trajektorie: Bild für Plattformbewegung
if debug_plot
  change_current_figure(100*robnr+24);clf;
  set(100*robnr+24, 'Name', sprintf('Rob%d_TrajX', robnr), 'NumberTitle', 'off')
  sgtitle('Trajektorie X (Details)');
  for i = 1:6
    subplot(3,2,i); hold on
    linhdl1=plot(t, reshape(XE_all(:,i,:),n,length(Namen_Methoden)));
    linhdl3=plot(t, X_t(:,i), ':');
    % Eckpunkte einzeichnen
    plot(t(IL), X_t(IL,i), 'ko');
    if i < 4, unit = 'm';
    else, unit = 'rad';
    end
    ylabel(sprintf('x %d in %s', i, unit));
    grid on
    leghdl1=line_format_publication(linhdl1, format_mlines);
  %     if i == 5, legend([leghdl1;linhdl3], [Namen_Methoden,{'Soll 3T3R'}]); end
  end
  linkxaxes
  if usr_save_figures
    saveas(100*robnr+24, fullfile(respath,sprintf( ...
      'Rob%d_%s_TrajX.fig',robnr, RP.mdlname)));
  end
end
%% Trajektorie: Vergleich Gelenkverlauf nach verschiedenen Methoden
if debug_plot
  for Dtype = 1:3
    change_current_figure(100*robnr+70+Dtype); clf;
    if     Dtype == 1, Dstring = ''; %#ok<ALIGN>
    elseif Dtype == 2, Dstring = 'D';
    elseif Dtype == 3, Dstring = 'DD'; end
    set(100*robnr+70+Dtype, 'Name', sprintf('Rob%d_Q%s', robnr, Dstring), 'NumberTitle', 'off');
    sgtitle(sprintf('Q %s', Dstring));
    hdl = NaN(RP.NLEG,RP.Leg(1).NJ,length(Namen_Methoden));
    axhdl = NaN(RP.NLEG,RP.Leg(1).NJ);
    for kk = 1:length(Namen_Methoden)
      Q = Q_t_all(:,:,kk);
      QD = QD_t_all(:,:,kk);
      QDD = QDD_t_all(:,:,kk);
      ii = 0;
      for i = 1:RP.NLEG
        for j = 1:RP.Leg(i).NJ
          ii = ii + 1;
          axhdl(i,j)=subplot(RP.NLEG,RP.Leg(1).NJ,sprc2no(RP.NLEG, RP.Leg(1).NJ, i, j));
          hold on; grid on;
          if Dtype == 1
            hdl(i,j,kk)=stairs(t, Q(:,ii));
            plot(t([1,end]), repmat(RP.Leg(i).qlim(j,:),2,1), 'r-');
            ylim([-0.1, +0.1]+qlim(ii,:));
          elseif Dtype == 2
            hdl(i,j,kk)=stairs(t, QD(:,ii));
            plot(t([1,end]), repmat(RP.Leg(i).qDlim(j,:),2,1), 'r-');
            ylim([-0.1, +0.1]+qDlim(ii,:));
          elseif Dtype == 3
            hdl(i,j,kk)=stairs(t, QDD(:,ii));
            plot(t([1,end]), repmat(RP.Leg(i).qDDlim(j,:),2,1), 'r-');
            ylim([-0.1, +0.1]+qDDlim(ii,:));
          end
          if j == 1, ylabel(sprintf('BK %d', i)); end
          if i == 1, title(sprintf('q%s %d', Dstring, j)); end
        end
      end
    end
    linkxaxes
    % Linien nachträglich neu formatieren (bessere Lesbarkeit)
    for i = 1:RP.NLEG
      for j = 1:RP.Leg(i).NJ
        leghdl=line_format_publication(hdl(i,j,:), format_mlines);
      end
    end
    remove_inner_labels(axhdl,1);
    legend(leghdl, Namen_Methoden);
    if usr_save_figures
      saveas(100*robnr+70+Dtype, fullfile(respath,sprintf('Rob%d_Q%s', robnr, Dstring)));
    end
  end
end

%% Animation des bewegten Roboters
if usr_create_anim
RP.I_EE_Task = I_EE_red; % zum Zeichnen, damit x(6) ignoriert wird
for kk = 1:length(Namen_Methoden)
for i_dur = 1:(1+usr_anim_realtime) % create two videos (two durations)
  if i_dur == 1
    maxduration_animation = 10; % Dauer der Animation als mp4 (in s)
  else
    maxduration_animation = t(end); % Video in Echtzeit
  end
  Q_t = Q_t_all(:,:,kk);
  Name = Namen_Methoden{kk};
  anim_filename = fullfile(respath, sprintf('nullspace_traj_M%d_%s_duration%1.0fs', ...
    kk, Namen_Methoden{kk}, maxduration_animation));
  s_anim = struct( 'mp4_name', [anim_filename,'.mp4'] );
  s_plot = struct( 'ks_legs', [], 'straight', 0, 'mode', 4);
  fhld_kk = change_current_figure(100*robnr+30);clf;hold all;
  if strcmp(get(fhld_kk, 'windowstyle'), 'docked')
    set(fhld_kk, 'windowstyle', 'normal');
  end
  set(fhld_kk, 'name', 'Anim', ...
    'color','w', 'NumberTitle', 'off', 'units','normalized',...
    'outerposition',[0 0 1 1]); % Vollbild, damit GIF größer wird
  view(3); axis auto; hold on; grid on;
  xlabel('x in m'); ylabel('y in m'); zlabel('z in m');
  zlim([min(Traj_W.X(:,3))-0.2, RP.T_W_0(3,4)+0.2])
  Traj_0 = struct('X', X_t);
  Traj_W = RP.transform_traj(Traj_0);
  plot3(Traj_W.X(:,1), Traj_W.X(:,2), Traj_W.X(:,3));
  t_Vid = (0:1/30*(t(end)/maxduration_animation):t(end))';
  I_anim = knnsearch( t , t_Vid ); % Berechne Indizes in Traj.-Zeitstempeln
  RP.anim( Q_t(I_anim,:), X_t(I_anim,:), s_anim, s_plot);
  fprintf('Animation der Bewegung gespeichert: %s\n', s_anim.mp4_name);
end
end
end

fprintf('Test für %s beendet\n',RP.mdlname);
