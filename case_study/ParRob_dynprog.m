% Perform nullspace motion of a hexapod robot with dynamic programming. 
% This creates the results of the LNEE paper and Fig. 1, 11, 12, 13.
% 
% See also: ParRob_nullspace_trajectory.m from this directory

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-02
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clear
clc

if isempty(which('serroblib_path_init.m'))
  error('The serial robot database is not initialized in Matlab.');
end
if isempty(which('parroblib_path_init.m'))
  error('The parallel robot database is not initialized in Matlab.');
end
rng(0); % Für Reproduzierbarkeit
%% User Settings and Initialization
%#ok<*UNRCH>
use_mex_functions = true; % use mex functions (much faster)
usr_recreate_mex = false; % recreate and recompile mex functions from templates
usr_short_traj = false; % Trajektorie stark abkürzen, um prinzipielle Funktionalität zu zeigen
usr_create_anim = false; % create an animation video of the robot motion
usr_anim_realtime = false; % save real-time animation (video as long as trajectory in seconds)
usr_highres_distrfig = true; % high resolution of the paper figure for performance criterion map
debug_plot = false;% Create debug plots
usr_plot_robot = false; % Robot image for paper.
usr_save_figures = false; % save figures to disk
usr_load_discretization = true; % load a previously computed performance map (if trajectory stays the same)
usr_load_dynprog = true;
usr_load_traj = true; % load a previously computed joint space trajectory
respath = fileparts(which('ParRob_dynprog.m'));
assert(~isempty(respath), 'Aktuelles Matlab-Skript muss im Pfad sein');
paperfig_path = fullfile(respath, '..', 'paper_LNEE', 'figures');
data_path = fullfile(respath, 'data_LNEE');
mkdirs(data_path);
usr_only_test_keypoints = false; % if true: Abort script after key point computation
usr_only_global_discretization = false; % if true: Abort script after performance map computation
usr_test_timing = false; % multiple calls of ik functions to determine runtime for evaluation in Sec. 6.1.3.
usr_trajnum = 1; % Switch between two examples in the paper (Fig. 11 vs Fig 13).

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
  mexerror = false;
  mexerror=mexerror||matlabfcn2mex({[RP.Leg(1).mdlname,'_invkin_eulangresidual']});
  mexerror=mexerror||matlabfcn2mex({[RP.Leg(1).mdlname,'_invkin_traj']});
  mexerror=mexerror||matlabfcn2mex({[RP.mdlname(1:end-6), '_invkin']});
  mexerror=mexerror||matlabfcn2mex({[RP.mdlname(1:end-6), '_invkin3']});
  mexerror=mexerror||matlabfcn2mex({[RP.mdlname(1:end-6), '_invkin_traj']});
  assert(~mexerror, 'Fehler beim Kompilieren')
else
  serroblib_update_template_functions({RP.Leg(1).mdlname});
  parroblib_update_template_functions({RP.mdlname(1:end-2)});
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
  % Setze Geschwindigkeit und Beschleunigung auf mittel-hohe Werte.
  % Sonst wird die Nullraumbewegung hiervon gebremst. Sollte nicht
  % langsamer als die Plattform-Drehung xDlim sein.
  RP.Leg(i).qDlim = repmat([-180, 180]*pi/180, RP.Leg(i).NQJ, 1); % 45deg/s is moderately high
  RP.Leg(i).qDlim(RP.Leg(i).MDH.sigma==1,:) = ...
    repmat([-2, 2], sum(RP.Leg(i).MDH.sigma==1), 1); % 2m/s
  % Set moderately high acceleration limits which could be feasible for a
  % high-performance parallel robot.
  RP.Leg(i).qDDlim = repmat([-20, 20], RP.Leg(i).NQJ, 1);
  RP.Leg(i).qDDlim(RP.Leg(i).MDH.sigma==1,:) = ...
    repmat([-20, 20], sum(RP.Leg(i).MDH.sigma==1), 1);
  % Setze Geschwindigkeit und Beschleunigung der Plattform ähnlich
  RP.xDlim = [NaN(5,2); [-1, +1]*pi]; % 180°/s
  RP.xDDlim = [NaN(5,2); [-1, +1]*pi/0.1]; % 180°/s in 0.1s
end

%% Definition der Kollisionskörper
collbodies_empty = struct( ...
        'link', zeros(0,2), ... % nx1 uint8, Nummer des zugehörigen Segments (0=Basis)
        'type', zeros(0,1), ... % nx1 uint8, Art des Ersatzkörpers
        'params', zeros(0,10)); % Parameter des jeweiligen Ersatzkörpers
% Kollisionskörper der Beinketten eintragen
for j = 1:RP.NLEG
  collbodies_j = collbodies_empty;
  for i = 1:RP.Leg(j).NJ
    % Prüfe, ob es überhaupt eine Verbindung gibt
    if RP.Leg(j).MDH.d(i)==0 && RP.Leg(j).MDH.sigma(i)==0 && RP.Leg(j).MDH.a(i)==0
      continue
    end
    % Schräge Verbindung mit Kapseln in Matlab-Klasse berechnen
    collbodies_j.link = [collbodies_j.link; uint8([i,i-1])];
    collbodies_j.type = [collbodies_j.type; uint8(6)];
    % Kapsel mit Radius 20mm
    collbodies_j.params = [collbodies_j.params; [20e-3,NaN(1,9)]];
  end
  RP.Leg(j).collbodies = collbodies_j;
end
% Kollisionskörper der Plattform eintragen
% Indizes der jeweiligen vorherigen benachbarten Beinkette
I1 = (1:RP.NLEG)'; I2 = [RP.NLEG, 1:RP.NLEG-1]';
% Variable für Kollisionsobjekte vorbereiten:
collbodies = collbodies_empty;
I_cb_base = 0;
% Kollisionsobjekte für die Plattform. Kapseln für jeden virtuellen
% Körper der Plattform-Koppelgelenke (auf Plattform-Seite). Kapsel als
% Verbindung zum jeweils vorherigen Koppelgelenk. Erzeugt Ring an der
% Plattform
collbodies.link = [collbodies.link; ...
  uint8([RP.I2L_LEG(I1)-(I1-1)-1, RP.I2L_LEG(I2)-(I2-1)-1])];
collbodies.type = [collbodies.type; repmat(uint8(6),RP.NLEG,1)];
collbodies.params = [collbodies.params; ... % Kapsel mit Radius 10mm
  [repmat(10e-3, RP.NLEG, 1), NaN(RP.NLEG, 9)]];
I_cb_platform = I_cb_base(end)+1:size(collbodies.type,1); % Ind. der Platf.-Koll.-körper
% Eintragen in Roboter-Klasse
RP.collbodies_nonleg = collbodies;
% Aktualisiere die Gesamt-Variable
RP.update_collbodies();
fprintf(['Collision bodies of the robot: %d in total, including %d of base ', ...
  'and platform\n'], size(RP.collbodies.link,1), size(RP.collbodies_nonleg.link,1));
RP.update_collchecks();
fprintf('Collision checks: %d\n', size(RP.collchecks,1));

% Startpose bestimmen
% Mittelstellung im Arbeitsraum
X0 = [ [0.00;0.00;0.5]; [0;0;0]*pi/180 ];
q0 = 0.5+rand(RP.NJ,1); % Startwerte für numerische IK (zwischen 0.5 und 1.5 rad)
q0(RP.MDH.sigma==1) = 0.5; % mit Schubaktor größer Null anfangen (damit Konfiguration nicht umklappt)
[qs, Phis, Tc_s] = RP.invkin_ser(X0, q0, struct('retry_on_limitviol',true));
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
  RP.Leg(i).qlim(1:2,:) = q0_i(1:2) + 100*repmat([-60, 60]*pi/180, 2, 1);
  RP.Leg(i).qlim(4:6,:) = q0_i(4:6) + 100*repmat([-60, 60]*pi/180, 3, 1);
end
% qlim für gesamte PKM festlegen
qlim   = cat(1, RP.Leg.qlim);
qDlim  = cat(1, RP.Leg.qDlim);
qDDlim = cat(1, RP.Leg.qDDlim);

% Roboter mit Objekten zeichnen
if false
  s_plot = struct( 'ks_legs', [], 'straight', 0, 'mode', 5, 'only_bodies', true);
  fhdl=change_current_figure(1);clf;set(fhdl,'Name','Collbodies','NumberTitle','off');
  hold on; grid on;
  xlabel('x in m'); ylabel('y in m'); zlabel('z in m');
  view(3);
  RP.plot(qs, X0, s_plot);
  title('Startpose mit Kollisionsmodell des Roboters');
end
[colldet_start, colldist_start] = check_collisionset_simplegeom_mex(RP.collbodies, ...
  RP.collchecks, Tc_s(:,4)', struct('collsearch', true));
assert(all(~colldet_start(:)), ['Es sollte keine Kollision in der Startpose ', ...
  'vorliegen (manuell geprüft)']);

%% Initialisierung Teil 2
% Roboter auf 3T2R einstellen
RP.update_EE_FG(I_EE_full, I_EE_red);
save(fullfile(data_path, 'robot_definition.mat'), 'RP');

%% Eckpunkte für Beispiel-Trajektorie bestimmen und IK prüfen
% IK-Grundeinstellungen
s = struct('Phit_tol', 1e-12, 'Phir_tol', 1e-12, ... % sehr genau rechnen
  'maxstep_ns', 1e-5, ... % Schrittweite für Nullraum-Inkremente gering halten
  'scale_lim', 0.7, ...
  'retry_limit', 0);
% Rechteck-Trajektorie
if usr_trajnum == 1
  d1=0.25; d2=0.2; % dimensions of the rectangle
  % First corner without tilting
  X1 = X0' + [-d1/5,+d2/5,+0.2,0,0,0];
  ta = pi/4; % tilting angle, point from outside on edge of rectangle
%   r_P_E = zeros(3,1);
elseif usr_trajnum == 2
  d1 = 0.4; d2=0.3; % dimensions of the rectangle
  X1 = X0' + [-d1/2,+d2/2,+0.2,0,0,0];
  ta = -pi/4;%pi/3;
%   r_P_E = [0;0;50e-3]; % 50mm langer Fräser inkl. Spindel
else
  error('Fall nicht definiert');
end
% RP.update_EE(r_P_E);
XL = [];
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
% Entferne z-Rotation
XL(:,6) = 0;
% Entferne doppelte, direkt aufeinanderfolgende Einträge
last_i = 1;
for i = 2:size(XL,1)
  if abs(XL(last_i,:)-XL(i,:)) < 1e-10
    XL(i,:) = NaN;
    continue
  end
  last_i = i;
end
XL = XL(all(~isnan(XL),2),:);

% Determine tilting angles
TA = NaN(size(XL,1), 1);
for i = 1:size(XL,1)
  R_i = eulxyz2r(XL(i,4:6)');
  TA(i) = acos(R_i(:,3)'*[0;0;1]);
end
% Debug:
% 180/pi*TA

%% Generate Trajectory
% Use moderate end effector velocity and acceleration to have reserves for
% nullspace motion
vmax = 50e-3; % 50mm/s
phiDmax = 10*pi/180; % 10°/s

[Y_t, YD_t, YDD_t, t, IL] = traj_trapez2_multipoint(XL(:,1:5), [vmax*[1;1;1];phiDmax*[1;1]]', ...
  0.01, 0.01, 1e-3, 0);
amax = vmax/0.01;
X_t = [Y_t, zeros(length(t),1)];
XD_t = [YD_t, zeros(length(t),1)];
XDD_t = [YDD_t, zeros(length(t),1)];
% Velocity profile to obtain rest-to-rest motion also in the nullspce
nullspace_maxvel_interp = nullspace_maxvel_from_tasktraj(t, ...
  IL, vmax/amax, 0.005 , 1e-3);

if usr_short_traj
  n = 2000;
else
  n = length(t);
end
IL = IL(IL<=length(t));
II = 1:n;
t = t(II);
X_t = X_t(II,:);
XD_t = XD_t(II,:);
XDD_t = XDD_t(II,:);
fprintf('Length of the trajectory: %1.1fmm, %1.1fdeg\n', ...
  1e3*sum(sum(abs(diff(X_t(:,1:3))))), 180/pi*sum(sum(abs(diff(X_t(:,4:5))))))
fprintf('Duration of the trajectory: %1.1fs\n', t(end));

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
    saveas(100, fullfile(data_path, sprintf('robot_with_traj%d.fig',usr_trajnum)));
    saveas(100, fullfile(data_path, sprintf('robot_with_traj%d.png',usr_trajnum)));
  end
end

%% Daten für die Tabellen
% Give table of support points of the trajectory
fprintf('Table of support points:\n');
for i = 1:size(XL,1)
  fprintf('%d & %1.0f & %1.0f & %1.0f & %1.0f & %1.0f\\\\\n', i-1, ...
    1e3*XL(i,1), 1e3*XL(i,2), 1e3*XL(i,3), 180/pi*XL(i,4), 180/pi*XL(i,5));
end

fprintf('Robot Data (Table I in paper):\n');
fprintf('base diameter: %1.1fmm\n', 2*RP.DesPar.base_par(1)*1e3);
fprintf('platform diameter: %1.1fmm\n', 2*RP.DesPar.platform_par(1)*1e3);
fprintf('platform pair joint distance: %1.1fmm\n', RP.DesPar.platform_par(2)*1e3);
fprintf('Stroke limits: %1.1f ... %1.1f mm\n', RP.Leg(1).qlim(3,1)*1e3, ...
  RP.Leg(1).qlim(3,2)*1e3);
fprintf('Diameter of collision bodies (legs): %1.1fmm\n', 2*RP.collbodies.params(end,1)*1e3);
fprintf('Diameter of collision bodies (platform): %1.1fmm\n', 2*RP.collbodies.params(1,1)*1e3);

%% Berechne IK zu den einzelnen Eckpunkten (zum Testen)
% Berechne IK zu den einzelnen Eckpunkten. Wenn das nicht geht, bringt die
% Trajektorie sowieso nichts. Benutze die IK mit Aufgabenredundanz 
s_ep = s; % Einstellungen für die Eckpunkte-IK
s_ep.wn = zeros(RP.idx_ik_length.wnpos, 1);
% s_ep.wn(RP.idx_ikpos_wn.qlim_par) = 1; % Nehme Zielfunktion 1. Damit Überschreitung der Ränder in Zwischenständen möglich
s_ep.wn(RP.idx_ikpos_wn.jac_cond) = 1; % Nehme Zielfunktion 1. Damit Überschreitung der Ränder in Zwischenständen möglich
s_ep.n_max = 5000; % Mehr Versuche (Abstände zwischen Punkten größer als bei Traj.-IK)
% s_ep.maxrelstep_ns = 0.05; % Große Werte, Größere Nullraumbewegung pro Zeitschritt
s_ep.retry_on_limitviol = true;
s_ep.retry_limit = 20; % Neuversuche erlauben (bei Einzelpunkt i.O.)
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
  XL(i,6) = x_i(6);
  [~,Phi_i_voll] = RP.constr1(q_i, x_i(:));
  assert(all(abs(Phi_i_voll)<1e-8), ...
    sprintf('Ergebnis der %s-IK (Parallel) ist falsch', I_EE_red_str));
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
[Y_tref,~,~,~,ILref] = traj_trapez2_multipoint(XL(:,1:5), 3, 0.01, 0.01, 1e-4, 0);
X_tref = [Y_tref, zeros(size(Y_tref,1),1)];
% Einstellungen vorbereiten
settings_perfmap = struct('settings_ik', s, ...
  'q0', QL(1,:)', 'I_EE_red', I_EE_red, 'map_phistart', XL(1,6), ...
  'maplim_phi', [-pi,pi]);
if usr_highres_distrfig % settings for high resulution of performance map
  % three values for horizontal resolution of the performance map
  settings_perfmap.mapres_thresh_eepos = 1e-3; % 1mm
  settings_perfmap.mapres_thresh_eerot = 3*pi/180; % 1deg
  settings_perfmap.mapres_thresh_pathcoordres = 0.01; % min 100 values between two trajectory key points
  % vertical resolution of the performance map. Resolution for the
  % redundant coordinate phi_z
  settings_perfmap.mapres_redcoord_dist_deg = 0.5; % deg
  highresstr = 'highres';
else % settings for low resulution of map
  settings_perfmap.mapres_thresh_eepos = 3e-3; % 3mm
  settings_perfmap.mapres_thresh_eerot = 3*pi/180; % 3deg
  settings_perfmap.mapres_thresh_pathcoordres = 0.05;% min 20 values between key points
  settings_perfmap.mapres_redcoord_dist_deg = 5; % deg
  highresstr = 'lowres';
end
% Load previously calculated performance map
filename_pre = sprintf('LNEE_pkm_traj%d_%dsamples', usr_trajnum, size(X_tref,1));
filename_perfmap = fullfile(data_path, [filename_pre, '_perfmap_', highresstr, '.mat']);
data_loaded_offline = false;
if usr_load_discretization && ~exist(filename_perfmap, 'file')
  warning('Unable to load discretization results from %s', filename_perfmap);
elseif usr_load_discretization && exist(filename_perfmap, 'file')
  d = load(filename_perfmap);
  if all(size(XL)==size(d.XL))
    test_XL = XL - d.XL;
  else
    test_XL = 1; % raises error
  end
  if ~isfield(d, 'collchecks') || size(d.collchecks,1)~=size(RP.collchecks,1) || ...
      any(any(d.collchecks-RP.collchecks))
    test_XL = 1; % raises error: Collision checks are different
  end
  if ~isfield(d, 'collbodies_params') || size(d.collbodies_params,1)~=size(RP.collbodies.params,1) || ...
      any(any(abs(d.collbodies_params-RP.collbodies.params) > 1e-10))
    test_XL = 1; % raises error: Collision checks are different
  end
  for f = fields(settings_perfmap)'
    if isa(settings_perfmap.(f{1}), 'double')
      delta_value = settings_perfmap.(f{1}) - d.settings_perfmap.(f{1});
      if any(abs(delta_value)>1e-6)
        warning('Setting %s is different in file. Delta=[%s]', f{1}, ...
          disp_array(delta_value(:)','%1.1e'));
        test_XL = 1; % raises error
        break;
      end
    end
  end
  if any(abs(test_XL(:)) > 1e-5)
    % Durch die Abhängigkeit der Werte von der IK oben sind die Ergebnisse
    % immer leicht anders
    warning('Data from file does not match the settings');
  else
    fprintf('Successfully loaded performance map from file\n');
    data_loaded_offline = true;
    % Load results
    H_all = d.H_all;
    Q_all = d.Q_all;
    phiz_range = d.phiz_range;
    s_ref = d.s_ref;
    s_tref = d.s_tref;
  end
end
if ~data_loaded_offline
  fprintf('Start discretization of the inverse kinematics\n');
  [H_all, Q_all, s_ref, s_tref, phiz_range] = RP.perfmap_taskred_ik( ...
    X_tref, ILref, settings_perfmap);
  fprintf('Finished discretization of the trajectory. Total %1.1fmin.\n', toc(t1)/60);
  collchecks = RP.collchecks;
  collbodies_params = RP.collbodies.params;
  save(filename_perfmap, 'H_all', 'Q_all', 's_ref', 's_tref', ...
    'phiz_range', 'XL', 'settings_perfmap', 'collchecks', 'collbodies_params');
end
if length(unique(phiz_range))~=length(phiz_range)
  error('Something went wrong when assembling phiz_range');
end
fprintf(['The performance map contains %d trajectory samples and %d ', ...
  'values for the redundant coordinates. %d evaluations in total\n'], ...
  length(s_ref), length(phiz_range), length(s_ref)*length(phiz_range));

% Funktion nochmal aufrufen und spätere Trajektorie gegen aktuelle
% abgleichen über die Bahnkoordinate
settings_perfmap.only_trajectory_normalization = true;
[~, ~, s_ref2, s_tref2] = RP.perfmap_taskred_ik( ...
  X_t, IL, settings_perfmap);
PM_s = NaN(length(s_ref),1);
PM_st = s_tref2;
s_err_max = 0;
for i = 1:length(s_ref)
  [s_err, I_best] = min(abs(s_ref(i)-s_ref2));
  PM_s(i) = s_ref2(I_best);
  s_err_max = max(s_err_max, s_err);
end
fprintf(['Redundanzkarte umgerechnet auf aktuelle Trajektorie. Größter ', ...
  'Fehler der Bahnkoordinate dabei %1.3f\n'], s_err_max);

%% Plot global distribution
% Kriterium für Singularität
abort_thresh_hpose = inf(RP.idx_ik_length.hnpos, 1);
abort_thresh_hpose(RP.idx_ikpos_hn.jac_cond) = 1e3;
abort_thresh_hpose(RP.idx_ikpos_hn.ikjac_cond) = 1e4;
if debug_plot && false % Distribution of PKM Jacobian condition number
  change_current_figure(3); clf; hold on;
  set(3, 'Name', 'Distr_PKMJacCrit', 'NumberTitle', 'off');
  settings_perfmapplot = struct( ...
    'markermindist', [0.20, 15], ... % nur alle s=... und phi=...° einen Marker
    'wn', zeros(RP.idx_ik_length.wnpos,1));
  settings_perfmapplot.wn(RP.idx_ikpos_wn.jac_cond) = 1;
  RP.perfmap_plot(H_all, phiz_range, s_ref, settings_perfmapplot); 
  if usr_save_figures
    saveas(3, fullfile(data_path, [filename_pre,'_pkmjac_condition.fig']));
    saveas(3, fullfile(data_path, [filename_pre,'_pkmjac_condition.png']));
  end
end
if debug_plot && false % Distribution of inverse kinematics Jacobian condition number
  change_current_figure(6); clf; hold on;
  set(6, 'Name', 'Distr_IKJacCrit', 'NumberTitle', 'off');
  settings_perfmapplot = struct( ...
    'markermindist', [0.20, 15], ...
    'wn', zeros(RP.idx_ik_length.wnpos,1));
  settings_perfmapplot.wn(RP.idx_ikpos_wn.ikjac_cond) = 1;
  RP.perfmap_plot(H_all, phiz_range, s_ref, settings_perfmapplot); 
  title(sprintf('IK Jacobian Condition Number Criterion'));
  if usr_save_figures
    saveas(6, fullfile(data_path, [filename_pre,'_ikjac_condition.fig']));
    saveas(6, fullfile(data_path, [filename_pre,'_ikjac_condition.png']));
  end
end
if debug_plot && false % Distribution of joint limits criterion
  PMjlfig = change_current_figure(5);clf; hold on;
  set(PMjlfig, 'Name', 'Distr_JointLimCrit', 'NumberTitle', 'off');
  settings_perfmapplot = struct( ...
    'markermindist', [0.20, 15], ...
    'wn', zeros(RP.idx_ik_length.wnpos,1));
  settings_perfmapplot.wn(RP.idx_ikpos_wn.qlim_hyp) = 1;
  RP.perfmap_plot(H_all, phiz_range, s_ref, settings_perfmapplot); 
  title('Joint Limit Criterion (1000=invalid)');
  if usr_save_figures
    saveas(PMjlfig, fullfile(data_path, [filename_pre,'_limitviol_crit.fig']));
    saveas(PMjlfig, fullfile(data_path, [filename_pre,'_limitviol_crit.png']));
  end
end
if debug_plot && false % Distribution of collision criterion
  PMcollfig = change_current_figure(9);clf; hold on;
  set(PMjlfig, 'Name', 'Distr_CollCrit', 'NumberTitle', 'off');
  settings_perfmapplot = struct( ...
    'markermindist', [0.20, 15], ...
    'wn', zeros(RP.idx_ik_length.wnpos,1));
  settings_perfmapplot.wn(RP.idx_ikpos_wn.coll_hyp) = 1;
  RP.perfmap_plot(H_all, phiz_range, s_ref, settings_perfmapplot); 
  title('Collision Criterion');
  if usr_save_figures
    saveas(PMcollfig, fullfile(data_path, [filename_pre,'_collision_crit.fig']));
    saveas(PMcollfig, fullfile(data_path, [filename_pre,'_collision_crit.png']));
  end
end
if debug_plot % Distribution of combined criterion (condition+limits)
  PMcombfig = change_current_figure(8);clf; hold on;
  set(PMcombfig, 'Name', 'Distr_CombCrit', 'NumberTitle', 'off');
  settings_perfmapplot = struct( ...
    'markermindist', [0.20, 15], ...
    'wn', zeros(RP.idx_ik_length.wnpos,1), ...
    'abort_thresh_h', abort_thresh_hpose);
  settings_perfmapplot.wn(RP.idx_ikpos_wn.jac_cond) = 1;
  settings_perfmapplot.wn(RP.idx_ikpos_wn.ikjac_cond) = 1;
  settings_perfmapplot.wn(RP.idx_ikpos_wn.qlim_hyp) = 1;
  settings_perfmapplot.wn(RP.idx_ikpos_wn.coll_hyp) = 1;
  [Hdl_all, s_pmp] = RP.perfmap_plot(H_all, phiz_range, s_ref, settings_perfmapplot); 
  title('Combined Criterion');
  if usr_save_figures
    saveas(PMcombfig, fullfile(data_path, [filename_pre,'_combined_crit.fig']));
    saveas(PMcombfig, fullfile(data_path, [filename_pre,'_combined_crit.png']));
  end
end

if false && debug_plot % For Debugging the global discretization
  i_traj = 1;
  [~,i_phi] = min(abs(phiz_range-0)); %#ok<ASGLU>
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
    PM_s(i_traj), i_traj));
  linkxaxes
end
%% Create paper figure with robot in initial pose of the trajectory
% This creates Fig. 1 and 12 of the LNEE paper
% Figures are updated manually for trajectory 2
if usr_plot_robot
RobPoseFig = change_current_figure(900);clf;
set(RobPoseFig,'Name','Rob','NumberTitle','off');
view(3);
% Plot arbitrary configuration to check if the platform rotation worked
if usr_trajnum == 1 % Fig. 1b of the LNEE paper
  i_traj = 1; % beginning of the trajectory
  [~,i_phi] = min(abs(phiz_range-0)); % value for phiz=0
else
  % Singularity: Has to be executed manually; Fig. 12a of the LNEE paper
%   [~, i_traj] = min(abs(PM_s-1.132)); % arbitrary value
%   [~,i_phi] = min(abs(phiz_range-36.268*pi/180)); % value for phiz=0
  % Collision: Has to be executed manually; Fig. 12c of the LNEE paper
%   [~, i_traj] = min(abs(PM_s-4.72)); % arbitrary value
%   [~,i_phi] = min(abs(phiz_range-(-70)*pi/180)); % value for phiz=0
  % Normal Case; Fig. 12b of the LNEE paper
  [~, i_traj] = min(abs(PM_s-3)); % arbitrary value
  [~,i_phi] = min(abs(phiz_range-0)); % value for phiz=0
end
% Gelenkwinkel für gesuchte Bahnkoordinate aus Redundanzkarte
q_test = Q_all(i_phi,:,i_traj)';
% EE-Trajektorie anhand der Bahnkoordinate aus Traj.-Variable bestimmen
[~, i_traj2] = min(abs(PM_st-PM_s(i_traj)));
% Test the loaded configuration
x_test = [Y_t(i_traj2,1:5), phiz_range(i_phi)]';
[~,Phi_test] = RP.constr1(q_test, x_test);
% Test nur mit grober Toleranz, da bei Umrechnung von Redundanzkarte und
% Trajektorie oben ein Fehler gemacht wird.
assert(all(abs(Phi_test)<1e-2), 'Kinematic constraints in loaded config do not match');
% title(sprintf('Robot in Selected Pose: traj idx %d, phiz=%1.0fdeg', ...
%   i_traj, 180/pi*x_test(6)));
hold on;grid on;
s_plot = struct(  'ks_legs', [], 'straight', 0, 'mode', 4);
% Modify the robot pose to avoid overlapping of the platform and the
% trajectory. Only minor change for visibility.
x_paperfig = x_test;
x_paperfig(1:3) = x_paperfig(1:3) + eulxyz2r(x_paperfig(4:6))*[0;0;-20e-3];
[q_paperfig, Phi_paperfig] = RP.invkin_ser(x_paperfig, q_test);
assert(all(abs(Phi_paperfig)<1e-6), 'IK for robot paper figure not successful');
% Eigenschaften nennen für Verwendung in Paper
filename_savepre = sprintf('pkm_traj%d_s%1.0f_pose_phi%1.0fdeg', usr_trajnum, PM_s(i_traj), 180/pi*phiz_range(i_phi));
if usr_save_figures
  fid = fopen(fullfile(paperfig_path, [filename_savepre, '.txt']), 'w');
  fprintf(fid, 'Trajectory %d:\n', usr_trajnum);
  fprintf(fid, 'Robot at s=%1.1f, phiz=%1.1f deg:\n', PM_s(i_traj), 180/pi*phiz_range(i_phi));
  fprintf(fid, 'cond(J) = %1.1f\n', H_all(i_traj, i_phi, end)); % Letzter Eintrag ist Kondition
  fprintf(fid, 'h_coll_hyp = %1.1f\n', H_all(i_traj, i_phi, RP.idx_ikpos_hn.coll_hyp));
  fprintf(fid, 'd_coll = %1.1fmm\n', 1e3*H_all(i_traj, i_phi, end-3));
  fclose(fid);
end
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
if usr_trajnum == 1
  view([140, 20]);
else
  view([-88, 20]); % TODO: -88?
end
set(RobPoseFig, 'windowstyle', 'normal');
set_size_plot_subplot(RobPoseFig, ...
  10,10,gca,...
  0,0,0,0,0,0)
drawnow();
% Bild speichern
if usr_save_figures
  exportgraphics(RobPoseFig, fullfile(paperfig_path, [filename_savepre,'.png']),'Resolution',600)
  saveas(RobPoseFig, fullfile(data_path, [filename_savepre,'.fig']));
end
if usr_only_global_discretization
  return
end
end
%% Inverse Kinematik zum Startpunkt der Trajektorie
% Inverse Kinematik berechnen;  Lösung der IK von oben als Startwert
% Wird nicht für die Trajektorien-IK benutzt, da die optimale Startkon-
% figuration von den benutzten Nebenbedingungen abhängt.
tic();
s_start = s;
s_start.wn = zeros(RP.idx_ik_length.wnpos, 1);
s_start.wn(RP.idx_ikpos_wn.qlim_par) = 1;
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
% Debug:
% Nehme ersten Punkt von oben. Dadurch genau gleiche Bedingungen wie für
% Redundanzkarte
% q1 = QL(1,:)';
% % Berechne Ersten Punkt der Trajektorie mit Aufgabenredundanz.
% % Dadurch bestmögliche Startkonfiguration
[q1, Psi_num1, ~, Stats1] = RP.invkin3(XL(1,:)', QL(1,:)', s_start);
if any(abs(Psi_num1) > 1e-4)
  error('IK konvergiert nicht für Startpunkt der Trajektorie');
end
x1 = RP.fkineEE_traj(q1')';
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

%% Compute trajektory with dynamic programming
for i_dpred = 1:5 % Different settings for DP
  overlap = false;
  freetransfer = false;
  stageopt_posik = false;
  if i_dpred == 1 % from Sect. 4.1
    % no redundancy in dynamic programming
    n_phi = 1+360/45; % wide discretization like SI-DP: 45°
    RP.update_EE_FG(I_EE_full, I_EE_full);
    suffix_red = 'nored';
  elseif i_dpred == 2 % from Sect. 4.1
    % no redundancy in dynamic programming
    n_phi = 1+360/15; % fine discretization necessary: 15°
    RP.update_EE_FG(I_EE_full, I_EE_full);
    suffix_red = 'nored';
  elseif i_dpred == 3 % from Sect. 4.2; with intervals
    n_phi = 1+360/45; % nur grobe Diskretisierung, sogar vorteilhaft
    RP.update_EE_FG(I_EE_full, I_EE_red);
    suffix_red = 'red';
  elseif i_dpred == 4 % from Sect. 4.3; with overlapping intervals
    n_phi = 1+360/60; % gröbere Diskretisierung (wegen Überlapp)
    RP.update_EE_FG(I_EE_full, I_EE_red);
    suffix_red = 'red';
    overlap = true;
    freetransfer = false; % TODO: In Auswertung diskutieren?
  elseif i_dpred == 5 % Additional evaluation for PhD thesis (stage optimization)
    n_phi = 1+360/45; % nur grobe Diskretisierung
    RP.update_EE_FG(I_EE_full, I_EE_full);
    suffix_red = 'nored'; % Nutze Redundanz nur mit Positions-IK auf Stufe
    stageopt_posik = true;
  else
    error('Fall nicht definiert');
  end
  fprintf('DP with setting "%s", overlap=%d and max. %d states:\n', suffix_red, overlap, n_phi);
  % String für Namen der zu speichernden Dateien
  wn_traj_default = zeros(RP.idx_ik_length.wntraj,1);
  wn_traj_default(RP.idx_iktraj_wnP.qlim_hyp) = 1; % K_P (hyperb. limit)
  wn_traj_default(RP.idx_iktraj_wnD.qlim_hyp) = 0.6; % K_D (limit)
  wn_traj_default(RP.idx_iktraj_wnP.qDlim_par) = 0.3; % K_v
  wn_traj_default(RP.idx_iktraj_wnP.jac_cond) = 1; % K_P (cond)
  wn_traj_default(RP.idx_iktraj_wnD.jac_cond) = 0.3; % K_D (cond)
  wn_traj_default(RP.idx_iktraj_wnP.ikjac_cond) = 1; % K_P (cond)
  wn_traj_default(RP.idx_iktraj_wnD.ikjac_cond) = 0.3; % K_D (cond)
  wn_traj_default(RP.idx_iktraj_wnP.coll_hyp) = 1; % K_P (hyperb. limit)
  wn_traj_default(RP.idx_iktraj_wnD.coll_hyp) = 0.6; % K_D (limit)
  % Umrechnen auf Positions-IK
  wnpos_dp = zeros(RP.idx_ik_length.wnpos,1);
  for f = fields(RP.idx_ikpos_wn)'
    wnpos_dp(RP.idx_ikpos_wn.(f{1})) = wn_traj_default(RP.idx_iktraj_wnP.(f{1}));
  end
  % Konditionszahl permanent optimieren. IK-Kondition nur, wenn sehr schlecht
  s_Traj = struct('cond_thresh_ikjac', 500, 'cond_thresh_jac', 1, ... 
    'thresh_ns_qa', 1); % immer Nullraumbewegung in Antriebskoordinaten. Vereinfacht das Paper
  if isinf(s_Traj.thresh_ns_qa)
    suffix_red = [suffix_red, '_nsactcoord']; %#ok<AGROW> 
  end
  % Abbruch, sobald die Jacobi-Matrix sehr schlecht wird oder andere Neben-
  % bedingungen verletzt werden
  abort_thresh_h = inf(RP.idx_ik_length.hntraj, 1);
  abort_thresh_h(RP.idx_iktraj_hn.jac_cond) = 1e3;
  abort_thresh_h(RP.idx_iktraj_hn.ikjac_cond) = 1e4;
  s_Traj.abort_thresh_h = abort_thresh_h;
  DP_settings = struct('phi_min', -pi, 'phi_max', pi, 'n_phi', n_phi, ...
    'abort_thresh_h', abort_thresh_h, ...
    'PM_H_all', H_all, 'PM_s_ref', PM_s, 'PM_s_tref', PM_st, 'PM_limit', true, ...
    'PM_phiz_range', phiz_range, 'IE', IL, 'PM_Q_all', Q_all, ... % 
    'wn', wnpos_dp, ...
    'use_free_stage_transfer', freetransfer, ...
    'overlap', overlap, ... % Überlappende Intervalle
    'stageopt_posik', stageopt_posik, ... % Nach-Optimierung mit Positions-IK
    'debug', true, ...
    'fastdebug', false, ... % true: damit Verzicht auf Prüfung und schnelleres Zeichnen
    'phi_lim_x0_dependent', false, ... % true: Grenzen abhängig von phi0 wählen
    'continue_saved_state', false, ...
    'settings_ik', s_Traj, ...
    'cost_mode', 'RMStraj', ... % Alternative: 'max', 'RMStime', 'average'
    'verbose', 2);
  if DP_settings.use_free_stage_transfer && i_dpred == strcmp(suffix_red, 'red')
    fststr = '_freetransfer';
  else         
    fststr = '';
  end
  if DP_settings.overlap && i_dpred == strcmp(suffix_red, 'red')
    fststr = [fststr, '_overlap']; %#ok<AGROW> 
  end
  suffix = ['cost', DP_settings.cost_mode, '_n', num2str(n_phi), '_', ...
    suffix_red, fststr];
  if DP_settings.use_free_stage_transfer
    suffix = [suffix, '_freetransfer']; %#ok<AGROW> 
  end
  if overlap
    suffix = [suffix, '_overlap']; %#ok<AGROW> 
  end
  if stageopt_posik
    suffix = [suffix, '_stageopt']; %#ok<AGROW> 
  end
  if ~DP_settings.phi_lim_x0_dependent
    suffix = [suffix, '_phi0fix']; %#ok<AGROW> 
  end
  DP_settings.debug_dir = fullfile(respath, sprintf('LNEE_Traj%d_DP_debug_%s', usr_trajnum, suffix));
  filename_dynprog= fullfile(data_path, [filename_pre, '_dynprog_', suffix, '.mat']);
  dynprog_loaded_offline = false;
  if usr_load_dynprog && (~exist(filename_dynprog, 'file') || ...
      ~exist(DP_settings.debug_dir, 'file'))
    fprintf('Unable to load discretization results from %s\n', filename_dynprog);
  elseif usr_load_dynprog && exist(filename_dynprog, 'file')
    d = load(filename_dynprog);
    % Test if settings from loaded results match
    if size(X_t,1) ~= size(d.X_t,1)
      warning('Loaded file %s does not match. Disregard', filename_dynprog);
    elseif any(DP_settings.wn ~= d.DP_settings.wn) || ...
        DP_settings.n_phi ~= d.DP_settings.n_phi
      warning('IK settings from loaded results are different. Disregard');
    elseif ~isfield(d, 'collchecks') || size(d.collchecks,1)~=size(RP.collchecks,1) || ...
        any(any(d.collchecks-RP.collchecks))
      warning('Collision checks loaded from results do not match. Disregard.');
    elseif ~isfield(d, 'collbodies_params') || size(d.collbodies_params,1)~=size(RP.collbodies.params,1) || ...
        any(any(abs(d.collbodies_params-RP.collbodies.params) > 1e-10))
      warning('Collision bodes loaded from results do not match. Disregard.');
    elseif ~isfield(d.DP_Stats, 'phi_range')
      warning('DP output DP_Stats has no field phi_range. Disregard.');
    else
      DP_XE = d.DP_XE;
      DP_Stats = d.DP_Stats;
      DP_TrajDetail = d.DP_TrajDetail;
      dynprog_loaded_offline = true;
    end
  end
  mkdirs(DP_settings.debug_dir);
  if ~exist(DP_settings.debug_dir, 'file')
    error('Debug-Ordner existiert nicht. Fehler mit Sym-Link?');
  end
  if ~dynprog_loaded_offline
    t1 = tic();
    fprintf('Start computation for dynamic programming of trajectory\n');
    [DP_XE, DP_Stats, DP_TrajDetail] = RP.dynprog_taskred_ik(X_t, XD_t, ...
      XDD_t, t, q1, DP_settings);
    fprintf('Dynamic programming computed in %1.1f min\n', toc(t1)/60);
    collchecks = RP.collchecks;
    collbodies_params = RP.collbodies.params;
    save(filename_dynprog, 'DP_XE', 'DP_Stats', 'DP_TrajDetail', ...
      'DP_settings', 'X_t', 'collchecks', 'collbodies_params');
  end
  fprintf(['Computation time for the DP trajectory: %1.1f min for %d ', ...
    'trajectory samples\n'], sum(DP_Stats.Stats_comptime_ik_all(:))/60, sum(DP_Stats.nt_ik_all(:)))
  ts_avg = sum(DP_Stats.Stats_comptime_ik_all(:))/sum(DP_Stats.nt_ik_all(:));
  fprintf('Average time per sample: %1.1fms\n', ts_avg*1e3);
  fprintf('Computed samples correspond to %1.2f x trajectory\n', ...
    sum(DP_Stats.nt_ik_all(:))/size(X_t,1));
  fprintf('Expected maximum %d x trajectory (n²)\n', n_phi^2);
  % Belege die Variablen mit den Ergebnissen
  if i_dpred == 1
    DP_TrajDetail_nored = DP_TrajDetail;
    DP_XE_nored = DP_XE;
  elseif i_dpred == 2
    DP_TrajDetail_nored_fine = DP_TrajDetail;
    DP_XE_nored_fine = DP_XE;
  elseif i_dpred == 3
    DP_TrajDetail_red = DP_TrajDetail;
    DP_XE_red = DP_XE;
  elseif i_dpred == 4
    DP_TrajDetail_redol = DP_TrajDetail;
    DP_XE_redol = DP_XE;
  elseif i_dpred == 5
    DP_TrajDetail_redso = DP_TrajDetail;
    DP_XE_redso = DP_XE;
  else
    error('Fall nicht definiert');
  end
end

%% IK für Trajektorie berechnen (Vorbereitung)

RP.update_EE_FG(I_EE_full, I_EE_red);
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

filename_traj = fullfile(data_path, filename_pre);
if ~DP_settings.phi_lim_x0_dependent
  filename_traj = [filename_traj, '_phi0fix'];
end
filename_traj = [filename_traj, '_traj.mat'];
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
  if length(d.Namen_Methoden) ~= length(Namen_Methoden)
    warning('Only %d trajectories in file, but %d expected.', length( ...
      d.Namen_Methoden), length(Namen_Methoden));
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
  Hcond_all = d.Hcond_all;
  Namen_Methoden = d.Namen_Methoden;
  Namen_Methoden_Leg_Paper = d.Namen_Methoden_Leg_Paper;
  fprintf('Loaded results for trajectory without calculation\n');
end

for kk = 1:length(Namen_Methoden)*(~usr_load_traj)
  set_kk = s_Traj;
  set_kk.debug = true; % check for errors in trajectory IK function
  set_kk.thresh_ns_qa = 1; % always use full joint projector
  for j = 1:RP.NLEG
    RP.Leg(j).qlim = qlim_backup(RP.I1J_LEG(j):RP.I2J_LEG(j),:);
    RP.Leg(j).qDlim = qDlim_backup(RP.I1J_LEG(j):RP.I2J_LEG(j),:);
    RP.Leg(j).qDDlim = qDDlim_backup(RP.I1J_LEG(j):RP.I2J_LEG(j),:);
  end
  phiz_0 = x1(6); % use optimal start configuration from above
  phiz_const = 0; % is ignored in 3T2R IK
  wn_traj = wn_traj_default;
  I_EE_Task_kk = I_EE_red; % use task redundancy by default
  % use resting velocity profile for nullspace motion
  set_kk.nullspace_maxvel_interp = nullspace_maxvel_interp;
  % Set tolerance band for redundant coordinate from dynamic programming
  use_dp_no_recalc = false;
  switch kk
    case 1
      % Ergebnis der dynamischen Programmierung ohne Redundanz.
      % Rast-zu-Rast-Bewegung.
      name_method=sprintf('DP_rest_nored45');
      name_method_leg = 'DP 45°';
      use_dp_no_recalc = true;
      DP_TrajDetail_kk = DP_TrajDetail_nored;
      set_kk.nullspace_maxvel_interp = zeros(2,0);
      set_kk.wn(:) = 0; % Keine zusätzliche Optimierung möglich
    case 2
      % Ergebnis der dynamischen Programmierung ohne Redundanz.
      % Rast-zu-Rast-Bewegung.
      name_method=sprintf('DP_rest_nored15');
      name_method_leg = 'DP 15°';
      use_dp_no_recalc = true;
      DP_TrajDetail_kk = DP_TrajDetail_nored_fine;
      set_kk.nullspace_maxvel_interp = zeros(2,0);
      set_kk.wn(:) = 0; % Keine zusätzliche Optimierung möglich
    case 3
      % Ergebnis der dynamischen Programmierung ohne Redundanz.
      % Rast-zu-Rast-Bewegung.
      name_method=sprintf('DP_rest_red');
      name_method_leg = 'SI 45°';
      use_dp_no_recalc = true;
      DP_TrajDetail_kk = DP_TrajDetail_red;
    case 4
      % Ergebnis der dynamischen Programmierung ohne Redundanz.
      % Rast-zu-Rast-Bewegung.
      name_method=sprintf('DP_rest_redol');
      name_method_leg = 'OL 90°';
      use_dp_no_recalc = true;
      DP_TrajDetail_kk = DP_TrajDetail_redol;
    case 5
      % Lokale Optimierung ohne Vorwissen
      % Rast-zu-Rast Bewegung (auch im Nullraum)
      name_method=sprintf('LocalOpt_rest');
      name_method_leg = 'NP';
    otherwise
      error('Fall %d noch nicht definiert', kk);
  end
  if isnan(phiz_0)
    error('phiz_0 is not set');
  end
  fprintf('Berechne Trajektorie %d (%s)\n', kk, name_method);
  set_kk.wn = wn_traj;
  % Positions-IK zum Startpunkt der Trajektorie mit genau den gleichen
  % Optimierungs-Nebenbedingungen wie in der Trajektorie. Dadurch keine
  % Nullraumbewegung am Anfang (Start in lokalem Optimum)
  s_pik_kk = struct('scale_lim',0.5); % Grenzen dürfen nicht überschritten werden
  % s_pik_kk.wn = s_kk.wn([1 2 5]); % Positions-Grenzen und Kondition
  % Wähle immer die gleichen Nebenbedingungen, damit alle mit gleicher
  % Konfiguration starten (besser vergleichbar)
  qs_kk = q1;
  RP.update_EE_FG(I_EE_full, I_EE_full); % Für Berechnung der Zwangsbed.
  X_t_kk = X_t; XD_t_kk = XD_t; XDD_t_kk = XDD_t;
  X_t_kk(:,6) = phiz_const; XD_t_kk(:,6) = 0; XDD_t_kk(:,6) = 0;
  assert(all(abs(RP.constr3(qs_kk,[X_t(1,1:5),phiz_0]'))<1e-10), ...
    'initial joint configuration does not match');
  qskknorm = (qs_kk-qlim(:,1)) ./ (qlim(:,2) - qlim(:,1));
  if any(qs_kk<qlim(:,1) | qs_kk>qlim(:,2))
    I_viol = find(qs_kk<qlim(:,1) | qs_kk>qlim(:,2));
    error('Startpunkt liegt außerhalb der Grenzen. Verletzt: [%s]', disp_array(I_viol','%d'));
  end
  xs_kk_test = RP.fkineEE_traj(qs_kk')';
  assert(abs(xs_kk_test(6)-phiz_0)< 1e-10, ...
    'redundant coordinate in starting pose is not as requested');
  RP.update_EE_FG(I_EE_full, I_EE_Task_kk);
  Namen_Methoden{kk} = name_method;
  Namen_Methoden_Leg_Paper{kk} = name_method_leg;
  assert(all(~isnan(X_t_kk(:))), 'Trajectory contains NaN. Syntax error');
  if ~use_dp_no_recalc
    % Grobe Abschätzung der Rechenzeit
    comptime_tmp = NaN(3,1);
    for jj = 1:3 % Erste Iterationen nur, um Funktion in Speicher zu laden, 2. für Differenzbildung
      t1 = tic();
      RP.invkin2_traj(X_t_kk(1:10*jj,:), XD_t_kk(1:10*jj,:), XDD_t_kk(1:10*jj,:), ...
        t(1:10*jj,:), qs_kk, set_kk);
      comptime_tmp(jj) = toc(t1);
    end
    fprintf(['Tested Trajectory execution for 10 samples. Avg %1.1fms per ', ...
      'sample. Estimated %1.1fmin for full trajectory of %d samples.\n'], ...
      1e3*diff(comptime_tmp([2;3]))/10, diff(comptime_tmp([2;3]))/10*length(t)/60, length(t));
    % Berechnung der vollständigen Trajektorie
    t1 = tic();
    [Q_t_kk, QD_t_kk, QDD_t_kk, Phi_t_kk,~,~,~,Stats_kk] = RP.invkin2_traj( ...
      X_t_kk, XD_t_kk, XDD_t_kk, t, qs_kk, set_kk);
    fprintf(['Traj.-IK Fall %d (%s) berechnet. Dauer: %1.1fs für %d/%d ', ...
      'Bahnpunkte. %1.1fms pro i.O.-Punkt\n']', kk, name_method, toc(t1), ...
      Stats_kk.iter, length(t), 1e3*toc(t1)/Stats_kk.iter);
  else
    % Trajektorie direkt aus gespeicherter DP laden
    Q_t_kk = DP_TrajDetail_kk.Q;
    QD_t_kk = DP_TrajDetail_kk.QD;
    QDD_t_kk = DP_TrajDetail_kk.QDD;
    Phi_t_kk = DP_TrajDetail_kk.PHI;
    Stats_kk = DP_TrajDetail_kk.Stats;
    fprintf('Trajektorie aus gespeicherten Daten geladen\n');
  end
  if Stats_kk.errorcode ~= 0
    warning(['Fehler in Trajektorie zu groß. Zuerst bei Zeitschritt %d/%d ', ...
      '(t=%1.3fs). Traj.-IK nicht vollständig berechenbar'], Stats_kk.iter, ...
      size(Q_t_kk,1), t(Stats_kk.iter));
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
  X_ist(:,4:6) = denormalize_angle_traj(X_ist(:,4:6));
  % Plattform-Pose aus direkter Kinematik abspeichern
  XE_all(:,:,kk) = X_ist(:,1:6); % Nehme Plattform-Pose von erster Beinkette gesehen.
  XDE_all(:,:,kk) = XD_ist(:,1:6); % Nehme Plattform-Pose von erster Beinkette gesehen.
  % Gelenkkoordinaten normieren
  Q_t_norm_all(:,:,kk) = (Q_t_all(:,:,kk) - repmat(qlim(:,1)',n,1)) ... % untere Grenze abziehen
                          ./ repmat(qlim(:,2)'-qlim(:,1)',n,1); % geteilt durch Spannweite
  % Werte eintragen
  Hcond_all(:, 1, kk) = Stats_kk.condJ(:,1); % IK-Jacobi
  Hcond_all(:, 2, kk) = Stats_kk.condJ(:,2); % PKM-Jacobi
end

if ~usr_load_traj
  save(filename_traj, 'Q_t_all', 'QD_t_all', 'QDD_t_all', 'XE_all', ...
    'XDE_all', 'Q_t_norm_all', 'Hcond_all', 'Namen_Methoden', 'Namen_Methoden_Leg_Paper');
end
% For Debugging
save(fullfile(data_path, [filename_pre, '_all_data.mat']));

%% Paper-Bild: Konditionszahl-Karte über redundante Koordinate und Traj
% This creates Fig.11, a and Fig. 13, a
% Umrechnung von abort_thresh_h auf hnpos
abort_thresh_hpos = NaN(RP.idx_ik_length.hnpos, 1);
for f = fields(RP.idx_ikpos_hn)'
  if isfield(RP.idx_iktraj_hn, f{1})
    abort_thresh_hpos(RP.idx_ikpos_hn.(f{1})) = ...
      DP_settings.abort_thresh_h(RP.idx_iktraj_hn.(f{1}));
  end
end
% Redundanzkarte zeichnen
settings_perfmapplot = struct( ...
  'markermindist', [0.20, 15], ... % nur alle s=0.1 und phi=10° einen Marker
  'wn', DP_settings.wn, ...
  'abort_thresh_h', abort_thresh_hpos);
fighdl = change_current_figure(2400);clf;hold on;
set(fighdl, 'Name', 'PaperFigure_PerfDistr', 'NumberTitle', 'off');
[Hdl_all, s_pmp] = RP.perfmap_plot(H_all, phiz_range, s_ref, settings_perfmapplot);
xlabel('Normalized trajectory progress $s$', 'interpreter', 'latex');
ylabel('Redundant coordinate $\varphi_z$ in deg', 'interpreter', 'latex');

plot_lw = 0.75;
plot_ms = 6; % default
format = {[0 200 200]/255,  '', '-', 12, plot_lw, plot_ms; ... % cyan-dunkel
          [0 80 155 ]/255, '^', '-', 16, plot_lw, plot_ms; ... %imesblau
          ...[231 123  41]/255, 'o', '-', 5, plot_lw, plot_ms; ... %imesorange
          [0 200 0]/255, 'o', '-', 10, plot_lw, plot_ms; ... %grün
          ...[200 211  23]/255, 'o', '-', 10, plot_lw, plot_ms; ... %imesgrün (dunkler)
          'k', 's',  '--', 9, plot_lw, plot_ms; ...
          'b', 'x', '--', 7, plot_lw, plot_ms; ...
          'g', 'v', '-', 34, plot_lw, plot_ms};

linhdl = NaN(length(Namen_Methoden), 1);
for kk = 1:length(Namen_Methoden)
  X_ist = XE_all(:,1:6,kk);
  % insert trajectory into plot (with less dense sampling for saving space)
  I = select_plot_indices_downsample_nonuniform(PM_st, 180/pi*X_ist(:,6), 0.02, 1);
  linhdl(kk) = plot(PM_st(I), 180/pi*X_ist(I,6), 'g-'); % style will be overwritten
end
xlim([0, ceil(max(PM_s))]);
set(gca, 'xtick', 0:7);
ylim([-70, 120]);
linhdl_leg = line_format_publication(linhdl, format);
cbyh = ylabel(Hdl_all.cb,'Performance criterion $h$ (condition.)', ...
  'Rotation',90, 'interpreter', 'latex');
figure_format_publication(fighdl);
set(gca, 'Box', 'off');
set(fighdl, 'windowstyle', 'normal');
set_size_plot_subplot(fighdl, ...
  8,7,gca,... % 12.2 according to llncs.cls
  0.12,0.21,0.1,0.12,... %l r u d
  0,0) % x y
drawnow();
% Legende
I_vmactive = [2 4]; % Manuelle Auswahl der aktiven Marker. Referenz: s_pmp.violation_markers
LegLbl = Namen_Methoden_Leg_Paper;
LegHdl = linhdl_leg;
legendflex(LegHdl, LegLbl, 'anchor', {'n','n'}, ...
  'ref', fighdl, ... % an Figure ausrichten (mitten oben)
  'buffer', [0 -1], ... % Kein Versatz notwendig, da mittig oben
  'ncol', 0, 'nrow', 1, ... % eine Zeile für Legende
  'fontsize', 8, ...
  'xscale', 0.4, ... % Kleine Symbole
  'padding', [-1,-1,1], ... % Leerraum reduzieren
  'box', 'on');
% adjust axes labels text
ylh = get(gca, 'ylabel');
[x_off, x_slope] = get_relative_position_in_axes(gca, 'x');
[y_off, y_slope] = get_relative_position_in_axes(gca, 'y');
set(ylh, 'Position', [x_off+x_slope*(-1.2), y_off+y_slope*(-0.2), 0]);
% text with subfig number (a)
thdl = text(0,0,'(a)');
set(thdl, 'Position', [x_off+x_slope*(-1.35), y_off+y_slope*(-1.2), 0]);
set(thdl, 'FontWeight', 'bold');
figure_format_publication(fighdl); % nach letzter Erstellung von Text
% color bar text may look different in pdf than in figure (if not using
% latex interpreter above). Move label to the left to avoid cutting it off.
set(Hdl_all.cb, 'Position', [0.81, 0.1, 0.05, 0.8]); % put cb to the right
set(cbyh, 'Position', [2.7, 500, 0]); % put cblabel to the left
% Zweite Legende für Nebenbedingungen unten ins Bild. Erst hier, damit sie
% nicht in den Hintergrund gerät
if usr_trajnum == 1 % nur bei erster. Später als bekannt vorraussetzen
LegLbl2 = s_pmp.violation_markers(1,I_vmactive);
LegHdl2 = Hdl_all.VM(I_vmactive);
LegLbl2{strcmp(LegLbl2,'jac_cond')} = 'sing.';
LegLbl2{strcmp(LegLbl2,'coll_hyp')} = 'coll.';
legendflex(LegHdl2, LegLbl2, 'anchor', {'n','n'}, ...
  'ref', gca, ... % an Redundanzkarte ausrichten (mitten unten, im schwarzen Bereich)
  'buffer', [0 2], ... % Leichter Versatz, damit oberhalb von Subplot-Rahmen
  'anchor', [6 6], ... % mitte-mitte unten
  'ncol', 0, 'nrow', 1, ... % eine Zeile für Legende
  'fontsize', 8, ...
  'xscale', 0.7, ... % Kleine Symbole
  'padding', [1,1,1], ... % Leerraum reduzieren, wenn Beschriftung lang ist
  'box', 'on');
end
if usr_save_figures% && false % Debug: Do not save due to computational requirements
t1 = tic();

figname = sprintf('nullspace_traj%d',usr_trajnum);
saveas(fighdl, fullfile(data_path, [figname,'.fig']));
% For this to work, the Java heap memory has to be high enough. For
% high-resolution image 2GB not enough, 4GB worked.
% https://de.mathworks.com/matlabcentral/answers/255083-matlab-and-robotic-system-toolbox-java-lang-outofmemoryerror-gc-overhead-limit-exceeded#answer_318119
exportgraphics(fighdl, fullfile(paperfig_path, [figname,'.pdf']),'ContentType','vector') % ,'Resolution','100'
fprintf('Exported performance map as vector graphics. Duration: %1.1fs\n', toc(t1));
% exportgraphics(fighdl, fullfile(paperfig_path, ['nullspace_traj','.png']),'Resolution','600')  
% export_fig(fullfile(paperfig_path, ['nullspace_traj_exportfig','.pdf']), fighdl)
cd(paperfig_path);
ghostscript(['-dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 ', ...
  '-dPDFSETTINGS=/prepress -sOutputFile=',figname,'_compressed.pdf ',figname,'.pdf']);
end

%% Paper-Bild: Verlauf der Konditionszahl für verschiedene Trajektorien
% This creates Fig.11, b and Fig. 13, b
% Benutze nur die Konditionszahl und nicht die sonstigen Nebenbedingungen
Hcond_PKMJac_all = reshape(Hcond_all(:,2,:), n, length(Namen_Methoden));
% Integriere die Verläufe
Hcond_PKMJac_int_all = NaN(size(Hcond_PKMJac_all));
LegLblInt = cell(1,length(Namen_Methoden));
for i = 1:length(Namen_Methoden)
  Hcond_PKMJac_int_all(:,i) = sqrt(cumtrapz(PM_st, Hcond_PKMJac_all(:,i).^2)./PM_st);
  if isnan(Hcond_PKMJac_int_all(end, i))
    LegLblInt{i} = 'N/A';
  else
    LegLblInt{i} = sprintf('%1.0f', Hcond_PKMJac_int_all(end, i));
  end
end
% So formatieren, dass es neben die Redundanzkarte passt.
fighdl2 = change_current_figure(2500);clf;
axhdl2 = subplot(1,1,1); hold on;
hdl = NaN(length(Namen_Methoden), 1);
for i = 1:length(Namen_Methoden)
  % insert trajectory into plot (with less dense sampling for saving space)
  I = select_plot_indices_downsample_nonuniform(PM_st, Hcond_PKMJac_all(:,i), 0.02, 5);
  hdl(i)=plot(PM_st(I), Hcond_PKMJac_all(I,i));
end
LegHdl2=line_format_publication(hdl, format);
grid on;
xlabel('Norm. traj. progr. $s$', 'interpreter', 'latex');
% ylabel('Jacobian condition number $\mathrm{cond}(\boldmath{J})$', 'interpreter', 'latex');
set(gca, 'xtick', 0:7);
set(gca, 'YScale', 'log')
xlim([0, max(PM_st)+1e-3]);
if usr_trajnum == 1
  set(gca, 'ytick', [50, 75, 100, 150, 200:100:600]);
  ylim([50, 700]);
else
  set(gca, 'ytick', [100:100:500, 700, 900]);
end
% text with subfig number (b)
thdl = text(0,0,'(b)');
[x_off, x_slope] = get_relative_position_in_axes(gca, 'x');
[y_off, y_slope] = get_relative_position_in_axes(gca, 'y');
if usr_trajnum == 1 % manuell Position einstellen
  set(thdl, 'Position', [x_off+x_slope*(-1.33), y_off+y_slope*(-0.80), 0]);
else
  set(thdl, 'Position', [x_off+x_slope*(-1.33), y_off+y_slope*(-1.2), 0]);
end
set(thdl, 'FontWeight', 'bold');
figure_format_publication(fighdl2);
set_size_plot_subplot(fighdl2, ...
  4,7,axhdl2,...
  0.15,0.02, 0.01,0.12,... %l r u d
  0,0) % x y
I_valid = ~strcmp(LegLblInt, 'N/A');
if usr_trajnum == 1
  leghdl = legendflex(LegHdl2, LegLblInt, 'anchor', {'n','n'}, ...
    'ref', axhdl2, ... % an Axis ausrichten (mitten oben)
    'buffer', [0 -5], ... % Versatz so, dass zwischen den Linien
    'ncol', 1, 'nrow', 0, ... % eine Spalte für Legende
    'fontsize', 8, ...
    'xscale', 0.4, ... % Kleine Symbole
    'padding', [1,1,1], ... % Leerraum reduzieren
    'title', 'RMS', ...
    'box', 'on');
else
  leghdl = legendflex(LegHdl2(I_valid), LegLblInt(I_valid), 'anchor', {'n','n'}, ...
    'ref', fighdl2, ... % an Figure ausrichten (mitten oben)
    'buffer', [-4 0], ... 
    'ncol', 1, 'nrow', 0, ... % eine Spalte für Legende
    'anchor', [3 3], ...
    'fontsize', 8, ...
    'xscale', 0.4, ... % Kleine Symbole
    'padding', [1,1,1], ... % Leerraum reduzieren
    'title', 'RMS: ', ...
    'box', 'on');
end
% legch = get(leghdl, 'Children');
% set(legch(1), 'Interpreter', 'latex'); % erster Eintrag ist Titel
drawnow();
if usr_save_figures
  figname = sprintf('nullspace_traj%d_condition', usr_trajnum);
  saveas(2500, fullfile(data_path, [figname,'.fig']));
  exportgraphics(fighdl2, fullfile(paperfig_path, [figname,'.pdf']),'ContentType','vector')
%   exportgraphics(gcf, fullfile(paperfig_path, [figname,'.png']),'Resolution','600')
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
  fhld_kk = change_current_figure(100+30);clf;hold all;
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
