% Perform nullspace motion of a hexapod robot. This creates the results of
% Sec. 6.1 of the paper.
% Figures are saved into the directory paper/figures and update latex doc.
% 
% This file is adapted from the example ParRob_nullspace_convergence_test.m
% located at examples_tests/ParRob in this Git repository:
% https://github.com/SchapplM/robotics-toolbox

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2021-03
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clear
clc

%% User Settings and Initialization
if isempty(which('serroblib_path_init.m'))
  error('The serial robot database is not initialized in Matlab.');
end
if isempty(which('parroblib_path_init.m'))
  error('The parallel robot database is not initialized in Matlab.');
end
% User settings
%#ok<*UNRCH>
use_mex_functions = true; % use mex functions (much faster)
usr_recreate_mex = false; % recreate and recompile mex functions from templates
usr_plot_debug = false; % Create debug plots
usr_plot_objfun = true; % plot the objective function over the redundant coordinate
usr_plot_robot = true; % plot a sketch of the robot
usr_save_figures = true; % save figures to disk
usr_save_anim = false; % create an animation video of the robot motion
usr_figures_paper = true; % create and save the paper figures
usr_fast_debug = false; % skip some calculations for fast debugging
usr_opt_gains = false; % optimize the gains via particle swarm optimization. Has to be postprocessed manually
usr_select_pose = [1 2]; % 1=normal pose (Sec. 6.1.1), 2=singular pose (Sec. 6.1.2) as starting values
usr_test_timing = true; % multiple calls of ik functions to determine runtime for evaluation in Sec. 6.1.3.

% Further initialization:
respath = fileparts(which('ParRob_nullspace_static_pose.m'));
addpath(fullfile(respath, '..', 'gain_optimization'));
paperfig_path = fullfile(respath, '..', 'paper', 'figures');
% Namen of optimization criteria of the IK (für Diagramme)
hnames = {'qlim squared', 'qlim hyperb.', 'cond IK-Jac.', 'cond PKM-Jac.', 'sum'};
I_wn_traj = [1 2 5 6]; % Zuordnung zwischen Nebenbedingungen der Einzelpunkt- und Traj.-IK
% Threshold for deactivitating the hyperbolic potential to push the robot
% away from the limits in the IK nullspace optimization
optimcrit_limits_hyp_deact = 0.9;
%% Initializiation of the Robot Model
% Hexapod Robot
RP = parroblib_create_robot_class('P6RRPRRR14V3G1P4A1', 0.6, [0.200;0.100]);
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
%% Eckpunkte bestimmen, die angefahren werden
% Startpose in Arbeitsraum-Mitte
X0 = [ [0;0;0.6]; [0;0;0]*pi/180 ];
% Einstellungen für Einzelpunkt-IK.
s_ep = struct( ...
  'n_max', 5000, 'Phit_tol', 1e-12, 'Phir_tol', 1e-12); % , 'finish_in_limits', true
q0_ik1 = -0.5+rand(RP.NJ,1);
q0_ik1(RP.I_qa) = 1.0; % damit bei 6UPS das Schubgelenk positiv bleibt
q0_ik_fix = q0_ik1;
q0_ik_fix(RP.I1J_LEG(2):end) = NaN; % damit Übernahme Ergebnis Beinkette 1
% inverse Kinematik in Startpose. Gelenkgrenzen danach anpassen
[q0, Phi,~,Stats] = RP.invkin_ser(X0, q0_ik_fix, s_ep);
if any(abs(Phi)>1e-9), error('IK does not match'); end

% Bild: Roboter
if usr_plot_robot
  change_current_figure(100);clf;
  set(100,'Name','Rob','NumberTitle','off');
  title(sprintf('Roboter in Mittelstellung'));
  hold on;grid on;
  xlabel('x in m');ylabel('y in m');zlabel('z in m');
  view(3);
  s_plot = struct(  'ks_legs', [], 'straight', 0, 'mode', 4);
  RP.plot( q0, X0, s_plot );
end
% Select pose with strong tilting angle
XL = [ [0.05;0.03;0.6]; [30;-30;0]*pi/180 ]';

% TODO: Gelenkgrenzen nicht anpassen. Ist wesentlich aufwändiger in IK
% mit strengen Grenzen. 
for i = 1:RP.NLEG
  q0_i = q0(RP.I1J_LEG(i):RP.I2J_LEG(i));
  % TODO: Begrenze die Winkel der Kugel- und Kardangelenke auf +/- 30°
  % TODO: Geht erst, wenn Grenzen richtig eingehalten werden können.
  % Aktuell zu starke Überschwinger. Ignoriere Grenzen zunächst.
  RP.Leg(i).qlim(1:2,:) = q0_i(1:2) + repmat([-360, 360]*pi/180, 2, 1);
  RP.Leg(i).qlim(4:6,:) = q0_i(4:6) + repmat([-360, 360]*pi/180, 3, 1);
end
qlim   = cat(1, RP.Leg.qlim);
qDlim  = cat(1, RP.Leg.qDlim);
qDDlim = cat(1, RP.Leg.qDDlim);

%% Assemble Data for Table I of the Paper
fprintf('Robot Data (Table I in paper):\n');
fprintf('base diameter: %1.1fmm\n', 2*RP.DesPar.base_par(1)*1e3);
fprintf('platform diameter: %1.1fmm\n', 2*RP.DesPar.platform_par(1)*1e3);
fprintf('platform pair joint distance: %1.1fmm\n', RP.DesPar.platform_par(2)*1e3);
fprintf('platform position: [%s]mm\n', disp_array(1e3*XL(1,1:3),'%1.0f'));
fprintf('platform orientation: [%s]deg\n', disp_array(180/pi*XL(1,4:5),'%1.0f'));
fprintf('Stroke limits: %1.1f ... %1.1f mm\n', RP.Leg(1).qlim(3,1)*1e3, ...
  RP.Leg(1).qlim(3,2)*1e3);
fprintf('Prismatic (active joint) velocity limits: %1.1f m/s\n', RP.Leg(1).qDlim(3,2));
fprintf('Revolute (passive joint) velocity limits: %1.1f deg/s\n', 180/pi*RP.Leg(1).qDlim(6,2));
fprintf('Prismatic (active joint) acceleration limits: %1.1f m/s^2\n', RP.Leg(1).qDDlim(3,2));
fprintf('Revolute (passive joint) acceleration limits: %1.1f deg/s^2\n', 180/pi*RP.Leg(1).qDDlim(6,2));

%% Inverse Kinematik berechnen (Rasterung)
% Einstellungen für Dummy-Berechnung ohne Änderung der Gelenkwinkel.
s_ep_dummy = s_ep;
s_ep_dummy.retry_limit = 0;
s_ep_dummy.wn = ones(RP.idx_ik_length.wnpos,1); % hierdurch werden die Kriterien berechnet
s_ep_dummy.K = zeros(RP.NJ,1); % hierdurch keine Bewegung und damit ...
s_ep_dummy.Kn = zeros(RP.NJ,1); % ... sofortiger Abbruch
s_ep_dummy.optimcrit_limits_hyp_deact = optimcrit_limits_hyp_deact; % Hyperbolische Funktion nur nahe an Grenzen
ii_restab = 0; ii_restab_start = 0;
for k = 1:size(XL,1) % Schleife über alle Eckpunkte
  x_k = XL(k,:)';
  %% Get Distribution of the Performance Criteria by Global Discretization
  x_test_ges = []; %#ok<NASGU>
  h_ges = []; %#ok<NASGU>
  q_ges = []; %#ok<NASGU>
  if usr_plot_objfun
    nn = 360;
    if usr_figures_paper
      nn = nn * 4; % feinere Diskretisierung
    end
    if usr_fast_debug
      nn = 4;
    end
    q_ges = NaN(nn,RP.NJ);
    % set of platform poses for testing. Create increasing values
    x_test_ges = repmat(x_k',nn,1);
    x_test_ges(:,6) = [linspace(0,pi,nn/2),linspace(-pi,0,nn/2)];
    h_ges = NaN(nn,4);
    q_jj_old = q0;
    t0 = tic();
    for jj = 1:nn % loop over all discretized platform poses
      RP.update_EE_FG(I_EE_full, I_EE_full); % set robot to 3T3R
      x_jj = x_test_ges(jj,:)';
      t1=tic();
      q_jj_best = NaN(RP.NJ,1);
      if jj > 1
        % Get first estimate for step in joint space since the step in task
        % space is known and only an increment in one dimension.
        % Use normalization for step from pi to -pi.
        delta_x = [zeros(5,1);normalize_angle(x_test_ges(jj,6)-x_test_ges(jj-1,6))];
        [~,Phi_q] = RP.constr4grad_q(q_jj_old);
        [~,Phi_x] = RP.constr4grad_x(x_test_ges(jj-1,:)');
        Jtilde_inv_x = -Phi_q\Phi_x; % Full coordinate Jacobian (equ. 17 in paper)
        q0_jj = q_jj_old + Jtilde_inv_x*delta_x;
      else
        q0_jj = q0; % take middle pose as initial value
      end
      % Define IK settings
      s_ep_jj = s_ep; % leave as is (temporarily permit limit violations; no scaling until limit)
      s_ep_jj.normalize = false; % no normalization (due to joint limits)
      s_ep_jj.retry_on_limitviol = true; % assure joint limits by retrying, only they match
      [q_jj, Phi, ~, Stats] = RP.invkin_ser(x_jj, q0_jj, s_ep_jj);
      if jj == 1
        fprintf(['Start discretization of the nullspace with %d IK calls. ', ...
          'Duration for one call: %1.1fms. total estimated: %1.1fs\n'], ...
          nn, toc(t1)*1e3, (nn-1)*toc(t1));
      end
      if any(abs(Phi) > 1e-8)
        warning('IK did not find a solution for phi_z=%1.1fdeg', 180/pi*x_jj(6));
        continue % this IK configuration did not work.
      end
      RP.update_EE_FG(I_EE_full,I_EE_red); % Roboter auf 3T2R einstellen
      % IK benutzen, um Zielfunktionswerte zu bestimmen (ohne Neuberechnung)
      [q_dummy, Phi,~,Stats_dummy] = RP.invkin4(x_jj, q_jj, s_ep_dummy);
      if any(abs(q_jj - q_dummy) > 1e-8)
        error('IK-Ergebnis hat sich bei Test verändert');
      end
      h_ges(jj,:) = Stats_dummy.h(Stats_dummy.iter+1,2:5);
      q_ges(jj,:) = q_jj;
      q_jj_old = q_jj;
    end
    fprintf(['Calculation of %d discrete configurations in %1.1fs ', ...
      'done. %d/%d successful\n'], nn, toc(t0), sum(~isnan(h_ges(:,1))), nn);
    % Sortiere die Daten (andere Reihenfolge wegen Anfangswerten und
    % 2pi-Periodizität). Entferne dabei auch das Duplikat bei 0.
    [~,II] = unique(x_test_ges(:,6));
    x_test_ges = x_test_ges(II,:);
    h_ges = h_ges(II,:);
    q_ges = q_ges(II,:);
    
    % Verdoppele die Daten aufgrund der Periodizität
    x_test_ges = [x_test_ges(:,1:5),x_test_ges(:,6)-2*pi; x_test_ges;
                  x_test_ges(:,1:5),x_test_ges(:,6)+2*pi]; %#ok<AGROW>
    h_ges = repmat(h_ges,3,1);
    q_ges = repmat(q_ges,3,1);
    % Debug: Verteilung der Zielfunktionen über die redundante Koordinate
    % plotten. Wird jetzt weiter unten gemacht.
    change_current_figure(10);clf;
    optimcritaxhdl = NaN(3,1);
    for i = 1:4
      optimcritaxhdl(i) = subplot(2,2,i); hold on;
      h_ges_i_iO = h_ges(:,i);
      h_ges_i_iO(isinf(h_ges(:,2))) = NaN;
      plot(x_test_ges(:,6)*180/pi, h_ges_i_iO, 'g^-');
      h_ges_i_niO = h_ges(:,i);
      h_ges_i_niO(~isinf(h_ges(:,2))) = NaN;
      plot(x_test_ges(:,6)*180/pi, h_ges_i_niO, 'rv-');
      ylabel(sprintf('h %d (%s)', i, hnames{i})); grid on;
      xlabel('x6 in deg');
    end
    legend({'limits valid', 'limits violated'});
    sgtitle(sprintf('Objective Function Point %d', k));
    linkxaxes
    if usr_plot_debug
      change_current_figure(11);clf;
      for i = 1:RP.NJ
        subplot(6,6,i); hold on;
        plot(x_test_ges([1; end],6)*180/pi, qlim(i,1)*[1;1], 'r-');
        plot(x_test_ges([1; end],6)*180/pi, qlim(i,2)*[1;1], 'r-');
        plot(x_test_ges(:,6)*180/pi, q_ges(:,i));
        ylabel(sprintf('q %d', i)); grid on;
        xlabel('x6 in deg');
      end
      linkxaxes
      sgtitle('Joint Positions over Redundant Coordinate');
    end
  end

  %% Posen-Optimierung durchführen
  I_sing = islocalmax( h_ges(:,4) );
  I_onerev = x_test_ges(:,6)>-pi & x_test_ges(:,6) < pi;
  I_sing2 = I_sing & I_onerev & h_ges(:,4) > 1e3;
  fprintf(['Numeric detection of singularities from global distribution: ', ...
    'at phi_z=[%s] deg\n'], disp_array(x_test_ges(I_sing2,6)'*180/pi, '%1.2f'));
  % Manually selected values based on the given pose X0
  x6_l_range = [0, 33.8]*pi/180;
  x6_l_range = unique(x6_l_range, 'stable');
  num_ik_qs_successfull = 0;
  for l = usr_select_pose % Schleife über Anfangswerte der IK
    x_l = x_k;
    x_l(6) = x6_l_range(l);
    if l == 1
      % Erstes Beispiel (normale Bewegung). Benutze Antriebs-Koordinaten
      thresh_ns_qa = 1e3;
    else
      % Zweites Beispiel (aus Singularität). Benutze vollständige
      % Koordinaten
      thresh_ns_qa = 1;
    end
    fprintf('Start Inverse Kinematics for Pose %d (phiz0=%1.1fdeg)\n', l, 180/pi*x_l(6));
    % Inverse Kinematik zum Startpunkt
    RP.update_EE_FG(I_EE_full,I_EE_full); % Roboter auf 3T3R einstellen
    s_ep_start = s_ep;
    s_ep_start.scale_lim = 0.5; % Grenzen nicht verlassen (lieber neu versuchen)
    [qs, Phi_ep_s, ~, Stats_ep_s] = RP.invkin2(x_l, q0_ik_fix, s_ep_start);
    if any(abs(Phi_ep_s) > 1e-8)
      % Mit Einhaltung der Grenzen keine Lösung gefunden. Dann eben mit
      % Verletzung der Grenzen. Führt auch zu Abbruch, ist aber für
      % Statistik wichtig, ob es rechnerisch überhaupt funktioniert.
      s_ep_start = s_ep;
      [qs, Phi_ep_s, ~, Stats_ep_s] = RP.invkin2(x_l, q0_ik_fix, s_ep_start);
    end
    xs = RP.fkineEE_traj(qs')';
    ii_restab_start = ii_restab_start + 1;
    if any(abs(Phi_ep_s) > 1e-8)
      error('Inverse Kinematik konnte in Startpose für Punkt %d / Ori. %d nicht berechnet werden', k, l);
      % muss immer möglich sein für Paper (Erzwungene Einhaltung der Grenzen)
    end
    if any(qs(RP.I_qa&RP.MDH.sigma==1)<0)
      error('Die Schubgelenke sind umgeklappt. Darf nicht sein.');
    end
    qs_norm = (qs - qlim(:,1)) ./ (qlim(:,2)-qlim(:,1));
    if any(qs<qlim(:,1)) || any(qs>qlim(:,2))
      % Dieser Fall sollte nicht mehr auftreten, da meistens scale_lim benutzt wird
      % Die Grenzen müssen beim Start eingehalten werden (sonst
      % destabiliseren die Korrekturterme die IK)
      error('Startpose für Punkt %d / Ori. %d verletzt Grenzen', k, l);
      change_current_figure(999);clf;
      set(999,'Name','Rob_Debug','NumberTitle','off');
      title(sprintf('Startpose Punkt %d / Ori. %d', k, l));
      hold on; grid on; view(3);
      xlabel('x in m');ylabel('y in m');zlabel('z in m');
      s_plot = struct( 'straight', 0, 'mode', 4);
      RP.plot( qs, xs, s_plot );
    end
    num_ik_qs_successfull = num_ik_qs_successfull + 1;
    % IK benutzen, um Zielfunktionswerte zu bestimmen (ohne Neuberechnung)
    RP.update_EE_FG(I_EE_full,I_EE_red); % Roboter dafür auf 3T3R einstellen
    [q_dummy, Phi_dummy,~,Stats_dummy] = RP.invkin4(x_l, qs, s_ep_dummy);
    if any(abs(qs - q_dummy) > 1e-8)
      error('IK-Ergebnis hat sich bei Test verändert');
    end
    hs = Stats_dummy.h(Stats_dummy.iter,2:5)';
    % Füge Datenpunkt zu Gesamt-Rasterung hinzu
    x_test_ges = [x_test_ges; xs']; %#ok<AGROW>
    h_ges = [h_ges; hs']; %#ok<AGROW>
    [~,Isort] = sort(x_test_ges(:,6));
    x_test_ges = x_test_ges(Isort,:);
    h_ges = h_ges(Isort,:);

    %% IK für Einzelpunkt berechnen
    % Roboter auf 3T2R einstellen
    RP.update_EE_FG(I_EE_full,I_EE_red);
    s_ep_nolimits = s_ep;
    s_ep_nolimits.scale_lim = 0;
    wn_traj = zeros(RP.idx_ik_length.wntraj,1);
    s_ep_nolimits.wn = zeros(RP.idx_ik_length.wnpos,1); 
    s_ep_nolimits.wn(RP.idx_ikpos_wn.jac_cond) = 1;
    [q_ep_nolimits, Phi,~,Stats_ep_nolimits] = RP.invkin4(x_l, qs, s_ep_nolimits);
    x_ep_nolimits = RP.fkineEE_traj(q_ep_nolimits')';
    assert(all(abs(Phi)<1e-9), 'IK does not converge in test');
    s_ep_withlimits = s_ep_nolimits;
    s_ep_withlimits.scale_lim = 0.9;
    [q_ep_withlimits, Phi,~,Stats_ep_withlimits] = RP.invkin4(x_l, qs, s_ep_withlimits);
    x_ep_withlimits = RP.fkineEE_traj(q_ep_withlimits')';
    assert(all(abs(Phi)<1e-9), 'IK does not converge in test');
    fprintf(['Best condition number in position-level IK: h=%1.1f (phiz=%1.1f deg). ', ...
      'While considering joint limits h=%1.1f (phiz=%1.1f deg). Start at h=%1.1f\n'], ...
      Stats_ep_nolimits.h(Stats_ep_nolimits.iter+1,1+4), ...
      180/pi*x_ep_nolimits(6), ...
      Stats_ep_withlimits.h(Stats_ep_withlimits.iter+1,1+4), ...
      180/pi*x_ep_withlimits(6), Stats_ep_withlimits.h(1,1+4));

    %% Optimiere die Verstärkungsfaktoren

    if usr_opt_gains
      options = optimoptions('particleswarm');
      options.MaxIterations = 15;
      options.SwarmSize = 20;
      n = 2000; dt = 5e-3; % Länge der Trajektorie (eher kurz hier)
      Traj_t = (0:dt:(n-1)*dt)';
      Traj_X = repmat(x_l', n, 1);
      Traj_XD = zeros(n,6); Traj_XDD = Traj_XD;
      optstruct = struct('RP', RP, 'debug', false);
      optstruct.verbose = 1;
      optstruct.hopt_groundtruth = Stats_ep_withlimits.h(Stats_ep_withlimits.iter+1,1+4);
      optstruct.Traj_X = Traj_X;
      optstruct.Traj_XD = Traj_XD;
      optstruct.Traj_XDD = Traj_XDD;
      optstruct.Traj_t = Traj_t;
      optstruct.thresh_ns_qa = thresh_ns_qa;
      optstruct.qs = qs;
      fitnessfcn = @(p)fitness_gainopt(p, optstruct, options);
      clear fitness_gainopt
      t1 = tic();
      p_example = [0.1;0.5;0.1;1;1];
      fitnessfcn(p_example);
      t_fitness = toc(t1);
      fprintf('Time for call of the fitness function for PSO: %1.1fs\n', t_fitness);

      lb = zeros(length(p_example),1);
      ub = ones(length(p_example),1);
      fprintf('Estimated total duration of the PSO: %1.1f min\n', ...
        t_fitness*options.MaxIterations*options.SwarmSize/60);
      clear fitness_gainopt
      t1 = tic();
      [p_opt, f_opt, exitcode, output] = particleswarm(fitnessfcn, length(p_example), lb, ub, options);
      t_pso = toc(t1);
      fprintf('PSO finished. Duration: %1.1fmin\n', t_pso/60);
      [~, resultstats] = fitness_gainopt(NaN(3,1), optstruct, options);
      % Konvertiere das Format der Ergebnisse in eine Tabelle
      f_all = NaN(size(resultstats.p,1)*size(resultstats.p,3),1);
      p_all = NaN(size(resultstats.p,1)*size(resultstats.p,3),size(resultstats.p,2));
      kk = 0;
      for i = 1:size(resultstats.f,1) % gen (iterations)
        for j = 1:size(resultstats.f,2) % ind (swarm)
          kk = kk + 1;
          f_all(kk) = resultstats.f(i,j);
          p_all(kk,:) = resultstats.p(j,:,i); 
        end
      end
      I_nan = isnan(f_all);
      f_all = f_all(~I_nan);
      p_all = p_all(~I_nan,:);
      [~,I] = sort(f_all, 'ascend');
      ResTab_GainOpt = array2table([p_all(I,:), f_all(I)]);
      ResTab_GainOpt.Properties.VariableNames = {'K_v', 'K_p_cond', ...
        'K_d_cond', 'K_p_limits', 'K_d_limits', 'fval'};
      % Speichere die Ergebnisse ab
      resdir_gainopt = fullfile(respath, '..', 'gain_optimization', 'results');
      mkdirs(resdir_gainopt);
      save(fullfile(resdir_gainopt, ['gainopt_res', datestr(now,'yyyymmdd_HHMMSS'), '.mat']), ...
        'resultstats', 'x_l', 'qs', 'RP', 'ResTab_GainOpt');
      fprintf('Saved result of optimization to %s\n', resdir_gainopt);
      p_opt2 = p_all(I(1),:)';
      f_opt2 = f_all(I(1))';
      assert(all(abs(p_opt2-p_opt(:)) < 1e-8), ['Aus Diagnose geladene ', ...
        'beste Parameter stimmen nicht mit PSO-Ergebnis überein']);
      clear fitness_gainopt
      optstruct_debug = optstruct;
      optstruct_debug.debug = true;
      f_opt2_test = fitness_gainopt(p_opt2, optstruct_debug, options);
      assert(all(abs(f_opt2 - f_opt2_test) < 1e-8), ['Aus Diagnose geladener ', ...
        'Fitness-Wert stimmt nicht mit Neuberechnung überein']);
    end

    %% IK für Trajektorie berechnen
    Q_ii_norm_all = cell(3,1);
    QD_ii_norm_all = cell(3,1);
    QDD_ii_norm_all = cell(3,1);
    Traj_t_all = cell(3,1);
    Stats_traj_h_all = cell(3,1);
    X_ii_all = cell(3,1);
    if l == 1 %#ok<IFBDUP> % normale Pose
      n_cases_ii = 2;
    else % Singuläre Pose
      n_cases_ii = 2;
    end
    for ii = 1:n_cases_ii % Schleife über verschiedene Parametrierungen
      % Kriterien zusammenstellen
      wn_traj = zeros(RP.idx_ik_length.wntraj,1);
      % Dämpfung der Geschwindigkeit immer einbauen. Bei symmetrischen
      % Grenzen entspricht das dem Standard-Dämpfungsterm aus der Literatur
      wn_traj(RP.idx_iktraj_wnP.qDlim_par) = 0.5; % K_v
      % Zusätzlich Dämpfung bei Überschreitung des Grenzbereichs zu den
      % Positions-Grenzen. Dadurch weniger Überschreitungen der Grenze.
      wn_traj(RP.idx_iktraj_wnP.qlim_hyp) = 1; % K_P (hyperb. limit)
      if l == 1 % normale Pose
        switch ii
          case 1
            wn_traj(RP.idx_iktraj_wnP.qDlim_par) = 0; % K_v
            wn_traj(RP.idx_iktraj_wnP.jac_cond) = 1; % K_P (cond)
            wn_traj(RP.idx_iktraj_wnD.jac_cond) = 0; % K_D (cond)
          case 2
            wn_traj(RP.idx_iktraj_wnP.qDlim_par) = 0.8; % K_v
            wn_traj(RP.idx_iktraj_wnP.jac_cond) = 1; % K_P (cond)
            wn_traj(RP.idx_iktraj_wnD.qlim_hyp) = 0.5; % K_D (limit)
            wn_traj(RP.idx_iktraj_wnD.jac_cond) = 0.5; % K_D (cond)
        end
      else % singuläre Pose
        switch ii
          case 1 % manuell getuned
            wn_traj(RP.idx_iktraj_wnP.qDlim_par) = 0.03; % K_v
            wn_traj(RP.idx_iktraj_wnP.jac_cond) = 0.05; % K_P (cond)
            wn_traj(RP.idx_iktraj_wnD.jac_cond) = 0.01; % K_D (cond)
          case 2 % aus PSO
            wn_traj(2) = 0; % K_P (limit)
            wn_traj(RP.idx_iktraj_wnP.qDlim_par) = 0.2; % K_v
            wn_traj(RP.idx_iktraj_wnP.jac_cond) = 0.5; % K_P (cond)
            wn_traj(RP.idx_iktraj_wnD.qlim_hyp) = 0.6; % K_D (limit)
            wn_traj(RP.idx_iktraj_wnD.jac_cond) = 0.03; % K_D (cond)
          otherwise
            warning('Fall nicht definiert');
            continue
        end
      end
      % Roboter auf 3T2R einstellen
      RP.update_EE_FG(I_EE_full,I_EE_red);

      % IK mit Einzelpunkt-Verfahren berechnen
      s_ep_ii = s_ep;
      if usr_save_anim % sehr feinschrittige Bewegungen (für flüssige Animation)
        s_ep_ii.maxrelstep = 0.005;
        s_ep_ii.n_max = s_ep_ii.n_max * 2;
      end
      s_ep_ii.retry_limit = 0; %nur ein Versuch. Sonst zufällig neue Gelenkwinkel.
      % Grenzen dürfen auch in Zwischenschritten nicht überschritten
      % werden.
      s_ep_ii.scale_lim = 0.9;
      s_ep_ii.wn = ones(RP.idx_ik_length.wnpos,1);
      s_ep_ii.wn(1:4) = wn_traj(I_wn_traj); % TODO: Indizes nicht hart-kodiert aufstellen
      s_ep_ii.wn(s_ep_ii.wn~=0) = 1; % Immer auf 1 setzen. Nur eine Zielfunktion betrachtet. Ist für Summe egal.
      s_ep_ii.optimcrit_limits_hyp_deact = optimcrit_limits_hyp_deact; % Hyperbolische Funktion nur nahe an Grenzen
      t_ep = NaN(51,1);
      if ii == 1 % nur ein mal als Vergleich
        for iit = 1:(1+50*usr_test_timing)
        t1 = tic();
        [q_ep_ii, Phi,~,Stats_ep] = RP.invkin4(x_l, qs, s_ep_ii);
        t_ep(iit) = toc(t1);
        end
        if usr_test_timing
          t_ps = t_ep / Stats_ep.iter;
          fprintf(['Timing Evaluation of position-level kinematics: in total: ', ...
            '\n\tavg %1.2fms, std %1.6fms, n=%d, data [%s]ms (ignore first); per sample: ', ...
            '\n\tavg %1.2fms, std %1.6fms, n=%d, data [%s]ms (ignore first); %d samples\n'], ...
            1e3*mean(t_ep(2:end)), 1e3*std(t_ep(2:end)), length(t_ep(2:end)), disp_array(1e3*t_ep(:)', '%1.2f'), ...
            1e3*mean(t_ps(2:end)), 1e3*std(t_ps(2:end)), length(t_ps(2:end)), disp_array(1e3*t_ps(:)', '%1.2f'), Stats_ep.iter);
          % re-generate function to always compute all joints of the robot
          fname = which('P6RRPRRR14V3_invkin3.m');
          fid  = fopen(fname,'r');
          f=fread(fid,'*char')';
          fclose(fid);
          % change line which determines which gradient calculation to use
          f = strrep(f,'if all(abs(Phi)<1e-3) && taskred_rotsym','if false % for testing the timing using always all joint variations');
          fid  = fopen(fname,'w');
          fprintf(fid,'%s',f);
          fclose(fid);
          % compile modified ik function
          matlabfcn2mex({'P6RRPRRR14V3_invkin3'});
          % Test timing again
          t_ep2 = NaN(50,1);
          for iit = 1:(1+49*usr_test_timing)
          t1 = tic();
          [q_ep_ii2, Phi2,~,Stats_ep2] = RP.invkin4(x_l, qs, s_ep_ii);
          t_ep2(iit) = toc(t1);
          % Do not test if the result stays the same. The calculation is
          % different, since the platform pose is not regarded in this case
%           assert(all(abs(q_ep_ii2-q_ep_ii)<1e-8), 'result using the different nullspaces has to be equal');
          end
          t_ps2 = t_ep2 / Stats_ep2.iter;
          fprintf(['Timing Evaluation of position-level kinematics without ', ...
            'using efficient form of gradient approximation: in total: ', ...
            '\n\tavg %1.2fms, std %1.6fms, n=%d, data [%s]ms; per sample: ', ...
            '\n\tavg %1.2fms, std %1.6fms, n=%d, data [%s]ms; %d samples\n'], ...
            1e3*mean(t_ep2), 1e3*std(t_ep2), length(t_ep2), disp_array(1e3*t_ep2(:)', '%1.2f'), ...
            1e3*mean(t_ps2), 1e3*std(t_ps2), length(t_ps2), disp_array(1e3*t_ps2(:)', '%1.2f'), Stats_ep2.iter);
          % Change function back to the original version
          parroblib_create_template_functions({RP.mdlname(1:end-2)}, false, false);
          matlabfcn2mex({[RP.mdlname(1:end-6), '_invkin3']});
        end
      end
      assert(all(abs(Stats_ep.PHI(1,RP.I_constr_red))<1e-9), ...
        'Residuum im ersten Schritt muss Null sein');
      assert(all(all(abs(Stats_ep.PHI(1:1+Stats_ep.iter,RP.I_constr_red))<1e-3)), ...
        'Residuum während der Annäherung an Ziel muss klein bleiben. Aktuell max. %1.1e', ...
        max(max(Stats_ep.PHI(1:1+Stats_ep.iter,RP.I_constr_red))));
      % IK mit Trajektorien-Verfahren berechnen. Setze virtuelle
      % Trajektorie, die keine Soll-Vorgaben hat. Dadurch entsteht eine
      % reine Nullraumbewegung.
      if usr_fast_debug
        n = 2000; dt = 5e-3; % Länge der Trajektorie
      else
        n = 10000; dt = 1e-3; % Länge der Trajektorie
      end
      Traj_t = (0:dt:(n-1)*dt)';
      Traj_X = repmat(x_l', n, 1);
      Traj_XD = zeros(n,6); Traj_XDD = Traj_XD;
      assert(length(Traj_t)==n, 'Zeit-Basis der virtuellen Trajektorie ist falsch');
      s_traj_ii = struct('wn', wn_traj);
      s_traj_ii.thresh_ns_qa = thresh_ns_qa;
      t_traj = NaN(51,1);
      for iit = 1:(1+50*usr_test_timing)
      t1 = tic();
      [Q_ii, QD_ii, QDD_ii, Phi_ii,~,~,~,Stats_traj] = RP.invkin2_traj(Traj_X, Traj_XD, Traj_XDD, Traj_t, qs, s_traj_ii);
      t_traj(iit) = toc(t1);
      end
      % Kürze die Trajektorie, falls Bewegung vorzeitig abgeklungen ist
      % (oder abgebrochen wurde)
      I_firstnan = find(any(isnan(QD_ii),2),1,'first');
      if ~isempty(I_firstnan)
        warning('Abbruch der Trajektorienberechnung bei I=%d', I_firstnan);
        I_finishacc = I_firstnan-1;
      else
        I_noacc = all(abs(QDD_ii)<1e-8,2);
        I_finishacc = find(I_noacc==0,1,'last');
      end
      if usr_test_timing
        I_lastcalc = min([I_firstnan, size(QD_ii,1)]);
        t_ps = t_traj / I_lastcalc;
        fprintf(['Timing Evaluation of trajectory kinematics:', ...
          '\n\tin total: avg %1.2fs, std %1.6fs, n=%d, data [%s]s (ignore first);', ...
          '\n\tper sample: avg %1.2fms, std %1.6fms, data [%s]ms (ignore first); %d samples\n'], ...
          mean(t_traj(2:end)), std(t_traj(2:end)), length(t_traj(2:end)), disp_array(t_traj(:)', '%1.2f'), ...
          1e3*mean(t_ps(2:end)), 1e3*std(t_ps(2:end)), disp_array(1e3*t_ps(:)', '%1.2f'), I_lastcalc);
        % Run the code without mex and with the profiler to assess where
        % the computation time is spent in the code.
        RP.fill_fcn_handles(false);
        profile on
        RP.invkin2_traj(Traj_X, Traj_XD, Traj_XDD, Traj_t, qs, s_traj_ii);
        profsave(profile('info'),fullfile(respath,sprintf('profile_results_case%d_set%d', l, ii)));
        if use_mex_functions, RP.fill_fcn_handles(true); end % set back to compiled functions
      end
      fprintf(['Pt. %d/ Ori. %d/ Case %d. IK computed: %d steps position-', ...
        'level-IK (%1.1fs; %1.1fms); %d steps Traj.-IK (%1.1fs; %1.1fms)\n'], k, l, ...
        ii, Stats_ep.iter, t_ep(1), 1e3*t_ep(1)/Stats_ep.iter, I_finishacc, t_traj(1), 1e3*t_traj(1)/I_finishacc);
      Q_ii = Q_ii(1:I_finishacc,:);
      QD_ii = QD_ii(1:I_finishacc,:);
      QDD_ii = QDD_ii(1:I_finishacc,:);
      Phi_ii = Phi_ii(1:I_finishacc,:);
      Traj_t = Traj_t(1:I_finishacc);
      Stats_traj.h = Stats_traj.h(1:I_finishacc,:);
      q_traj_ii = Q_ii(end,:)'; % Ergebnis der Trajektorien-IK
      % Traj-Ergebnis kürzen (wenn Ende mit Toleranz erreicht ist)
      I_finalvalue = all(abs(repmat(q_traj_ii',I_finishacc,1)-Q_ii) < 1e-10,2) & ...
                     abs(repmat(Stats_traj.h(end,1),I_finishacc,1)-Stats_traj.h(:,1)) < 1e-10;
      I_finish = find(I_finalvalue==0,1,'last');
      Traj_t = Traj_t(1:I_finish,:);
      Q_ii = Q_ii(1:I_finish,:);
      QD_ii = QD_ii(1:I_finish,:);
      QDD_ii = QDD_ii(1:I_finish,:);
      Phi_ii = Phi_ii(1:I_finish,:);
      Stats_traj.h = Stats_traj.h(1:I_finish,:);

      % Speichere EE-Pose resultierend aus Gelenkwinkeln aus dem
      % Konvergenzverlauf. Benutze bei Einzelpunkt-IK nur i.O.-Posen
      % (keine Posen, bei denen die Ketten nicht geschlossen sind).
      X_ii = RP.fkineEE2_traj(Q_ii);
      Stats_ep.X = RP.fkineEE2_traj(Stats_ep.Q);
      I_invalid = any(abs(Stats_ep.PHI(:,RP.I_constr_red))>1e-3,2);
      Stats_ep.X(I_invalid,:) = NaN;

      % TODO: Um 2pi verschoben, so dass am nächsten an Startwert
      % (aktuell springt das Ergebnis im Bild, falls es über +/- pi geht.
      x_traj_ii = RP.fkineEE_traj(q_traj_ii')';
      x_ep_ii = RP.fkineEE_traj(q_ep_ii')';

      % Summe der positionsbezogenen Leistungsmerkmale der Traj. am Ende.
      % Schließe geschwindigkeitsbezogene Merkmale aus Vergleich aus.
      % Dadurch ist der Vergleich von Positions- und Traj.-IK möglich.
      h_traj_ii = Stats_traj.h(I_finish,:);
      h_traj_ii_sum = sum(Stats_traj.h(I_finish,1+I_wn_traj)*s_ep_ii.wn(1:4)); % TODO: Keine hart kodierten Indizes
      % Leistungsmerkmale der Positions-IK am Ende
      h_ep_ii = Stats_ep.h(Stats_ep.iter+1,:);
      % Bilde Summe der Merkmale neu (eigentlich nur notwendig, wenn Ge-
      % wichtungen zwischen Trajektorie und Einzelpunkt unterschiedlich
      % sind aber hier trotzdem die Ergebnisse verglichen werden sollen.
      h_ep_ii_sum = sum(Stats_ep.h(Stats_ep.iter,2:1+length(s_ep_ii.wn))'.*s_ep_ii.wn);%wn_traj(I_wn_traj));
      % Ergebnis prüfen
      reserr_q = q_traj_ii - q_ep_ii;
      reserr_h_sum = h_traj_ii_sum - h_ep_ii_sum;
      I_wnact = s_ep_ii.wn ~= 0; % Indizes der aktiven Nebenbedingungen
      hs_ii = Stats_traj.h(1,1); % Startwert der Nebenbedingungen (bei Traj.)
      step_h_traj = h_traj_ii_sum - hs_ii; % Verbesserung der Nebenbedingungen bei Traj.
      step_h_ep = h_ep_ii_sum - hs_ii; % Verbesserung der Nebenbedingungen bei Einzelpunkt
      % Bestimme relativen Fehler (bezogen auf die Methode mit dem
      % besseren Schritt; Schritte sind negativ, da Minimierung)
      reserr_h_rel = reserr_h_sum/min(step_h_traj,step_h_ep);
      fprintf(['Point %d/ Orientation %d (x6=%1.1f deg). Result from Traj.-IK is ', ...
        '%1.1e worse than position-level IK (%1.5f vs %1.5f). Improvement ', ...
        'vs start: %1.3f bzw. %1.3f\n'], ...
        k, l, 180/pi*x_l(6), h_traj_ii_sum - h_ep_ii(1), ...
        h_traj_ii_sum, Stats_ep.h(Stats_ep.iter,1), -step_h_traj, -step_h_ep);

      % Prüfe hier, ob ein Fehler vorliegt und breche dann ab, nachdem
      % die Bilder gezeichnet wurden.
      raise_error_h = abs(reserr_h_rel)>5e-2;
      if isnan(raise_error_h)
        warning('Ungültiges Ergebnis');
        ResStat.Error(ii_restab) = 5;
        raise_error_h = true;
      end
      if any(abs(reserr_h_sum) > 1e-3) && ...
          ~(l==1 && ii==1) % do not show warning in case of voluntary oscillations
        warning('Zielfunktion weicht absolut bei beiden Methoden ab. Aber kein Fehler, da eventuell auch Verzweigung der Lösung.');
      end
      if raise_error_h
        warning('Zielfunktion weicht bei beiden Methoden zu stark voneinander ab (relativer Fehler %1.1f%%)', reserr_h_rel*100);
      end
      if step_h_traj >= 0
        % In Traj.-IK Schwingungen, wenn kein PD-Regler benutzt wird.
        warning('Nebenbedingungen haben sich bei Traj.-IK verschlechtert');
        raise_error_h = true;
      end
      if step_h_ep >= 0
        warning('Nebenbedingungen haben sich bei Einzelpunkt.-IK verschlechtert');
        raise_error_h = true;
      end
      % Normiere die Größen
      Q_ii_norm = (Q_ii - repmat(qlim(:,1)',size(Q_ii,1),1)) ./ ...
                  repmat(qlim(:,2)'-qlim(:,1)',size(Q_ii,1),1);
      QD_ii_norm = (QD_ii - repmat(qDlim(:,1)',size(QD_ii,1),1)) ./ ...
                  repmat(qDlim(:,2)'-qDlim(:,1)',size(QD_ii,1),1);
      QDD_ii_norm = (QDD_ii - repmat(qDDlim(:,1)',size(QDD_ii,1),1)) ./ ...
                  repmat(qDDlim(:,2)'-qDDlim(:,1)',size(QDD_ii,1),1);
      Stats_ep.Q_norm = (Stats_ep.Q - repmat(qlim(:,1)',size(Stats_ep.Q,1),1)) ./ ...
                  repmat(qlim(:,2)'-qlim(:,1)',size(Stats_ep.Q,1),1);
      % Prüfe, ob Grenzen eingehalten werden. Ist teilweise nicht anders
      % möglich, wenn Grenzen in mehreren Gelenken gleichzeitig verletzt
      % werden.
      if any(QDD_ii_norm(:) < -1) || any(QDD_ii_norm(:)>2)
        warning(['Gelenkbeschleunigungsgrenzen werden in Traj. massiv ', ...
          'verletzt (qDD norm: %1.1f bis %1.1f)'], min(QDD_ii_norm(:)), ...
          max(QDD_ii_norm(:)));
      end
      if any(QD_ii_norm(:) < -0.01) || any(QD_ii_norm(:)>1.01)
        warning(['Gelenkgeschwindigkeitsgrenzen werden in Traj. ', ...
          'verletzt (qD norm: %1.4f bis %1.4f)'], min(QD_ii_norm(:)), ...
          max(QD_ii_norm(:)));
      end
      if any(Q_ii_norm(:) < 0) || any(Q_ii_norm(:)>1)
        warning(['Gelenkpositionsgrenzen werden in Traj. ', ...
          'verletzt (q norm: %1.4f bis %1.4f)'], min(Q_ii_norm(:)), ...
          max(Q_ii_norm(:)));
      end
      % Ergebnisse speichern
      Traj_t_all{ii} = Traj_t;
      Q_ii_norm_all{ii} = Q_ii_norm;
      QD_ii_norm_all{ii} = QD_ii_norm;
      QDD_ii_norm_all{ii} = QDD_ii_norm;
      Stats_traj_h_all{ii} = Stats_traj.h;
      X_ii_all{ii} = X_ii;
    end % ii (Berechnung der Nullraumbewegung)
    save(fullfile(respath, sprintf('pkm_nullspace_results_case%d.mat', l)));
    %% Create Debug plots
    % Create a virtual time base for the position-level IK
    % Simple case: Equidistant and normalized to trajectory end time. This
    % works well for high velocities in the trajectory.
%     Stats_ep.T = (0:1:Stats_ep.iter)'/(Stats_ep.iter)*Traj_t(end);
    % use maximal joint velocity to create the time base
    Q_ep = Stats_ep.Q(1:1+Stats_ep.iter,:);
    Qdiff_ep = diff(Q_ep);
    Tdiff_ep_allj = abs(Qdiff_ep)./repmat(abs(qDlim(:,2))', size(Qdiff_ep,1), 1);
    Tdiff_ep = max(Tdiff_ep_allj, [], 2);
    % sample to 1ms steps
    Tdiff_ep_sampled = max([Tdiff_ep, repmat(1e-3, size(Tdiff_ep,1), 1)], [], 2);
    Tdiff_ep_sampled = ceil(Tdiff_ep_sampled/1e-3)*1e-3;
    Stats_ep.T = [0;cumsum(Tdiff_ep_sampled)];
    Stats_ep.QD = [zeros(1,RP.NJ); diff(Stats_ep.Q(1:Stats_ep.iter+1,:))...
                  ./repmat(diff(Stats_ep.T), 1, RP.NJ); ...
                  zeros(size(Stats_ep.Q,1)-Stats_ep.iter-1,RP.NJ)];
    % append time steps until the end (for markers in later plots)
    T_add = Stats_ep.T(end)+1e-3:1e-3:(Stats_ep.T(end)+ size(Stats_ep.Q,1)*1e-3);
    n_rest = size(Stats_ep.Q,1)-Stats_ep.iter-1;
    Stats_ep.Q(1+Stats_ep.iter+1:end,:) = repmat(Stats_ep.Q(1+Stats_ep.iter,:), n_rest, 1);
    Stats_ep.h(1+Stats_ep.iter+1:end,:) = repmat(Stats_ep.h(1+Stats_ep.iter,:), n_rest, 1);
    Stats_ep.T(1+Stats_ep.iter+1:size(Stats_ep.Q,1),:) = T_add(1:n_rest);
    legend_entries = {'ii=1', 'ii=2', 'ii=3'};
    for ii = 1:n_cases_ii % Plotten (Debug)
      filename_pre = sprintf('pkm_nullspace_motion_ori%d_case%d_', l, ii);
      Traj_t = Traj_t_all{ii};
      if isempty(Traj_t), continue; end % keine Daten vorliegend.
      Q_ii_norm = Q_ii_norm_all{ii};
      QD_ii_norm = QD_ii_norm_all{ii};
      QDD_ii_norm = QDD_ii_norm_all{ii};
      Stats_traj_h = Stats_traj_h_all{ii};
      X_ii = X_ii_all{ii};
      % Vergleich
      if usr_plot_debug
        % Bild: Gelenkpositionen
        change_current_figure(1);set(1,'Name','q','NumberTitle','off');
        if ii == 1, clf; end
        for i = 1:RP.NJ
          subplot(6,6,i); hold on;
          if ii == 1
            q_i_int = Stats_ep.Q(1,i)+cumtrapz(Stats_ep.T, Stats_ep.QD(:,i));
            q_i_int_norm = (q_i_int-qlim(i,1))/(qlim(i,2)-qlim(i,1));
            plot(Stats_ep.T(1:Stats_ep.iter+1), Stats_ep.Q_norm(1:Stats_ep.iter+1,i));
            plot(Stats_ep.T, q_i_int_norm);
          end
          plot(Traj_t, Q_ii_norm(:,i));
          ylabel(sprintf('q %d (norm)', i)); grid on;
          xlabel('Time in s');
        end
        linkxaxes
        legend(['Pos.-Level','Pos.-Level (int.)',legend_entries(1:ii)]);
        sgtitle('Joint Positions');
        if usr_save_figures && ii == 2
          saveas(1, fullfile(respath, [filename_pre,'jointpositions.fig']));
          saveas(1, fullfile(respath, [filename_pre,'jointpositions.png']));
        end

        % Bild: Gelenk-Geschwindigkeiten
        change_current_figure(11);set(11,'Name','qD','NumberTitle','off');
        if ii == 1, clf; end
        for i = 1:RP.NJ
          subplot(6,6,i); hold on;
          if ii == 1
            qD_ep_i = Stats_ep.QD(:,i);
            qD_ep_i_norm = (qD_ep_i-qDlim(i,1))/(qDlim(i,2)-qDlim(i,1));
            plot(Stats_ep.T, qD_ep_i_norm);
          end
          plot(Traj_t, QD_ii_norm(:,i));
          ylabel(sprintf('qD %d (norm)', i)); grid on;
          xlabel('Time in s');
        end
        linkxaxes
        legend(['Pos.-Level',legend_entries(1:ii)]);
        sgtitle('Joint Velocities (Traj.-IK)');
        if usr_save_figures && ii == 2
          saveas(11, fullfile(respath, [filename_pre,'jointvelocities.fig']));
          saveas(11, fullfile(respath, [filename_pre,'jointvelocities.png']));
        end
        
        % Bild: Gelenk-Beschleunigungen
        change_current_figure(12);set(12,'Name','qDD','NumberTitle','off');
        if ii == 1, clf; end
        for i = 1:RP.NJ
          subplot(6,6,i); hold on;
          plot(Traj_t, QDD_ii_norm(:,i));
          ylabel(sprintf('qDD %d (norm)', i)); grid on;
          xlabel('Time in s');
        end
        linkxaxes
        legend(legend_entries(1:ii));
        sgtitle('Joint Accelerations (Traj.-IK)');
        if usr_save_figures && ii == 2
          saveas(12, fullfile(respath, [filename_pre,'jointaccelerations.fig']));
          saveas(12, fullfile(respath, [filename_pre,'jointaccelerations.png']));
        end

        % Bild: Zielkriterien (Zeitverlauf)
        change_current_figure(2);set(2,'Name','h','NumberTitle','off');
        if ii == 1, clf; end
        for i = 1:4
          subplot(2,2,i); hold on;
          if ii == 1
            plot(Stats_ep.T(1:Stats_ep.iter), Stats_ep.h(1:Stats_ep.iter,1+i));
          end
          plot(Traj_t, Stats_traj_h(:,1+I_wn_traj(i)));
          ylabel(sprintf('h %d (%s) (wn=%1.1f)', i, hnames{i}, s_ep_ii.wn(i))); grid on;
          xlabel('Time in s');
        end
        linkxaxes
        legend(['Pos.-Level',legend_entries(1:ii)]);
        sgtitle('Objectives of IK');
        if usr_save_figures && ii == 2
          saveas(2, fullfile(respath, [filename_pre,'objectives.fig']));
          saveas(2, fullfile(respath, [filename_pre,'objectives.png']));
        end

        % Bild: Redundante Koordinate
        change_current_figure(25);
        if ii == 1, clf; end
        if ii == 1
          plot(Stats_ep.T(1:Stats_ep.iter+1), 180/pi*Stats_ep.X(1:Stats_ep.iter+1,6));
        end
         hold on;
        plot(Traj_t, 180/pi*X_ii(:,6));
        legend(['Pos.-Level',legend_entries(1:ii)]); grid on; % TODO: passt nicht
        xlabel('Time in s');
        ylabel('Redundant coordinate x6 in deg');
        set(25,'Name','x6','NumberTitle','off');
        if usr_save_figures && ii == 2
          saveas(25, fullfile(respath, [filename_pre,'redcoordX.fig']));
          saveas(25, fullfile(respath, [filename_pre,'redcoordX.png']));
        end
      end

      % In Zielfunktions-Bild eintragen
      if usr_plot_objfun
        change_current_figure(20);set(20,'Name','h_x6','NumberTitle','off');
        if ii == 1
          clf;
          h_ges_jj_sum = zeros(size(h_ges,1),1);
          for jj = 1:5
            subplot(2,3,jj); hold on; grid on;
            if jj < 5
              h_ges_jj = h_ges(:,jj);
              h_ges_jj_sum = h_ges_jj_sum + h_ges_jj*s_ep_ii.wn(jj);
            else
              h_ges_jj = h_ges_jj_sum;
            end
            plot(x_test_ges(:,6)*180/pi, h_ges_jj);
            ylabel(sprintf('h %d (%s)', jj, hnames{jj})); grid on;
            xlabel('x6 in deg');
          end
        end
        for jj = find(s_ep_ii.wn(1:4)~=0)'
          subplot(2,3,jj); hold on;
          switch ii
            case 1
              marker = 'r+-';
            case 2
              marker = 'bo--';
            case 3
              marker = 'cd:';
          end
          % Debug: Ortskurve der Zwischenzustände einzeichnen
          if ii == 1
            plot(Stats_ep.X(:,6)*180/pi,Stats_ep.h(:,1+jj), 'gx-');
          end
          plot(X_ii(:,6)*180/pi,Stats_traj_h(:,1+I_wn_traj(jj)), marker);
        end
        legend(['Global', 'Pos.-Level',legend_entries(1:ii)]); grid on;
        linkxaxes
        if usr_save_figures
          saveas(20, fullfile(respath, [filename_pre,'Ortskurve_Ziel_RedKoord.fig']));
          saveas(20, fullfile(respath, [filename_pre,'Ortskurve_Ziel_RedKoord.png']));
        end
      end
    end
    %% Create Animation Video
    if usr_save_anim
      for ii = 0:n_cases_ii % create animation for all methods/parameterizations
        maxduration_animation = 5; % Dauer der Animation als mp4 (in s)
        if ii == 0 % position-level IK
          t_Vid = (0:1/30*(Stats_ep.iter/maxduration_animation):Stats_ep.iter)';
          I_anim = knnsearch( (1:Stats_ep.iter)' , t_Vid ); % Berechne Indizes in Traj.-Zeitstempeln
          Q_anim = Stats_ep.Q(I_anim,:);
          X_anim = Stats_ep.X(I_anim,:);
          anim_name = 'Animation_Position';
        else % trajectory IK with different settings
          Traj_t = Traj_t_all{ii};
          if isempty(Traj_t), continue; end % keine Daten vorliegend.
          Q_ii_norm = Q_ii_norm_all{ii};
          Q_ii = repmat(qlim(:,1)', size(Q_ii_norm,1),1) + Q_ii_norm.*...
            (repmat(qlim(:,2)'-qlim(:,1)', size(Q_ii_norm,1),1));
          Stats_traj_h = Stats_traj_h_all{ii};
          
          X_ii = X_ii_all{ii};
          t_Vid = (0:1/30*(Traj_t(end)/maxduration_animation):Traj_t(end))';
          I_anim = knnsearch( Traj_t , t_Vid );
          Q_anim = Q_ii(I_anim,:);
          X_anim = X_ii(I_anim,:);
          anim_name = sprintf('Animation_Traj_Case%d', ii);
        end
        
        change_current_figure(200);clf;hold all;
        if strcmp(get(200, 'windowstyle'), 'docked'), set(200, 'windowstyle', 'normal'); end
        set(200, 'name', 'Anim', ...
          'color','w', 'NumberTitle', 'off', 'units','normalized',...
          'outerposition',[0 0 1 1]); % Vollbild
        view(3); axis auto; hold on; grid on;
        xlabel('x in m');ylabel('y in m');zlabel('z in m');
        if length(I_anim) > max(I_anim)
          warning(['Es gibt nicht genug Zwischenschritte für ein ', ...
            'flüssiges Video (%d Video-Bilder, %d IK-Schritte)'], length(I_anim), max(I_anim));
        end
        s_plot = struct( 'ks_legs', [], 'straight', 0, 'mode', 4);
        s_anim = struct('mp4_name', [fullfile(respath, [filename_pre, anim_name]),'.mp4'] );
        RP.anim( Q_anim, X_anim, s_anim, s_plot);
      end
    end

    %% Create Robot Plot for Paper: Fig. 1 (Sub-Figures)
    if usr_figures_paper
      for ii = 1:2
        change_current_figure(100);clf;
        set(100,'Name','Rob','NumberTitle','off');
        hold on;grid on;
        xlabel('x in m');ylabel('y in m');zlabel('z in m');
        set(gca,'XTICKLABEL',{});set(gca,'YTICKLABEL', {});set(gca,'ZTICKLABEL',{});
        set(gca,'xtick',[],'ytick',[],'ztick',[]);
        set(get(gca, 'XAxis'), 'visible', 'off');
        set(get(gca, 'YAxis'), 'visible', 'off');
        set(get(gca, 'ZAxis'), 'visible', 'off');
        view([-66.8, 6.8]); % manuelle Einstellung
        s_plot = struct(  'ks_legs', [], 'straight', 0, 'mode', 4);
        if ii == 1
          x = Stats_ep.X(1,:)';
          q = Stats_ep.Q(1,:)';
        else
          x = Stats_ep.X(Stats_ep.iter+1,:)';
          q = Stats_ep.Q(Stats_ep.iter+1,:)';
        end
        % Set EE transformation to avoid occlusion of the platform frame
        % by the platform visualisation
        r_P_E_backup = RP.r_P_E;
        RP.update_EE([0;0;10e-3]);
        RP.plot( q, x, s_plot );
        RP.update_EE(r_P_E_backup);
        ch = get(gca, 'Children');
        for jj = 1:length(ch)
          % KS-Texte entfernen
          if strcmp(ch(jj).Type, 'hgtransform')
            chjjch = ch(jj).Children;
            % Entferne den KS-Text wieder
            delete(chjjch(1:4));
          end
        end
        figure_format_publication(gca);
        set(gca, 'Box', 'off');
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
            % remove the cuboid representation of prismatic joints
            if ~isempty(regexp(get(ch(jj), 'DisplayName'), '^Leg_([\d]+)_Joint_3$', 'match'))
              delete(ch(jj));
            end
          end
        end
        
        set(100, 'windowstyle', 'normal');
        set_size_plot_subplot(100, ...
          6,6,gca,...
          0,0,0,0,0,0)
        drawnow();
        % Bild speichern
        if usr_save_figures
          if ii == 1
            filename_pre = sprintf('pkm_case%d_start_pose', l);
          else
            filename_pre = sprintf('pkm_case%d_end_pose', l);
          end
          exportgraphics(gcf, fullfile(paperfig_path, [filename_pre,'.png']),'Resolution',600)
          saveas(100, fullfile(respath, [filename_pre,'.fig']));
        end
      end
    end
    %% Create Plots for Paper: Fig. 2 and Fig. 3
    colors = [[0 80 155 ]/255; ...%imesblau
              [231 123 41 ]/255;... %imesorange
              [0 0 0];... % schwarz
              [200 211 23 ]/255];... %imesgrün
    format = cell(4,5);
    markerlist = {'^', 'o', 'x', 's'};
    linewidth = 1.5;
    markernumber = [12, 25, 18, 45];
    for fff = 1:1+n_cases_ii
      format{fff,1} = colors(fff,:);
      format{fff,2} = markerlist{fff};
      format{fff,3} = '-';
      format{fff,4} = markernumber(fff);
      format{fff,5} = linewidth;
    end

    linhdl1 = NaN(1+n_cases_ii,1);
    for ii = 1:n_cases_ii % Plotten (Paper)
      Traj_t = Traj_t_all{ii};
      Q_ii_norm = Q_ii_norm_all{ii};
      QD_ii_norm = QD_ii_norm_all{ii};
      QDD_ii_norm = QDD_ii_norm_all{ii};
      Stats_traj_h = Stats_traj_h_all{ii};
      X_ii = X_ii_all{ii};
      %% Kombiniertes Bild: 3 Teile: h/t, phiz/t, h/phiz
      if usr_figures_paper
        fighdl = change_current_figure(1000); set(1000, 'Name', sprintf('case%d_overview', l), 'numbertitle', 'off');
        if strcmp(get(fighdl, 'windowstyle'), 'docked'), set(fighdl, 'windowstyle', 'normal'); end
        if ii == 1, clf; end
        axhdl = NaN(1,3);
        axhdl(1)=subplot(1,3,1); hold on; % h/t
        if ii == 1
          linhdl1(1) = plot(Stats_ep.T, Stats_ep.h(:,1+4), ...
            'Color', colors(1,:), 'linewidth', linewidth);
        end
        linhdl1(1+ii) = plot(Traj_t, Stats_traj_h(:,1+I_wn_traj(4)), ...
          'Color', colors(1+ii,:), 'linewidth', linewidth);
        ylabel(sprintf('$h=\\mathrm{cond}$(\\boldmath${J}_{x}$)'), ...
          'interpreter', 'latex');
        grid on;
        xlabel('time in s');
        if l == 1
          set(axhdl(1), 'xlim', [0, 3]);
          set(axhdl(1), 'ylim', [53, 110]);
        else
          set(axhdl(1), 'xlim', [0, 3]);
          set(axhdl(1), 'ylim', [47, 210]);
        end
        if ii == n_cases_ii
          line_format_publication(linhdl1, format);
        end

        if ii == 1, linhdl=NaN(n_cases_ii+1,1); end
        axhdl(2)=subplot(1,3,2); hold on; % phiz/t
        if ii == 1
          linhdl(1) = plot(Stats_ep.T, 180/pi*Stats_ep.X(:,6), ...
            'Color', colors(1,:), 'linewidth', linewidth);
        end
        linhdl(1+ii) = plot(Traj_t, 180/pi*X_ii(:,6), ...
          'Color', colors(1+ii,:), 'linewidth', linewidth);
        ylabel(sprintf('$\\varphi_z$ in deg'), 'interpreter', 'latex');
        grid on;
        xlabel('time in s');
        set(axhdl(2), 'xlim', get(axhdl(1), 'xlim'));
        if l == 1
          set(axhdl(2), 'ylim', [-30, 5]);
        else
          set(axhdl(2), 'ylim', [-30, 35]);
        end
        if ii == n_cases_ii
          linleghdl = line_format_publication(linhdl, format);
        end
        axhdl(3)=subplot(1,3,3); hold on; % h/phiz
        if ii == 1
          linhdl3 = NaN(2,1);
          % Startwert
          linhdl3(1) = plot(x_l(6)*180/pi, Stats_traj_h(1,1+I_wn_traj(4)), 'rs');
          % Gesamter Verlauf
%             I_iO = ~isinf(h_ges(:,2));
          linhdl3(2) = plot(x_test_ges(:,6)*180/pi, h_ges(:,4), 'r--', 'linewidth', linewidth);
%             linhdl3(3) = plot(x_test_ges(:,6)*180/pi, h_ges(:,4), 'r--', 'linewidth', linewidth);
        end
        if ii == 1
          plot(180/pi*Stats_ep.X(1:Stats_ep.iter+1,6), ...
            Stats_ep.h(1:Stats_ep.iter+1,1+4), 'Color', colors(1,:), 'linewidth', linewidth);
          plot(180/pi*Stats_ep.X(Stats_ep.iter+1,6), ...
            Stats_ep.h(Stats_ep.iter+1,1+4), markerlist{1}, 'Color', colors(1,:), 'linewidth', linewidth);
        end
        plot(180/pi*X_ii(end,6), Stats_traj_h(end,1+I_wn_traj(4)), ...
          markerlist{1+ii}, 'Color', colors(1+ii,:), 'linewidth', linewidth);
        % https://de.mathworks.com/matlabcentral/answers/100493-how-do-i-obtain-bold-formatted-mathematical-expressions-in-matlab-using-latex
        ylabel(sprintf('$h=\\mathrm{cond}$(\\boldmath${J}_{x}$)'), ...
          'interpreter', 'latex');
        xlabel(sprintf('$\\varphi_z$ in deg'), 'interpreter', 'latex');
        grid on;
        set(axhdl(3), 'xlim', get(axhdl(2), 'ylim'));
        set(axhdl(3), 'YScale', 'log')
        if l == 1
          set(axhdl(3), 'xlim', [-200, 200]);
          set(axhdl(3), 'xtick', [-180, -90, 0, 90, 180]);
          xtickangle(0) % werden Standardmäßig schräg gestellt. Wieder gerade machen.
          set(axhdl(3), 'ylim', [min(h_ges(:,4))-1, 2e4]);%max(h_ges(:,4))]);
        else
          set(axhdl(2), 'ytick', [-30, 0, 30]);
          set(axhdl(3), 'xtick', [-30, 0, 30]);
          set(axhdl(3), 'ylim', [min(h_ges(:,4))-1, max(h_ges(:,4))*1.3]);%max(h_ges(:,4))]);
        end
        set(axhdl(3), 'ytick', [1e2, 1e3, 1e4, 1e5]);
        if ii == n_cases_ii
          % text with subfig number (a/b/c)
          thdl = NaN(3,1);
          for jth = 1:3
            axes(axhdl(jth)); %#ok<LAXES>
            thdl(jth) = text(0,0,sprintf('(%s)', char(96+jth)));
          end
          for jth = 1:3
            [~, ~, x_pos] = get_relative_position_in_axes(axhdl(jth), 'x', -1.4);
            [~, ~, y_pos] = get_relative_position_in_axes(axhdl(jth), 'y', -1.25);
            set(thdl(jth), 'Position', [x_pos, y_pos, 0]);
            set(thdl(jth), 'FontWeight', 'bold');
          end
          % create zoom box into the existing axis object. Has to be after
          % resizing the plot, because this changes the z order
          if l == 1
            zoomhdl1 = add_zoom_axis(axhdl(1), [1.5, 2.8; 55.5, 58.5]);
            set(zoomhdl1, 'position', [0.15, 0.44, 0.17, 0.45]);
            arrhdl1 = annotation(fighdl,'arrow');
            set(arrhdl1, 'position', [0.25, 0.43, -0.00, -0.17]);
            zoomhdl2 = add_zoom_axis(axhdl(2), [1.5, 2.8; -28, -20]);
            set(zoomhdl2, 'position', [0.47, 0.49, 0.17, 0.40]);
            arrhdl2 = annotation(fighdl,'arrow');
            set(arrhdl2, 'position', [0.58, 0.48, -0.00, -0.1]);
          end
          if l == 2
            zoomhdl1 = add_zoom_axis(axhdl(1), [1, 2.5; 51, 70]);
            set(zoomhdl1, 'position', [0.18, 0.45, 0.14, 0.45]);
            arrhdl1 = annotation(fighdl,'arrow');
            set(arrhdl1, 'position', [0.24, 0.43, -0.0, -0.14]);
            zoomhdl2 = add_zoom_axis(axhdl(2), [1.6, 2.2; -28, -19]);
            set(zoomhdl2, 'position', [0.50, 0.50, 0.15, 0.40]);
            arrhdl2 = annotation(fighdl,'arrow');
            set(arrhdl2, 'position', [0.58, 0.49, -0.00, -0.19]);
          end
          drawnow();
          figure_format_publication(axhdl);
          set_size_plot_subplot(fighdl, ...
            15.8, 5.5, ...
            axhdl, ...
            0.07, 0.01, 0.09, 0.15, ... % l r u d
            0.07, 0.011) % x y
          drawnow();
          if l == 1
            lh = legend(linleghdl, {'position-level IK', 'trajectory IK with $K_\mathrm{D}=K_\mathrm{v}=0$', ...
              'trajectory IK with tuned gains'}, 'interpreter', 'latex');
          else
            lh = legend(linleghdl, {'position-level IK', 'trajectory IK with gain set 1', ...
              'trajectory IK with gain set 2'}, 'interpreter', 'latex');
          end
          set(lh, 'orientation', 'horizontal', 'position', ...
            [0.07,0.94,0.92,0.05]); % x y b h
          % Do not plot an additional legend for the performance criterion
          if l == 1
            lh3 = legend(linhdl3, {'initial value', 'computed offline'});
            set(lh3, 'position', [0.7831    0.7134    0.1922    0.1383], 'orientation', 'vertical'); % x y b h
          end
        end
        if ii == n_cases_ii && usr_save_figures
          filename_pre = sprintf('pkm_nullspace_case%d',l);
          saveas(fighdl, fullfile(respath, [filename_pre,'_overview.fig']));
          exportgraphics(fighdl, fullfile(paperfig_path, ...
            [filename_pre,'_overview.pdf']),'ContentType','vector');
          fprintf('Figure saved: %s\n', filename_pre);
        end
      end
    end % for ii (Plotten Testfälle der IK)
    fprintf('Case %d successfully generated\n', l);
  end % for l (Drehung des Endeffektors)
end % for k (Punkte)
