% Example Figure for State-Inteval Dynamic Programming (Sect. 4.2 and 4.3 of LNEE paper)
% This creates Fig. 5-8 and 9-10 of the paper
% 
% Preliminaries:
% Run case_study/ParRob_dynprog.m to create necessary files.
% 
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-02
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear
% User settings
usr_highres_distrfig = true; % high resolution of the paper figure for performance criterion map
usr_overlapmode = false; % Switch between Sect. 4.2 and 4.3
%% Initialize
% Initialize Robot
% Default output directory (for paper)
paperfig_path = fileparts(which('fig_dynprog_interv_example.m'));
this_dir = fileparts(which('fig_dynprog_interv_example.m'));
assert(~isempty(this_dir), 'The script currently run has to be in the PATH');
data_path = fullfile(this_dir, '..', '..', 'case_study', 'data_LNEE');
d = load(fullfile(data_path, 'robot_definition.mat'));
RP = d.RP;
parroblib_addtopath({RP.mdlname});
% Load performance map
if usr_highres_distrfig, resstr = '_highres';
else,                    resstr = '_lowres'; end
filename_pre = sprintf('LNEE_pkm_traj1_%dsamples', 12263); % has to match trajectory from other file
filename_perfmap = fullfile(data_path, [filename_pre, '_perfmap', resstr, '.mat']);
assert(exist(filename_perfmap, 'file'), sprintf(['performance map file ', ...
  'does not exist in %s'], data_path));
d = load(filename_perfmap);
H_all = d.H_all;
XL = d.XL;
phiz_range = d.phiz_range;
s_ref = d.s_ref;
s_tref = d.s_tref;
if ~usr_overlapmode % Standard-Modus ohne überlappende Intervalle
  filename_post = ['_dynprog_costRMStraj_n9_red', '', '.mat'];
else % Mit überlappenden Intervallen
  filename_post = ['_dynprog_costRMStraj_n7_red_overlap', '', '.mat'];
end
filename_dynprog= fullfile(data_path, [filename_pre, filename_post]);
assert(exist(filename_dynprog, 'file'), 'dynamic programming results file does not exist');
d = load(filename_dynprog);
DP_XE = d.DP_XE;
DP_Stats = d.DP_Stats;
DP_TrajDetail = d.DP_TrajDetail;
DP_settings = d.DP_settings;
overlapstr = '';
if ~usr_overlapmode 
  debugfoldername_full = 'LNEE_Traj1_DP_debug_costRMStraj_n9_red';
else
  debugfoldername_full = 'LNEE_Traj1_DP_debug_costRMStraj_n7_red_overlap';
  overlapstr= '_overlap';
end
dpres_dir = fullfile(this_dir, '..', '..', 'case_study', debugfoldername_full);
% Use saved data in project folder
% dpres_dir = fullfile(fileparts(which('robsynth_projektablage_path.m')), ...
%   '06_Publikationen/2022_LNEE_3T2R_DynProg/DP_Ergebnisse/LNEE_DP_debug_costRMStraj_n9_red');
% Use network share directly
% dpres_dir = ['/home/moritz/SeaDrive/Für mich freigegeben/imes-projekt-dfg_robotersynthese/', ...
%     '06_Publikationen/2022_LNEE_3T2R_DynProg/DP_Ergebnisse/LNEE_DP_debug_costRMStraj_n9_red'];
assert(exist(dpres_dir, 'file'), 'directory with debug information for DP does not exist');

% Eigenschaften der Intervalle aus Ergebnis ausgeben lassen. Siehe dynprog_taskred_ik
delta_phi = DP_Stats.delta_phi;
phi_range = DP_Stats.phi_range;
fprintf('%d reference interaval centers actually used with step size %1.1f°: [%s]°\n', ...
  length(phi_range), 180/pi*delta_phi, disp_array(phi_range*180/pi, '%1.1f'));
%% Prepare performance map plot
% Umrechnung auf hnpos
abort_thresh_hpos = NaN(RP.idx_ik_length.hnpos, 1);
for f = fields(RP.idx_ikpos_hn)'
  if isfield(RP.idx_iktraj_hn, f{1})
    abort_thresh_hpos(RP.idx_ikpos_hn.(f{1})) = ...
      DP_settings.abort_thresh_h(RP.idx_iktraj_hn.(f{1}));
  end
end
settings_perfmap = struct( ...
  'markermindist', [0.10, 10], ... % nur alle s=... und phi=...° einen Marker
  'phiz_range_plot', pi/180*[-80, 150], ... % nur das Plotten, was später auch ins Paper kommt. Sonst lange Rechenzeit
  'extend_map', false, ... % sowieso nicht sichtbar
  'wn', DP_settings.wn, ...
  'abort_thresh_h', abort_thresh_hpos);
d_init = load(fullfile(dpres_dir, 'dp_init.mat'));
%% Plot performance map with first three states of the first stage (Fig. 5+6 and 9)
fprintf('Zeichne Bilder für SI-DP Einzeltransfer\n');
for i_stage1 = [1 2] % create one figure for each of the first two stages
  if usr_overlapmode && i_stage1 == 2, break; end % zweites Bild nicht zeichnen
  pmfhdl = change_current_figure(i_stage1); clf;
  set(pmfhdl, 'Name', sprintf('SIDP_singletransfer_%d%s', i_stage1, overlapstr), 'NumberTitle', 'off');
  axhdl = NaN(1,3);
  DP_hdl = NaN(4,1); % Handle für die verschiedenen Linien (für Legende)
  i_ax = 0;
  % Load data for stage
  d_state = load(fullfile(dpres_dir, sprintf('dp_stage%d_final.mat', i_stage1)));
  if i_stage1 == 1
    if ~usr_overlapmode
      i_state2_range = [1 4 7];
    else % Wähle Zustände so aus, dass man die Überlappung sieht
      i_state2_range = [3 4 9]; % 6 liegt zwischen 2 und 3 (da von 5 an neu beginnend)
    end
    i_state1 = 1;
  else
    i_state2_range = [3 4 5];
    % Select manually from existing transfers in debug folder
    i_state1 = 4;
  end
  for i_state2 = i_state2_range % select illustrative state transfers
    i_ax = i_ax + 1;
    axhdl(i_ax) = subplot(1,3,i_ax); hold on;
    s_stage = DP_settings.PM_s_tref(DP_settings.IE(i_stage1):DP_settings.IE(i_stage1+1));
    t_stage = d_init.t(DP_settings.IE(i_stage1):DP_settings.IE(i_stage1+1));
    settings_perfmap_stage = settings_perfmap;
    settings_perfmap_stage.s_range_plot = minmax2(s_stage') + [-0.1, 0.1];
    t1 = tic();
    [Hdl_all, s_pmp] = RP.perfmap_plot(H_all, phiz_range, s_ref, settings_perfmap_stage);
    fprintf('Redundanzkarte gezeichnet: %1.1fs\n', toc(t1));
    VM_hdl = Hdl_all.VM;
    % load data for transfer
    d_ii = load(fullfile(dpres_dir, sprintf('dp_stage%d_state%d_to%d_result.mat', ...
      i_stage1, i_state1, i_state2)));
    % Plot Reference trajectory
    Iplot = select_plot_indices_downsample_nonuniform(...
      s_stage, d_ii.X6_traj_ref, 0.02, 1*pi/180);
    DP_hdl(1) = plot(s_stage(Iplot), 180/pi*d_ii.X6_traj_ref(Iplot), 'k--'); % TODO: _ref
    % Plot tolerance band for optimization variable
    DP_hdl(2) = plot(s_stage(Iplot), 180/pi*(d_ii.X6_traj_ref(Iplot)+...
      interp1(d_ii.xlim6_interp(1,:), d_ii.xlim6_interp(2,:), t_stage(Iplot), 'spline')), 'k-');
    plot(s_stage(Iplot), 180/pi*(d_ii.X6_traj_ref(Iplot)+...
      interp1(d_ii.xlim6_interp(1,:), d_ii.xlim6_interp(3,:), t_stage(Iplot), 'spline')), 'k-');
    % Plot actual trajectory
    Iplot = select_plot_indices_downsample_nonuniform(...
      s_stage(1:length(d_ii.X6_traj)), d_ii.X6_traj, 0.02, 1*pi/180);
    hdl = plot(s_stage(Iplot), 180/pi*d_ii.X6_traj(Iplot), 'k-');
    set(hdl, 'LineWidth', 2);
    if ~isinf(d_state.F_stage(i_state1, i_state2))
      set(hdl, 'Color', 'c');
    end
    % Dummy-Plots for legend
    DP_hdl(3) = plot(NaN,NaN,'k-', 'LineWidth', 2);
    DP_hdl(4) = plot(NaN,NaN,'c-', 'LineWidth', 2);
    ylabel('red. coord. $\varphi_z$ in deg', 'interpreter', 'latex');
    % Put text for state interval
    plot(i_stage1, 180/pi*phi_range, 'k_', 'MarkerSize', 8);
    plot(i_stage1, 180/pi*(phi_range+delta_phi/2), 'k_', 'MarkerSize', 5);
    for kk = 1:length(phi_range)
      if 180/pi*phi_range(kk) < -50 || 180/pi*phi_range(kk) > 100
        continue % nur Zustände einzeichnen, die im Plot-Bereich liegen
      end
      if usr_overlapmode == 0 % Nur Zahl rechts an Achse schreiben
        text(i_stage1+0.06, 180/pi*phi_range(kk)-2, sprintf('%d', kk), ...
          'interpreter', 'latex');
      else % Vollständige Bezeichnung (da genug Platz da ist)
        text(i_stage1+0.09, 180/pi*phi_range(kk)-30, sprintf('$[x_{\\mathrm{ref},%d}]$', kk), ...
          'Rotation', 90, 'interpreter', 'latex');
      end

    end
    quiver(i_stage1*ones(1,length(phi_range)), 180/pi*phi_range, ...
      zeros(1,length(phi_range)),  180/pi*0.9*delta_phi/2*ones(1,length(phi_range)), 'off', 'k.', 'LineWidth', 1);
    quiver(i_stage1*ones(1,length(phi_range)), 180/pi*phi_range, ...
      zeros(1,length(phi_range)), -180/pi*0.9*delta_phi/2*ones(1,length(phi_range)), 'off', 'k.', 'LineWidth', 1);

    xlim(minmax2(s_stage'));
    ylim([-70, 120]);
    if i_ax ~= 2
      xlabel('');
    else
      xlabel('Normalized trajectory progress $s$', 'interpreter', 'latex');
    end
    if i_ax > 1
      ylabel('');
    end
    if i_ax < 3
      delete(colorbar());
    end
    if i_state2 <= length(phi_range) % Normale Intervalle
      targetinterval = sprintf('[x_{\\mathrm{ref},%d}]', i_state2);
    else % Überlappendes Intervall, andere Bezeichnung
      targetinterval = sprintf('[x_{\\mathrm{add},%d}]', i_state2-length(phi_range));
    end
    title(sprintf('\\textbf{(%s)} $f(x_%d,u_%d) {\\in} %s$', ...
      char(96+i_ax), i_stage1-1, i_stage1, targetinterval), 'interpreter', 'latex');
  end
  % Format plot and save pdf for paper. Distances different in PDF than in
  % Matlab figure window!
  set(Hdl_all.cb, 'Position', [0.91, 0.1, 0.02, 0.7]); % put cb to the right
  cyh=ylabel(Hdl_all.cb, 'Perf. crit.  $h$ (cond.)', 'Rotation', 90, 'interpreter', 'latex');
  set(cyh, 'Position', [3.4, 500, 0]); % put cblabel to the left
  figure_format_publication(pmfhdl);
  remove_inner_labels(axhdl, 2);
  set_size_plot_subplot(pmfhdl, ...
    12.2,4,axhdl,... % 12.2 according to llncs.cls
    0.08,0.14,0.23,0.19,... %l r u d
    0.04,0) % x y
  drawnow();
  % Legende
  I_vmactive = [2 4]; % Manuelle Auswahl der aktiven Marker. Referenz: s_pmp.violation_markers
  LegHdl = [DP_hdl; VM_hdl(I_vmactive)];
  LegLbl = ['reference', 'border', 'invalid', 'valid', s_pmp.violation_markers(1,I_vmactive)];
  % LegLbl{strcmp(LegLbl,'qlim_hyp')} = 'joint limit';
  LegLbl{strcmp(LegLbl,'jac_cond')} = 'singularity';
  LegLbl{strcmp(LegLbl,'coll_hyp')} = 'collision';
  legendflex(LegHdl, LegLbl, 'anchor', {'n','n'}, ...
    'ref', pmfhdl, ... % an Figure ausrichten (mitten oben)
    'buffer', [0 -1], ... % Kein Versatz notwendig, da mittig oben
    'ncol', 0, 'nrow', 1, ... % eine Zeile für Legende
    'fontsize', 8, ...
    'xscale', 0.6, ... % Kleine Symbole
    ... 'padding', [-3,-3,3], ... % Leerraum reduzieren
    'box', 'on');
  pdfname = sprintf('dp_interv%s_stage%d_3cases', overlapstr, i_stage1);
  exportgraphics(pmfhdl, fullfile(paperfig_path, [pdfname, '.pdf']),'ContentType','vector');
  fprintf('Bild gespeichert: %s\n', pdfname);
  cd(paperfig_path);
  ghostscript(['-dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 ', ...
    '-dPDFSETTINGS=/prepress -sOutputFile=',pdfname,'_compressed.pdf ',pdfname,'.pdf']);
end

%% Plot performance map for the first two stages (Fig. 7 and 10)
fprintf('Zeichne Redundanzkarte für erste 2 Stufen\n');
pmfhdl = change_current_figure(1); clf;
set(pmfhdl, 'Name', sprintf('SIDP_Stage1and2%s',overlapstr), 'NumberTitle', 'off');
axhdl = NaN(1,2);
DP_hdl = NaN(3,1); % Handle für die verschiedenen Linien (für Legende)
d_final = load(fullfile(dpres_dir, 'dp_final.mat'), 'I_all', 'I_best');
for i_stage1 = 1:length(axhdl)
  axhdl(i_stage1) = subplot(1,2,i_stage1); hold on;
  s_stage = DP_settings.PM_s_tref(DP_settings.IE(i_stage1):DP_settings.IE(i_stage1+1));
  settings_perfmap_stage = settings_perfmap;
  settings_perfmap_stage.s_range_plot = minmax2(s_stage') + [-0.1, 0.1];
  t1 = tic();
  [Hdl_all, s_pmp] = RP.perfmap_plot(H_all, phiz_range, s_ref, settings_perfmap_stage);
  fprintf('Redundanzkarte gezeichnet: %1.1fs\n', toc(t1));
  VM_hdl = Hdl_all.VM;
  % Plot lines for first stage
  d_state = load(fullfile(dpres_dir, sprintf('dp_stage%d_final.mat', i_stage1)));
  tfn = dir(fullfile(dpres_dir, sprintf('dp_stage%d_state*_to*_result.mat', i_stage1)));
  % Gehe alle Transfers durch, die gespeichert sind
  hdl_all = NaN(length(tfn), 2); % Erste Spalte Handle, zweite Spalte Kategorie als Zahl
  for ii = 1:length(tfn)
    d_ii = load(fullfile(dpres_dir, tfn(ii).name));
    i_state1 = d_ii.k;
    i_state2 = d_ii.l;
    % Abtastung der Datenpunkte für Plot spärlicher für geringe Dateigröße
    Iplot = select_plot_indices_downsample_nonuniform(...
      s_stage(1:length(d_ii.X6_traj)), d_ii.X6_traj, 0.03, 2*pi/180);
    hdl_all(ii,1) = plot(s_stage(Iplot), 180/pi*d_ii.X6_traj(Iplot), 'k-');
    set(hdl_all(ii,1), 'LineWidth', 1);
    if usr_overlapmode && i_state2 > length(phi_range)
      set(hdl_all(ii,1), 'LineStyle', '--'); % Teil der überlappenden Linien
    end
    % Prüfe ob es sich um die optimale Teil-Politik handelt
    if d_final.I_all(i_stage1+1,i_state2) == i_state1
      hdl_all(ii,2) = 1;
      set(hdl_all(ii,1), 'Color', 'm');
      if i_state2 <= length(phi_range)
        DP_hdl(3) = hdl_all(ii,1); % optimal transmission (to this stage)
      end
    elseif ~isinf(d_state.F_stage(i_state1, i_state2))
      hdl_all(ii,2) = 2;
      set(hdl_all(ii,1), 'Color', 'c');
      if i_state2 <= length(phi_range)
        DP_hdl(2) = hdl_all(ii,1); % line for valid transition
      end
    else
      hdl_all(ii,2) = 3;
      DP_hdl(1) = hdl_all(ii,1); % line for invalid transition
    end
  end
  % Gehe alle Handles durch. Zuerst alle cyan nach hinten, dann alle
  % schwarzen, dann Farbkarte
  for ii = find(hdl_all(:,2)==2)' % 2=cyan
    ZOrderSet(hdl_all(ii,1), 0); % Setze Linie ganz nach unten (damit optimale immer oben sind)
  end
  for ii = find(hdl_all(:,2)==3)' % 3=black
    ZOrderSet(hdl_all(ii,1), 0); % Setze Linie ganz nach unten (damit i.O. immer oben sind)
  end
  % Farbkarte muss danach wieder nach unten gesetzt werden
  ZOrderSet(Hdl_all.surf, 0);
  % Put text for state interval
  plot(i_stage1, 180/pi*phi_range, 'k_', 'MarkerSize', 8);
  plot(i_stage1, 180/pi*(phi_range+delta_phi/2), 'k_', 'MarkerSize', 5);
  for kk = 1:length(phi_range)
    if 180/pi*phi_range(kk) < -50 || 180/pi*phi_range(kk) > 100
      continue % nur Zustände einzeichnen, die im Plot-Bereich liegen
    end
    if usr_overlapmode == 0 % Zustände zu eng. Zu wenig Platz für volle Beschriftung
      text(i_stage1+0.06, 180/pi*phi_range(kk)-2, sprintf('%d', kk), ...
        'interpreter', 'latex');
    else % Volle Beschriftung
      text(i_stage1+0.08, 180/pi*phi_range(kk)-20, sprintf('$[x_{\\mathrm{ref},%d}]$', kk), ...
        'Rotation', 90, 'interpreter', 'latex'); % Legende für rechte Achse
    end
  end
  if usr_overlapmode == 0 % Achsbeschriftung rechts
%     text(i_stage1+0.07, 180/pi*phi_range(8)-35, '$[x_{\mathrm{ref},8}]$', ...
%       'Rotation', 90, 'interpreter', 'latex'); % Legende für rechte Achse
  end
  quiver(i_stage1*ones(1,length(phi_range)), 180/pi*phi_range, ...
    zeros(1,length(phi_range)),  180/pi*0.9*delta_phi/2*ones(1,length(phi_range)), 'off', 'k.', 'LineWidth', 1);
  quiver(i_stage1*ones(1,length(phi_range)), 180/pi*phi_range, ...
    zeros(1,length(phi_range)), -180/pi*0.9*delta_phi/2*ones(1,length(phi_range)), 'off', 'k.', 'LineWidth', 1);

  xlim(minmax2(s_stage'));
  ylim([-70, 140]);
%   if i_stage1 ~= 2
%     xlabel('');
%   end
  if i_stage1 > 1
    ylabel('');
  else
    ylabel('Redundant coordinate $\varphi_z$ in deg', 'interpreter', 'latex');
  end
  if i_stage1 < length(axhdl)
    delete(colorbar());
  end
  title(sprintf('\\textbf{(%s)} stage %d', char(96+i_stage1), i_stage1), ...
    'interpreter', 'latex');
  xlabel('Norm. traj. progress $s$', 'interpreter', 'latex');
end

set(Hdl_all.cb, 'Position', [0.90, 0.05, 0.02, 0.85]); % put cb to the right
cyh=ylabel(Hdl_all.cb, 'Performance criterion $h$ (cond.)', 'Rotation', 90, 'interpreter', 'latex');
set(cyh, 'Position', [3.6, 500, 0]); % put cblabel to the left

figure_format_publication(pmfhdl);

remove_inner_labels(axhdl, 2);
set_size_plot_subplot(pmfhdl, ...
  12.2,6,axhdl,... % 12.2 according to llncs.cls
  0.08,0.15,0.16,0.16,... %l r u d
  0.06,0) % x y
drawnow();
% Legende
I_vmactive = [2 4]; % Manuelle Auswahl der aktiven Marker. Referenz: s_pmp.violation_markers
if usr_overlapmode == 0
  LegHdl = [DP_hdl(1:3); VM_hdl(I_vmactive)];
  LegLbl = ['invalid', 'valid', 'optimal', s_pmp.violation_markers(1,I_vmactive)];
  LegLbl{strcmp(LegLbl,'jac_cond')} = 'singularity';
  LegLbl{strcmp(LegLbl,'coll_hyp')} = 'collision';
else
  LegHdl = [DP_hdl; VM_hdl(I_vmactive)];
  LegLbl = ['invalid', 'valid', 'optimal', 'add. overlap', s_pmp.violation_markers(1,I_vmactive)];
  LegLbl{strcmp(LegLbl,'jac_cond')} = 'sing.'; % Shorter text
  LegLbl{strcmp(LegLbl,'coll_hyp')} = 'coll.';
end
% LegLbl{strcmp(LegLbl,'qlim_hyp')} = 'Joint Limit';

legendflex(LegHdl, LegLbl, 'anchor', {'n','n'}, ...
  'ref', pmfhdl, ... % an Figure ausrichten (mitten oben)
  'buffer', [-4 -1], ... % Kein Versatz notwendig, da mittig oben
  'ncol', 0, 'nrow', 1, ... % eine Zeile für Legende
  'fontsize', 8, ...
  'xscale', 0.6, ... % Kleine Symbole
  ... 'padding', [-3,-3,3], ... % Leerraum reduzieren
  'box', 'on');
pdfname = ['dp_interv', overlapstr, '_stage1_to_2'];
exportgraphics(pmfhdl, fullfile(paperfig_path, [pdfname, '.pdf']),'ContentType','vector');
fprintf('Bild gespeichert: %s\n', pdfname);
cd(paperfig_path);
ghostscript(['-dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 ', ...
  '-dPDFSETTINGS=/prepress -sOutputFile=',pdfname,'_compressed.pdf ',pdfname,'.pdf']);

%% Plot complete performance map (Fig. 8)
fprintf('Zeichne gesamte Redundanzkarte für SI-DP\n');
% Debug: Load previous folder
% dpres_dir = fullfile(data_path, '..', 'LNEE_DP_debug_red');
settings_perfmap_complete = settings_perfmap;
settings_perfmap_complete.markermindist = [0.2, 15]; % Größerer Markerabstand, da kleiner im Bild
pmfhdl = change_current_figure(3); clf; hold on;
set(pmfhdl, 'Name', ['SIDP_Perfmap_All',overlapstr], 'NumberTitle', 'off');
DP_hdl = NaN(3,1); % Handle für die verschiedenen Linien (für Legende)
t1 = tic();
[Hdl_all, s_pmp] = RP.perfmap_plot(H_all, phiz_range, s_ref, settings_perfmap_complete);
fprintf('Redundanzkarte gezeichnet: %1.1fs\n', toc(t1));
VM_hdl = Hdl_all.VM;
% Plot all valid transfers from stage to stage
for i_stage1 = 1:length(DP_settings.IE)-1
  d_state = load(fullfile(dpres_dir, sprintf('dp_stage%d_final.mat', i_stage1)));
  tfn = dir(fullfile(dpres_dir, sprintf('dp_stage%d_state*_to*_result.mat', i_stage1)));
  s_stage = DP_settings.PM_s_tref(DP_settings.IE(i_stage1):DP_settings.IE(i_stage1+1));
  % Gehe alle Transfers durch, die gespeichert sind
  for ii = 1:length(tfn)
    d_ii = load(fullfile(dpres_dir, tfn(ii).name));
    i_state1 = d_ii.k;
    i_state2 = d_ii.l;
    if isinf(d_state.F_stage(i_state1, i_state2))
      % Nicht zeichnen wenn Transfer nicht erfolgreich (sonst unübersichtlich)
      continue
    end
    % Abtastung der Datenpunkte für Plot spärlicher
    Iplot = select_plot_indices_downsample_nonuniform(...
      s_stage(1:length(d_ii.X6_traj)), d_ii.X6_traj, 0.05, 3*pi/180);
    hdl = plot(s_stage(Iplot), 180/pi*d_ii.X6_traj(Iplot), 'c-', 'LineWidth', 1);
    if usr_stageoptmode && i_state2 > length(phi_range)
      set(hdl, 'LineStyle', '--'); % Teil der neu optimierten Linien
      DP_hdl(4) = hdl;
    end
    if d_final.I_best(i_stage1) == i_state1 && ...
        d_final.I_best(i_stage1+1) == i_state2
      set(hdl, 'Color', 'b');
      set(hdl, 'LineWidth', 1.5);
      DP_hdl(3) = hdl; % optimal transmission (globally)
    elseif i_state2 <= size(d_final.I_all,2) &&  d_final.I_all(i_stage1+1,i_state2) == i_state1
      DP_hdl(2) = hdl; % optimal transmission (to this stage)
      set(DP_hdl(2), 'Color', 'm');
    else
      DP_hdl(1) = hdl;
    end
  end
end
% Plot nachbearbeiten
axch = get(gca, 'children');
% Setze alle nicht-optimalen Linien nach hinten
for i = 1:length(axch)
  if strcmp(get(axch(i),'type'), 'line') && all(get(axch(i), 'color') == [0 1 1]) % cyan
    ZOrderSet(axch(i), 0);
  end
end
% Setze die global optimale Linie ganz nach vorne
for i = 1:length(axch)
  if strcmp(get(axch(i),'type'), 'line') && all(get(axch(i), 'color') == [0 0 1]) % blau
    ZOrderSet(axch(i), 1);
  end
end
ZOrderSet(Hdl_all.surf,0); % Farbkarte wieder nach hinten
% Hier nicht die neu berechnete Trajektorie, sondern die gestückelte aus
% den einzelnen Teilschritten. Sonst eventuell verwirrend/angreifbar.
% Daher folgenden Zeilen auskommentiert:
% d_output = load(fullfile(dpres_dir, 'dp_output.mat'), 'TrajDetail');
% Iplot = select_plot_indices_downsample_nonuniform(...
%   DP_settings.PM_s_tref, d_output.TrajDetail.X6, 0.05, 3*pi/180);
% DP_hdl(3) = plot(DP_settings.PM_s_tref, 180/pi*d_output.TrajDetail.X6, 'b-', 'LineWidth', 2);
xlabel('Normalized trajectory progress $s$', 'interpreter', 'latex');
ylabel('Redundant coordinate $\varphi_z$ in deg', 'interpreter', 'latex');
xlim(minmax2(DP_settings.PM_s_tref'));
ylim([-70, 140]);
set(Hdl_all.cb, 'Position', [0.90, 0.1, 0.02, 0.85]); % put cb to the right
cyh=ylabel(Hdl_all.cb, 'Performance criterion $h$', 'Rotation', 90, 'interpreter', 'latex');
set(cyh, 'Position', [3.3, 500, 0]); % put cblabel to the left
figure_format_publication(pmfhdl);
set_size_plot_subplot(pmfhdl, ...
  12.2,5,gca,... % 12.2 according to llncs.cls
  0.08,0.12,0.12,0.16,... %l r u d
  0.00,0) % x y
drawnow();
% Legende
I_vmactive = [2 4]; % Manuelle Auswahl der aktiven Marker. Referenz: s_pmp.violation_markers
LegHdl = [DP_hdl; VM_hdl(I_vmactive)];
LegLbl = ['valid', 'stage opt.', 'global opt.',s_pmp.violation_markers(1,I_vmactive)];
% LegLbl{strcmp(LegLbl,'qlim_hyp')} = 'joint limit';
LegLbl{strcmp(LegLbl,'jac_cond')} = 'singularity';
LegLbl{strcmp(LegLbl,'coll_hyp')} = 'collision';
legendflex(LegHdl, LegLbl, 'anchor', {'n','n'}, ...
  'ref', pmfhdl, ... % an Figure ausrichten (mitten oben)
  'buffer', [0 -1], ... % Kein Versatz notwendig, da mittig oben
  'ncol', 0, 'nrow', 1, ... % eine Zeile für Legende
  'fontsize', 8, ...
  'xscale', 0.4, ... % Kleine Symbole
  ...'padding', [1,1,2], ... % Leerraum reduzieren
  'box', 'on');
pdfname = ['dp_interv', overlapstr, '_result'];
exportgraphics(pmfhdl, fullfile(paperfig_path, [pdfname,'.pdf']),'ContentType','vector'); % ,'Resolution','100'
fprintf('Bild gespeichert: %s\n', pdfname);
cd(paperfig_path);
ghostscript(['-dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 ', ...
  '-dPDFSETTINGS=/prepress -sOutputFile=',pdfname,'_compressed.pdf ',pdfname,'.pdf']);
