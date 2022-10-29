% Example Figure for Discrete Dynamic Programming (Sect. 4.1 of LNEE paper)
% This creates Fig. 3 and 4 of the paper
% 
% Preliminaries:
% Run case_study/ParRob_dynprog.m to create necessary files.
% 
% 
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-02
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear
% User settings
usr_highres_distrfig = true; % high resolution of the paper figure for performance criterion map
usr_stageoptmode = false;
%% Initialize
% Initialize Robot
% Default output directory (for paper)
paperfig_path = fileparts(which('fig_dynprog_discr_example.m'));
this_dir = fileparts(which('fig_dynprog_discr_example.m'));
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
if ~usr_stageoptmode
  filename_post = ['_dynprog_costRMStraj_n9_', 'nored', '', '.mat'];
else
  filename_post = ['_dynprog_costRMStraj_n9_', 'nored', '_stageopt', '.mat'];
end
filename_dynprog= fullfile(data_path, [filename_pre, filename_post]);
assert(exist(filename_dynprog, 'file'), 'dynamic programming results file does not exist');
d = load(filename_dynprog);
DP_XE = d.DP_XE;
DP_Stats = d.DP_Stats;
DP_TrajDetail = d.DP_TrajDetail;
DP_settings = d.DP_settings;
stageoptstr = '';
if ~usr_stageoptmode 
  debugfoldername_full = 'LNEE_Traj1_DP_debug_costRMStraj_n9_nored';
else
  debugfoldername_full = 'LNEE_Traj1_DP_debug_costRMStraj_n9_nored_stageopt';
  stageoptstr= '_stageopt';
end
dpres_dir = fullfile(this_dir, '..', '..', 'case_study', debugfoldername_full);
assert(exist(dpres_dir, 'file'), 'directory with debug information for DP does not exist');

% Eigenschaften der Intervalle aus Ergebnis ausgeben lassen. Siehe dynprog_taskred_ik
delta_phi = DP_Stats.delta_phi;
phi_range = DP_Stats.phi_range;
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
d_final = load(fullfile(dpres_dir, 'dp_final.mat'), 'I_all', 'I_best');
%% Plot performance map with three first stages (Figure 3 of the paper)
fprintf('Zeichne Bilder für DP Einzeltransfer\n');
pmfhdl = change_current_figure(1); clf;
axhdl = NaN(1,3);
DP_hdl = NaN(3,1); % Handle für die verschiedenen Linien (für Legende)
for i_stage1 = 1:3
  s_stage = DP_settings.PM_s_tref(DP_settings.IE(i_stage1):DP_settings.IE(i_stage1+1));
  axhdl(i_stage1) = subplot(1,3,i_stage1); hold on;
  settings_perfmap_stage = settings_perfmap;
  settings_perfmap_stage.s_range_plot = minmax2(s_stage') + [-0.1, 0.1];
  t1 = tic();
  [Hdl_all, s_pmp] = RP.perfmap_plot(H_all, phiz_range, s_ref, settings_perfmap_stage);
  fprintf('Redundanzkarte gezeichnet: %1.1fs\n', toc(t1));
  VM_hdl = Hdl_all.VM;
  % Plot lines for first stage
%   i_stage1 = 1;
  d_state = load(fullfile(dpres_dir, sprintf('dp_stage%d_final.mat', i_stage1)));
  tfn = dir(fullfile(dpres_dir, sprintf('dp_stage%d_state*_to*_result.mat', i_stage1)));
  % Gehe alle Transfers durch, die gespeichert sind
  hdl_all = NaN(length(tfn), 2); % Erste Spalte Handle, zweite Spalte Kategorie als Zahl 
  for ii = 1:length(tfn)
    d_ii = load(fullfile(dpres_dir, tfn(ii).name));
    i_state1 = d_ii.k;
    i_state2 = d_ii.l;
    % Abtastung der Datenpunkte für Plot spärlicher, damit Dateigröße gering
    % bleibt. Abstand wie hier 5% bzw. 3° macht Änderung von 1,6MB zu 100kB
    Iplot = select_plot_indices_downsample_nonuniform(...
      s_stage(1:length(d_ii.X6_traj)), d_ii.X6_traj, 0.05, 3*pi/180);
%     hdl = plot(s_stage(1:length(d_i.X6_traj)), 180/pi*d_i.X6_traj, 'k--'); % TODO: _ref
    hdl_all(ii,1) = plot(s_stage(Iplot), 180/pi*d_ii.X6_traj(Iplot), 'k-');
    set(hdl_all(ii,1), 'LineWidth', 1);
    if usr_stageoptmode && i_state2 > size(d_final.I_all,2)
      set(hdl_all(ii,1), 'LineStyle', '--'); % Teil der neu optimierten Linien
    end
    % Prüfe ob es sich um die optimale Teil-Politik handelt
    if i_state2 <= size(d_final.I_all,2) && d_final.I_all(i_stage1+1,i_state2) == i_state1
      hdl_all(ii,2) = 1;
      set(hdl_all(ii,1), 'Color', 'm');
      DP_hdl(3) = hdl_all(ii,1); % optimal transmission (to this stage)
    elseif ~isinf(d_state.F_stage(i_state1, i_state2))
      hdl_all(ii,2) = 2;
      set(hdl_all(ii,1), 'Color', 'c');
      DP_hdl(2) = hdl_all(ii,1); % line for valid transition
    else
      hdl_all(ii,2) = 3;
      if i_state2 <= size(d_final.I_all,2)
        DP_hdl(1) = hdl_all(ii,1); % line for invalid transition
      else
        DP_hdl(4) = hdl_all(ii,1); % for stage-optimization legend entry
      end
    end
  end
  if isnan(DP_hdl(2))
    % Es gibt keine "valid"-Linie, da alle Lösungen stage-optimal sind
    DP_hdl(2) = plot(NaN, NaN, 'c-', 'LineWidth', 1);
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
  xlim(minmax2(s_stage'));
  ylim([-70, 130]);

  % Put text for state interval
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
  if i_stage1 ~= 2
    xlabel('');
  else
    xlabel('Normalized trajectory progress $s$', 'interpreter', 'latex');
  end
end
set(Hdl_all.cb, 'Position', [0.90, 0.03, 0.02, 0.95]); % put cb to the right
cyh=ylabel(Hdl_all.cb, 'Performance criterion $h$ (cond.)', 'Rotation', 90, 'interpreter', 'latex');
set(cyh, 'Position', [3.5, 500, 0]); % put cblabel to the left and down
figure_format_publication(pmfhdl);

remove_inner_labels(axhdl, 2);
set_size_plot_subplot(pmfhdl, ...
  12.2,5,axhdl,... % 12.2 according to llncs.cls
  0.08,0.12,0.19,0.16,... %l r u d
  0.03,0) % x y
drawnow();
% Legende
I_vmactive = [2 4]; % Manuelle Auswahl der aktiven Marker. Referenz: s_pmp.violation_markers
if usr_stageoptmode == 0
  LegHdl = [DP_hdl(1:3); VM_hdl(I_vmactive)];
  LegLbl = ['invalid', 'valid', 'optimal', s_pmp.violation_markers(1,I_vmactive)];
else
  LegHdl = [DP_hdl; VM_hdl(I_vmactive)];
  LegLbl = ['invalid', 'valid', 'optimal', 'add. stage opt.', s_pmp.violation_markers(1,I_vmactive)];
end
% LegLbl{strcmp(LegLbl,'qlim_hyp')} = 'Joint Limit';
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
pdfname = sprintf('dp_discr%s_stage1_to_3', stageoptstr);
exportgraphics(pmfhdl, fullfile(paperfig_path, [pdfname, '.pdf']),'ContentType','vector') 
fprintf('Bild gespeichert: %s\n', pdfname);
% Zur Kompression des Bildes (bringt ca. 30% bei hoher Auflösung):
cd(paperfig_path);
ghostscript(['-dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 ', ...
  '-dPDFSETTINGS=/prepress -sOutputFile=',pdfname,'_compressed.pdf ',pdfname,'.pdf']);

%% Plot complete performance map (Figure 4 of the paper)
fprintf('Zeichne gesamte Redundanzkarte für DP\n');
settings_perfmap_complete = settings_perfmap;
settings_perfmap_complete.markermindist = [0.2, 15]; % Größerer Markerabstand, da kleiner im Bild
pmfhdl = change_current_figure(3); clf; hold on;
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
% Nehme nicht die global-optimale Trajektorie aus Ergebnis-Datei, sondern
% die von oben (sollte keinen Unterschied darstellen).
% d_output = load(fullfile(dpres_dir, 'dp_output.mat'), 'TrajDetail');
% Iplot = select_plot_indices_downsample_nonuniform(...
%   DP_settings.PM_s_tref, d_output.TrajDetail.X6, 0.05, 3*pi/180);
% DP_hdl(3) = plot(DP_settings.PM_s_tref, 180/pi*d_output.TrajDetail.X6, 'b-', 'LineWidth', 2);
xlabel('Normalized trajectory progress $s$', 'interpreter', 'latex');
ylabel('Redundant coordinate $\varphi_z$ in deg', 'interpreter', 'latex');
xlim(minmax2(DP_settings.PM_s_tref'));
ylim([-70, 130]);
set(Hdl_all.cb, 'Position', [0.90, 0.03, 0.02, 0.9]); % put cb to the right
cyh=ylabel(Hdl_all.cb, 'Performance criterion $h$ (cond.)', 'Rotation', 90, 'interpreter', 'latex');
set(cyh, 'Position', [3.3, 500, 0]); % put cblabel to the left
figure_format_publication(pmfhdl);
set_size_plot_subplot(pmfhdl, ...
  12.2,5,gca,... % 12.2 according to llncs.cls
  0.09,0.12,0.14,0.16,... %l r u d
  0.00,0) % x y
drawnow();
% Legende
I_vmactive = [2 4]; % Manuelle Auswahl der aktiven Marker. Referenz: s_pmp.violation_markers
LegHdl = [DP_hdl; VM_hdl(I_vmactive)];
LegLbl = ['valid', 'stage opt.', 'global opt.', s_pmp.violation_markers(1,I_vmactive)];
% LegLbl{strcmp(LegLbl,'qlim_hyp')} = 'joint limit';
LegLbl{strcmp(LegLbl,'jac_cond')} = 'singularity';
LegLbl{strcmp(LegLbl,'coll_hyp')} = 'collision';
legendflex(LegHdl, LegLbl, 'anchor', {'n','n'}, ...
  'ref', pmfhdl, ... % an Figure ausrichten (mitten oben)
  'buffer', [-3 -1], ... % Kein Versatz notwendig, da mittig oben
  'ncol', 0, 'nrow', 1, ... % eine Zeile für Legende
  'fontsize', 8, ...
  'xscale', 0.5, ... % Kleine Symbole
  'padding', [1,1,0], ... % Leerraum reduzieren
  'box', 'on');
pdfname = sprintf('dp_discr%s_result', stageoptstr);
% TODO: Beide Export-Befehle führen bei Windows-Laptop zu Fehler.
exportgraphics(pmfhdl, fullfile(paperfig_path, [pdfname, '.pdf']),'ContentType','vector') % ,'Resolution','100'
export_fig(pmfhdl, fullfile(paperfig_path, [pdfname, '2.pdf']));
fprintf('Bild gespeichert: %s\n', pdfname);
cd(paperfig_path);
ghostscript(['-dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 ', ...
  '-dPDFSETTINGS=/prepress -sOutputFile=',pdfname,'_compressed.pdf ',pdfname,'.pdf']);
