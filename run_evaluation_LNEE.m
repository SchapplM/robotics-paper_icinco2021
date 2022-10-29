% Run all scripts related to the LNEE paper to reproduce data and figures

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-02
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

this_path = fileparts( mfilename('fullpath') );
addpath(this_path);
run(fullfile(fileparts(which('run_evaluation_LNEE.m')), 'case_study', ...
  'ParRob_dynprog.m')); close all;

%% Auswertung für diskreten Modus
matlabscript_discr = fullfile(fileparts(which('run_evaluation_LNEE.m')), 'paper_LNEE', ...
  'figures', 'fig_dynprog_discr_example.m');
% Schalte den Modus der Datei um, um direkt beide Bilder zu generieren
f = strrep(fileread(matlabscript_discr), 'usr_stageoptmode = false;', ...
  'usr_stageoptmode = true;');
matlabscript_discr2 = strrep(matlabscript_discr, '.m', '_copy_modechange.m');
fid  = fopen(matlabscript_discr2,'w'); fprintf(fid,'%s',f); fclose(fid);
run(fullfile(fileparts(which('run_evaluation_LNEE.m')), 'paper_LNEE', ...
  'figures', 'fig_dynprog_discr_example.m')); close all;
run(fullfile(fileparts(which('run_evaluation_LNEE.m')), 'paper_LNEE', ...
  'figures', 'fig_dynprog_discr_example_copy_modechange.m')); close all;
%% Auswertung für Intervall-Modus
matlabscript_interv = fullfile(fileparts(which('run_evaluation_LNEE.m')), 'paper_LNEE', ...
  'figures', 'fig_dynprog_interv_example.m');
% Schalte den Modus der Datei um, um direkt beide Bilder zu generieren
f = strrep(fileread(matlabscript_interv), 'usr_overlapmode = false;', ...
  'usr_overlapmode = true;');
matlabscript_interv2 = strrep(matlabscript_interv, '.m', '_copy_modechange.m');
fid  = fopen(matlabscript_interv2,'w'); fprintf(fid,'%s',f); fclose(fid);
run(fullfile(fileparts(which('run_evaluation_LNEE.m')), 'paper_LNEE', ...
  'figures', 'fig_dynprog_interv_example.m')); close all;
run(fullfile(fileparts(which('run_evaluation_LNEE.m')), 'paper_LNEE', ...
  'figures', 'fig_dynprog_interv_example_copy_modechange.m')); close all;