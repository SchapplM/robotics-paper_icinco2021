% Run all scripts related to the LNEE paper to reproduce data and figures

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-02
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

this_path = fileparts( mfilename('fullpath') );
addpath(this_path);
run(fullfile(fileparts(which('run_evaluation_LNEE.m')), 'case_study', ...
  'ParRob_dynprog.m')); close all;
run(fullfile(fileparts(which('run_evaluation_LNEE.m')), 'paper_LNEE', ...
  'figures', 'fig_dynprog_discr_example.m')); close all;
run(fullfile(fileparts(which('run_evaluation_LNEE.m')), 'paper_LNEE', ...
  'figures', 'fig_dynprog_interv_example.m')); close all;