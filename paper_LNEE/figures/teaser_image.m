% Erzeuge das kleine Vorschaubild auf der ersten Seite des Papers

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-02
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear
paperfig_path = fileparts(which('teaser_image.m')); % this file's directory
assert(~isempty(paperfig_path), 'directory of this file has to be on the PATH');
data_path = fullfile(paperfig_path, '..', '..', 'case_study', 'data_LNEE');
% Lade Bild der Redundanzkarte und speichere erneut mit neuer Formatierung

uiopen(fullfile(data_path, ['nullspace_traj1','.fig']), 1);
fighdl = gcf();
fch = get(fighdl, 'children');
cbh = fch(strcmp(get(fch, 'type'), 'colorbar')); % vor dem Löschen machen
delete(fch(strcmp(get(fch, 'type'), 'legend')));
delete(fch(3)); % legendflex-Objekt
fch = get(fighdl, 'children');

linhdl = get(gca, 'children');
% TODO: Alle Linien löschen, bis auf eine gute

% Wird in InkScape eingefügt
set(gca, 'xtick', 0:7);
xlabel(''); % Trajectory progress');
ylabel(''); % Redundant coordgiinate');
% Colorbar neu formatieren
set(cbh, 'Orientation', 'Horizontal');
set(cbh, 'Position', [0.6 0.87 0.35 0.05]);
set(cbh, 'TickLabels', {});
cbylh = get(cbh, 'ylabel');
set(cbylh, 'String', '');
set_size_plot_subplot(fighdl, ...
  7,5,gca,... % 12.2 according to llncs.cls
  0.12,0.02,0.15,0.20,... %l r u d
  0,0) % x y

saveas(fighdl, fullfile(data_path, ['perfmap_teaser','.fig']));
exportgraphics(fighdl, fullfile(paperfig_path, ['perfmap_teaser','.pdf']),'ContentType','vector') % ,'Resolution','100'
cd(paperfig_path);
ghostscript(['-dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 ', ...
  '-dPDFSETTINGS=/prepress -sOutputFile=','perfmap_teaser','_compressed.pdf ','perfmap_teaser','.pdf']);