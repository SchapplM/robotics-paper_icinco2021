% Fitness function for gain optimization via particle swarm optimization
% 
% Input:
% p
%   optimization parameters: Gains of the controller to be optimized
% optstruct
%   manually created options structure for this function
% options
%   option set for particle swarm
% 
% Output:
% f
%   fitness value. Degree of convergence to given ground truth final value
% resultstats_output
%   Statistics. Stored intermediate variables using persistent variables

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2021-05
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function [f, resultstats_output] = fitness_gainopt(p, optstruct, options)
t1 = tic();
% Initialize persistent variables to store results (not possible in normal
% particle swarm optimization
persistent resultstats
f = NaN;
size_data = [options.MaxIterations+1, options.SwarmSize];
if isempty(resultstats)
  resultstats = struct( ...
    'f', NaN(size_data), ...
    'p', NaN(size_data(2), length(p), size_data(1)));
end
if nargout == 2
  resultstats_output = resultstats;
  return
end
while true % dummy loop to be able to break
  if all(p==0)
    % alle Parameter sind Null. Das führt dazu, dass in der
    % Trajektorien-Funktion keine Leistungsmerkmale berechnet werden.
    f = inf;
    if optstruct.verbose
      fprintf('fitness_gainopt: all parameters are zero\n');
    end
    break;
  end
  if all(p(2:3)==0)
    % alle h-Parameter sind Null. Das führt dazu, dass keine Bewegung
    % stattfindet. Es gibt eine Dämpfung, aber keine zu dämpfende Bewegung.
    f = 1e8;
    if optstruct.verbose
      fprintf('fitness_gainopt: only damping parameter is non-zero\n');
    end
    break;
  end
  % Initialize robot
  RP = optstruct.RP;
  % Ground truth value should be obtained via position-level kinematics
  % (which is more robust against oscillations)
  hopt_groundtruth = optstruct.hopt_groundtruth;
  % Trajectory
  Traj_X = optstruct.Traj_X;
  Traj_XD = optstruct.Traj_XD;
  Traj_XDD = optstruct.Traj_XDD;
  Traj_t = optstruct.Traj_t;
  qs = optstruct.qs; % Starting joint configuration

  % initialize weighting factors for the IK objectives
  wn_traj = zeros(8,1);
  wn_traj([3 6 10]) = p(1:3); % Kv, Kp (cond), Kd (cond)
  % hyperbolische Abstoßung von Grenzen
  wn_traj([2 8]) = p(4:5); % Kp and Kd (hyperbolic limit avoidance)

  s_traj_ii = struct('wn', wn_traj);
  s_traj_ii.thresh_ns_qa = optstruct.thresh_ns_qa;
  %% Calculate Trajektory
  [Q_ii, QD_ii, QDD_ii, Phi_ii,~,~,~,Stats_traj] = RP.invkin2_traj( ...
    Traj_X, Traj_XD, Traj_XDD, Traj_t, qs, s_traj_ii);
  %% Gütefunktion bestimmen
  % Bestimme normierte Überschreitung der Grenzen
  qlim = cat(1,RP.Leg(:).qlim);
  qDlim = cat(1,RP.Leg(:).qDlim);
  qDDlim = cat(1,RP.Leg(:).qDDlim)*1.1; % Numerische Toleranz berücksichtigen
  Q_ii_norm = (Q_ii-repmat(qlim(:,1)', size(Q_ii,1),1)) ./ ...
              repmat(qlim(:,2)'-qlim(:,1)', size(Q_ii,1),1);
  QD_ii_norm = (QD_ii-repmat(qDlim(:,1)', size(QD_ii,1),1)) ./ ...
              repmat(qDlim(:,2)'-qDlim(:,1)', size(QD_ii,1),1);
  QDD_ii_norm = (QDD_ii-repmat(qDDlim(:,1)', size(QDD_ii,1),1)) ./ ...
              repmat(qDDlim(:,2)'-qDDlim(:,1)', size(QDD_ii,1),1);
  [maxviol_q_norm, I_maxviol] = max(max([Q_ii_norm-1, -Q_ii_norm],[],2));
  % Test: [max(Q_ii_norm(:)), min(Q_ii_norm(:))]
  % Test: find(Q_ii_norm(I_maxviol,:)>1|Q_ii_norm(I_maxviol,:)<0)
  if maxviol_q_norm > 0
    % normiere auf 8e3 ... 9e3
    f = 1e3*(8+2/pi*atan(maxviol_q_norm));
    if optstruct.verbose
      fprintf(['fitness_gainopt: position limits exceeded by %1.1f (norm.). ', ...
        'Duration: %1.1fs\n'], maxviol_q_norm, toc(t1));
    end
    break
  end
  maxviol_qD_norm = max([QD_ii_norm(:)-1; -QD_ii_norm(:)]);
  % Test: [max(QD_ii_norm(:)), min(QD_ii_norm(:))]
  if maxviol_qD_norm > 0.1 % kann numerisch nict exakt eingehalten werden
    % normiere auf 7e3 ... 8e3
    f = 1e3*(7+2/pi*atan(maxviol_qD_norm));
    if optstruct.verbose
      fprintf(['fitness_gainopt: velocity limits exceeded by %1.1f (norm.)\n', ...
        'Duration: %1.1fs\n'], maxviol_qD_norm, toc(t1));
    end
    break
  end
  maxviol_qDD_norm = max([QDD_ii_norm(:)-1; -QDD_ii_norm(:)]);
  % Test: [max(QDD_ii_norm(:)), min(QDD_ii_norm(:))]
  if maxviol_qDD_norm > 0.1
    % normiere auf 6e3 ... 7e3
    f = 1e3*(6+2/pi*atan(maxviol_qDD_norm));
    if optstruct.verbose
      fprintf(['fitness_gainopt: acceleration limits exceeded by %1.1f (norm.)\n', ...
        'Duration: %1.1fs\n'], maxviol_qDD_norm, toc(t1));
    end
    break
  end
  
  finalerror = Stats_traj.h(end,1+6) - hopt_groundtruth;
  RMSE = sqrt(mean((Stats_traj.h(:,1+6) - hopt_groundtruth).^2));
  if optstruct.verbose
    fprintf(['fitness_gainopt: success: condition %1.1f (vs %1.1f ground ', ...
      'truth), rmse=%1.1f. Duration: %1.1fs\n'], Stats_traj.h(end,1+6), ...
      hopt_groundtruth, RMSE, toc(t1));
  end
  % Normiere Ergebnis auf 0..1e3
  % determine quality of convergence regarding proximity to ground truth
  % result and RMSE denoting oscillations
  f = 1e3* 2/pi*atan((finalerror + RMSE)/100);
  break % end dummy loop
end
%% Save output
data_transp = resultstats.f';
k=find(isnan(data_transp(:)), 1, 'first'); % 1D-Index in Matrix
if isempty(k)
  fprintf('Keine Speicherung des Ergebnisses möglich (nicht innerhalb des PSO)\n');
else
  [j,i] = ind2sub(fliplr(size_data),k); % Umrechnung in 2D-Indizes. i=Generation, j=Individuum
  resultstats.f(i,j) = f;
  resultstats.p(j,:,i) = p;
end
if isinf(f)
  % Kein Debuggen möglich
  return
end
if ~optstruct.debug
  return
end
%% Debug
figure(200); clf;
subplot(2,3,1);
plot(Traj_t, Q_ii_norm);
grid on; ylabel('q (norm)');
subplot(2,3,2);
plot(Traj_t, QD_ii_norm);
grid on; ylabel('qD (norm)');
subplot(2,3,3);
plot(Traj_t, QDD_ii_norm);
grid on; ylabel('qDD (norm)');
subplot(2,3,4); hold all
plot(Traj_t, Stats_traj.h(:,1+6));
plot(Traj_t([1 end]), hopt_groundtruth*[1;1], 'r--');
grid on; ylabel('h=cond');
subplot(2,3,5); hold all
plot(Traj_t, Stats_traj.h(:,1+2));
grid on; ylabel('h=hyperbol. limit');
subplot(2,3,6); hold all
plot(Traj_t, Stats_traj.h(:,1));
plot(Traj_t, wn_traj(6)*Stats_traj.h(:,1+6), '--');
grid on; ylabel('h sum');
legend({'h sum', 'wn*h (cond)'});
linkxaxes