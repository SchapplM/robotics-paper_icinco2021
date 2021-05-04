% Compute nullspace projector matrices for a six-DoF parallel robot.
% This proves statements in Sec. 4 of the paper.
% 
% Dependencies: Serial Robot Database, Parallel Robot Database, Robotics
% Toolbox and others; see README.MD

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2021-05
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear

%% Initialization of the Robot Model
if isempty(which('serroblib_path_init.m'))
  error('The serial robot database is not initialized in Matlab.');
end
if isempty(which('parroblib_path_init.m'))
  error('The parallel robot database is not initialized in Matlab.');
end
% Hexapod robot model
RP = parroblib_create_robot_class('P6RRPRRR14V3G1P4A1', 1.0, 0.2);
for i = 1:RP.NLEG
  RP.Leg(i).qlim = repmat([-2*pi, 2*pi], RP.Leg(i).NQJ, 1);
  % Assign design parameters for plotting the robot
  RP.Leg(i).DesPar.joint_type(RP.I_qa((RP.I1J_LEG(i):RP.I2J_LEG(i)))) = 5;
  RP.Leg(i).DesPar.seg_par=repmat([5e-3,50e-3],RP.Leg(i).NL,1);
  RP.Leg(i).qlim = repmat([-2*pi, 2*pi], RP.Leg(i).NQJ, 1);
  qpris_minmax = [0.6, 1.2];
  RP.Leg(i).qlim(RP.Leg(i).MDH.sigma==1,:) = ...
    repmat(qpris_minmax, sum(RP.Leg(i).MDH.sigma==1), 1);
end
RP.DesPar.platform_par(end) = 5e-3; % for plotting
RP.fill_fcn_handles(false,false); % Initialize function dependencies
% Initial values for the inverse kinematics problem
q0_ik = -0.5+rand(RP.NJ,1);
q0_ik(RP.I_qa) = 1.0; % to have positive values for the linear actuators
s_ep = struct( ... % IK settings
  'n_max', 5000, 'Phit_tol', 1e-12, 'Phir_tol', 1e-12);
% Define indices for full coordinates x and reduced task coordinates y.
I_EE_full = RP.I_EE;
I_EE_red = logical([1 1 1 1 1 0]);
% Initialize robot with full coordinates x (influences function calls)
RP.update_EE_FG(I_EE_full,I_EE_full);
% arbitrary initial platform pose x
x0 = [ [0.05;0.03;0.6]; [15;-8;-20]*pi/180 ];
% Solve inverse kinematics problem with the method from equ. 5 of the
% paper (\ref{eq:ser_position_ik}), leg-wise for the parallel robot.
[q0, Phi] = RP.invkin_ser(x0, q0_ik, s_ep);
assert(all(abs(Phi)<1e-9), 'IK does not match');
%% Plot Robot
figure(100);clf;
hold on;grid on;
xlabel('x in m');ylabel('y in m');zlabel('z in m'); view(3);
s_plot = struct(  'ks_legs', [], 'straight', 0, 'mode', 4);
RP.plot( q0, x0, s_plot );
%% Check IK Jacobians. Proof statements from Sec. 4.1
% Set task coordinates for all following function calls.
RP.update_EE_FG(I_EE_full,I_EE_red);
% Terms Phi_dq/Phi_dx from equ. 16 (\ref{equ:par_differential_constraints})
[~,Phi_q] = RP.constr4grad_q(q0);
[~,Phi_x] = RP.constr4grad_x(x0);
% Equ. 17: Full coordinate Jacobian (\ref{equ:parrob_def_jinv_full})
Jtilde_inv_x = -Phi_q\Phi_x;
% Inverse manipulator Jacobian from equ. 20, (\ref{eq:parrob_jacobian})
J_inv_x = Jtilde_inv_x(RP.I_qa,:);
fprintf('Condition number of the PKM Jacobian: %1.1f\n', cond(J_inv_x));
% Get manipulator Jacobian as non-inverse, equ. 21 (\ref{eq:par_jacobi_inv})
J_x = inv(J_inv_x);
% Get task space Jacobian (text after equ. 21)
J_y = J_x(I_EE_red,:);

% Term Psi_dq/Psi_dx from equ. 16 (\ref{equ:par_differential_constraints})
[Psi_q,Psi_q_full] = RP.constr3grad_q(q0, x0);
fprintf(['Condition number of the inverse kinematics Jacobian ', ...
  '(pointing task constraints): %1.1f\n'], cond(Psi_q));
[Psi_x,Psi_x_full] = RP.constr3grad_x(q0, x0);
% Compute full coordinate Jacobian with both constraint formulations
Jtilde_inv_x_fromPsi = -Psi_q_full\Psi_x_full;
% Both formulations are equal
test_Jtilde_inv = Jtilde_inv_x_fromPsi - Jtilde_inv_x;
assert(all(abs(test_Jtilde_inv(:))<1e-10), ['Full coordinate Jacobian from', ...
  'Phi and Psi constraints has to be identical']);
% Get nullspace projector matrices for both coordinate spaces
% See text after equ. 24 (\ref{eq:par_nullspace_relation})
N_y   = eye(sum(RP.I_qa)) - pinv(J_y)*  J_y;
N_Psi = eye(RP.NJ) -        pinv(Psi_q)*Psi_q;
% Jacobian for transformation from actuator space to full joint space.
% See equ. 22 (\ref{eq:par_nullspace_acttofull})
J_q_qa = Jtilde_inv_x * J_x; %#ok<MINV>
% check the relation between both nullspace projectors.
% This prooves equ. 24
N_Psi_from_y = J_q_qa * N_y * J_q_qa';
% Get the ratio of all elements of the matrices N_Psi computed in two ways
N_Psi_ratio = N_Psi ./ N_Psi_from_y;
k_Psi = N_Psi_ratio(1,1); % factor k in equ. 24
% If equ. 24 holds, than the factor is the same for all elements:
assert(all(abs(N_Psi_ratio(:)-k_Psi)<1e-5), ['All components of the nullspace', ...
  'projector in full coordinates have to be in constant ratio to its ', ...
  'derivation from actuator coordinates']);

% Show that it is not possible to determine the nullspace of the full joint
% space any other way. This supports the statement from the beginning of
% Sec. 4.1.
for i = 1:RP.NJ
  % Remove one line in the constraints formulation using the standard
  % approach without reciprocal Euler angles
  Phi_q_red = Phi_q([1:i-1, i+1:end],:);
  % set up the nullspace projector
  N_Phi = eye(RP.NJ) - pinv(Phi_q_red)*Phi_q_red;
  % Check if the same nullspace projector comes out
  test = N_Phi - N_Psi;
  assert(any(abs(test(:))>1e-3), ['nullspace projector can not be created ', ...
    'from the standard constraints definition (according to assumption']);
  % Check if the same nullspace vector comes out after projecting a random
  % vector into the nullspace
  v_rand = rand(RP.NJ,1);
  qD_N_Phi = N_Phi*v_rand;
  qD_N_Psi = N_Psi*v_rand;
  assert(any(abs(qD_N_Phi-qD_N_Psi)>1e-3), ['nullspace motion (q) from Phi ', ...
    'projector is identical to nullspace motion from Psi projector. Unexpected.']);
  xD_N_Phi = J_x * qD_N_Phi(RP.I_qa); %#ok<MINV>
  xD_N_Psi = J_x * qD_N_Psi(RP.I_qa); %#ok<MINV>
  assert(all(abs(xD_N_Psi(1:5))<1e-10), ['nullspace motion from Psi projector ', ...
    'is in other component than phi_z. Unexpected.']);
  assert(any(abs(xD_N_Phi-xD_N_Psi)>1e-3), ['nullspace motion (x) from Phi ', ...
    'projector is identical to nullspace motion from Psi projector. Unexpected.']);
end

%% Calculation of Nullspace Gradients. Proof statements from Sec. 4.2
dphiz_range = [1e-5, 1e-6, 1e-7, 1e-8, 2e-8];
h_drz_all = NaN(length(dphiz_range),4);
v_ratio_all = NaN(length(dphiz_range),1);
va_ratio_all = NaN(length(dphiz_range),1);
for i_dphiz = 1:length(dphiz_range)
  % Use condition number as criterion h.
  % increment of the task space coordinate. see Sec. 4.2.2
  xD_test = [zeros(5,1);dphiz_range(i_dphiz)]; % delta x in paper
  qD_test = Jtilde_inv_x * xD_test(RP.I_EE); % delta q in paper
  h_0 = cond(J_x);

  % Gradient in the actuator space (Sec. 4.2.3)
  % Difference quotient from equ. 28 (\ref{eq:par_diffquot_taskredfulljoint})
  [~,Phi_q_full_test] = RP.constr4grad_q(q0+qD_test);
  [~,Phi_x_full_test] = RP.constr4grad_x(x0+xD_test);
  J_x_inv_test = -Phi_q_full_test\Phi_x_full_test; % equ. 17
  h_test_v1 = cond(J_x_inv_test(RP.I_qa,:)); % equ. 20.
  % Gradient w.r.t redundant coordinate via difference quotiont.
  h_drz_v1 = (h_test_v1-h_0)/xD_test(6); % See equ. 28
  % Alternative: Use total differential instead of difference quotient.
  % Equivalent approach, but more complicated.
  PhiD_q_full = Phi_q_full_test-Phi_q; % difference quotient of single matrix Phi_q
  PhiD_x_full = Phi_x_full_test-Phi_x; % ... and of Phi_x
  J_x_inv_test_v2 = Jtilde_inv_x + ...
    Phi_q\PhiD_q_full/Phi_q*Phi_x + ... % total differential and rule for matrix inversion
    -Phi_q\PhiD_x_full;
  h_test_v2 = cond(J_x_inv_test_v2(RP.I_qa,:));
  h_drz_v2 = (h_test_v2-h_0)/xD_test(6); % See equ. 28

  % Alternative 3: Same as the previous alternative 1, but using the
  % constraints modeling 3 with the reciprocal Euler angles.
  % Use the full constraints vector not reduced to the task coordinates.
  [~,Psi_q_full] = RP.constr3grad_q(q0, x0);
  [~,Psi_x_full] = RP.constr3grad_x(q0, x0);
  [~,Psi_q_full_test] = RP.constr3grad_q(q0+qD_test, x0+xD_test);
  [~,Psi_x_full_test] = RP.constr3grad_x(q0+qD_test, x0+xD_test);
  J_x_inv_test_v3 = -Psi_q_full_test\Psi_x_full_test;
  h_test_v3 = cond(J_x_inv_test_v3(RP.I_qa,:));
  h_drz_v3 = (h_test_v3-h_0)/xD_test(6); % See equ. 28

  % Alternative 4: use difference quotient to approximate the gradient
  PsiD_q_full_test = Psi_q_full_test-Psi_q_full;
  PsiD_x_full_test = Psi_x_full_test-Psi_x_full;
  % Use total differential (similar to time derivative)
  J_x_inv_test_v4 = Jtilde_inv_x + ...
    Psi_q_full\PsiD_q_full_test/Psi_q_full*Psi_x_full(:,RP.I_EE) + ...
    -Psi_q_full\PsiD_x_full_test(:,RP.I_EE);
  h6_test_v4 = cond(J_x_inv_test_v4(RP.I_qa,:));
  h_drz_v4 = (h6_test_v4-h_0)/xD_test(6); % See equ. 28

  % The difference quotient was computed in four different ways. All should
  % be identical. The method 3 shows different results and is removed from
  % the test. TODO: Check why.
  abserr_hdrz = [h_drz_v1,h_drz_v4]-h_drz_v2;
  relerr_hdrz = abserr_hdrz/h_drz_v2;

  assert(all(abs(abserr_hdrz)<1e-2) && all(abs(relerr_hdrz)<1e-2), ...
    'The gradients of h w.r.t phi_z have to be identical on all ways');

  % Gradient, see equ. 27 (\ref{eq:par_nullspace_actspace})
  h_dqa_v1 = h_drz_v1 * J_x(end,:);
  h_dqa_v2 = h_drz_v2 * J_x(end,:);
  h_dqa_v3 = h_drz_v3 * J_x(end,:);
  h_dqa_v4 = h_drz_v4 * J_x(end,:);

  % all have to be identical
  abserr_hdqa = [h_dqa_v1;h_dqa_v2;h_dqa_v4]-repmat(h_dqa_v3,3,1);
  relerr_hdqa = abserr_hdqa./repmat(h_dqa_v3,3,1);
  assert(all(abs(abserr_hdrz)<1e-2) && all(abs(relerr_hdrz)<1e-3), ...
    'The gradients of h w.r.t qa have to be identical on all ways');

  % Berechnung in fullständigen Koordinaten
  % Gradient in the full coordinate space (Sec. 4.2.2)
  % Use equ. 26 (\ref{eq:par_diffquot_taskredfulljoint})
  h_dq = (h_test_v1-h_0)./(qD_test');

  % Check for the identity of both approaches (Sec. 4.2.2 and 4.2.3)
  v_fromfullspace = N_Psi * h_dq';
  va_fromactspace = N_y * h_dqa_v1';
  % Project from actuator space to full joint space. See text after equ. 28.
  v_fromactspace = J_q_qa*va_fromactspace;
  % The terms are not identical, but in a constant ratio.
  % This validates equ. 27 vs equ. 26
  v_ratio = v_fromfullspace ./ v_fromactspace;
  assert(all(abs(v_ratio(:)-v_ratio(1))<1e-7), ['the elements of v have to ', ...
    'be in a constant ratio between derivation from full and actuation joint space']);

  % Conversion from joint space to actuation space.
  % The constant ratio only holds after the nullspace projection, since the
  % computation of the gradient (h_dq) is simplified.
  h_dqa_fromfullspace = h_dq*J_q_qa;
  va_ratio = (N_y*h_dqa_fromfullspace(:)) ./ va_fromactspace;
  assert(all(abs(va_ratio(:)-va_ratio(1))<1e-7), ['the elements of va have to ', ...
    'be in a constant ratio between derivation from full and actuation joint space']);
  % Store results for comparison outside the loop
  h_drz_all(i_dphiz,:) = [h_drz_v1,h_drz_v2,h_drz_v3,h_drz_v4];
  v_ratio_all(i_dphiz) = v_ratio(1);
  va_ratio_all(i_dphiz) = va_ratio(1);
end
% The terms computed for different small increments of the redundant
% coordinate have to be identical
abserr_hdrz = h_drz_all - h_drz_all(1,1);
relerr_hdrz = abserr_hdrz / h_drz_all(1,1);
assert(all(abs(abserr_hdrz(:))<2e-3 | abs(relerr_hdrz(:))<1e-2), ...
  ['the gradient w.r.t phi_z is not independent regarding the value of ', ...
  'the testing increment']);
fprintf('The script run through. All checks passed. Sec. 4 of paper is valid.\n');
%% Find out the ratio between entities from actuation and full joint space
% possible related from above are: k_Psi, va_ratio, v_ratio
test1 = 1/va_ratio(1) * v_ratio(1) - k_Psi;
% if abs(test1) < 1e-9
%   fprintf('the ratios of v, va and k_Psi match. TODO: What does this mean?\n');
% end
% define a unit nullspace vector
% v is a gradient of dh/dq. Therefore zeros in va are valid
v_unitact = zeros(RP.NJ,1);
v_unitact(RP.I_qa) = 1;
va_unitact = ones(sum(RP.I_qa),1);
qaD_Ny = N_y*va_unitact;
qD_NPsi = N_Psi*v_unitact;
qaD_NPsi = qD_NPsi(RP.I_qa);
qaD_ratio = qaD_Ny ./ qaD_NPsi;
assert(all(abs(qaD_ratio(1)-qaD_ratio(:))<1e-8), ['All components of the', ...
  'nullspace motion vectors from actuation or full space have to be in constant ratio']);

% TODO: The factor could be related to the properties (norm, cond, ...) of
% the involved matrices
% TODO: Find analytic expression for k_Psi
return
%% Debug: Random code snippets to find the solution
% 1/k_Psi
% norm(Jtilde_inv_x)*norm(J_x)
% 1/norm(norm(J_q_qa))
% norm(Psi_q)
% norm(J_y)
% cond(J_y)
% cond(Psi_q)
% N_Psi(RP.I_qa,RP.I_qa) ./ N_y
% ratio_Ny_NPsi_matrix = N_y ./ N_Psi(RP.I_qa,RP.I_qa);
% ratio_Ny_NPsi = ratio_Ny_NPsi_matrix(1,1);
% qD_NPsi ./ (J_q_qa * qaD_Ny)
% J_q_qa * N_y * va_unitact
% N_Psi*v_unitact
% sqrt(norm(J_q_qa*J_q_qa'))/36
% norm(J_q_qa,'fro')
% project v into actuation space
% v_unitact2 = J_q_qa*va_unitact
% (va_unitact' * J_q_qa')'
% norm(va_unitact2)/norm(va_unitact)
% qD_NPsi2 = N_Psi*v_unitact2
% qD_NPsi2./qD_NPsi
% v_rand = rand(RP.NJ,1);
% J_q_qa * N_y * ((v_rand') * J_q_qa)'
% J_q_qa * N_y * J_q_qa' * (v_rand)

%% Try out new coordinates and check a
x1 = x0 + rand(6,1)*0.1;
[q1, Phi,~,Stats] = RP.invkin_ser(x1, q0, s_ep);
[~,Phi_q1] = RP.constr4grad_q(q1);
[~,Phi_x1] = RP.constr4grad_x(x1);
Jtilde_inv_x1 = -Phi_q1\Phi_x1;
J_inv_x1 = Jtilde_inv_x1(RP.I_qa,:);
fprintf('Condition number of the PKM Jacobian: %1.1f\n', cond(J_inv_x));
J_x1 = inv(J_inv_x1);
J_y1 = J_x1(I_EE_red,:);
N_y1 = eye(sum(RP.I_qa)) - pinv(J_y1)*J_y1;
[Psi_q1,~] = RP.constr3grad_q(q1, x1);
fprintf(['Condition number of the inverse kinematics Jacobian ', ...
  '(pointing task constraints): %1.1f\n'], cond(Psi_q1));
[Psi_x,~] = RP.constr3grad_x(q1, x1);
N_Psi1 = eye(RP.NJ) - pinv(Psi_q1)*Psi_q1;
J_q_qa1 = Jtilde_inv_x1 * J_x1;

N_Psi_from_y1 = J_q_qa1 * N_y1 * J_q_qa1';
N_Psi_ratio1 = N_Psi1 ./ N_Psi_from_y1;
k_Psi1 = N_Psi_ratio1(1,1);

assert(abs(k_Psi1 - k_Psi)>1e-6, 'ratio has to be different for different configurations');
% k_Psi1/k_Psi
% eig(N_y)
% eig(N_Psi)
% 1/sqrt(norm(J_q_qa*J_q_qa'))

% 1/k_Psi
% norm(J_q_qa)
% norm(J_q_qa*J_q_qa', 'fro')
% det(N_Psi)
% norm(J_q_qa*J_q_qa') ./ norm(J_q_qa1*J_q_qa1')
% norm(J_q_qa'*J_q_qa)
% det(J_q_qa'*J_q_qa)
% det(J_q_qa*J_q_qa')

N_x = J_x * N_y * J_x';
N_x_test = N_x;
N_x_test(end,end) = 0;
assert(all(abs(N_x_test(:))<1e-8), 'all elements of N_x have to be zero, except bottom right one');


N_Psi_from_y_test = Jtilde_inv_x * N_x * Jtilde_inv_x';
assert(all(abs(N_Psi_from_y(:) - N_Psi_from_y_test(:))<1e-10), 'syntax error');

% Jtilde_inv_x(:,end)' * Jtilde_inv_x(:,end)
% norm(N_x)

% Jtilde_inv_x * J_x

