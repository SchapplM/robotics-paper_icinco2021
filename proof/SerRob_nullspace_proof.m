% Compute nullspace projector matrices for a six-axis industrial robot.
% This proves statements in Sec. 3 of the paper.
% 
% Dependencies: Serial Robot Database, Robotics Toolbox; see README.MD
% 
% Reference:
% [SchapplerTapOrt2019a] Resolution of Functional Redundancy for 3T2R Robot
% Tasks using Two Sets of Reciprocal Euler Angles, Proc. IFToMM W.C., 2019.

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2021-05
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clear
clc

%% Initialization of the Robot Model
if isempty(which('serroblib_path_init.m'))
  error('The serial robot database is not initialized in Matlab.');
end
SName='S6RRRRRR10V2'; % RRRRRR Industrial Robot
RName='S6RRRRRR10V2_KUKA1'; % KUKA KR 30
RS = serroblib_create_robot_class(SName, RName);
RS.fill_fcn_handles(false, false);
% Set Task coordinates y to 3T2R (see Sec. 2 of the paper)
RS.I_EE_Task = logical([1 1 1 1 1 0]);
% define end effector simulating a tool (e.g. milling spindle).
% The transformation has to be set, otherwise the nullspace is only the
% last axis.
RS.update_EE([-0.2, 0, -0.2]', [0,-pi/2,0]');
% Define initial pose and random component. This generalizes the results
% and should avoid special cases.
q0 = [0;90+30;-60;0;-30;0]*pi/180 + rand(RS.NJ,1)*0.3;

%% Plot Robot
s_plot = struct( 'ks', RS.NJ+2, 'straight', 0);
figure(1);clf;
hold on; grid on; view(3);
xlabel('x in m'); ylabel('y in m'); zlabel('z in m');
RS.plot( q0, s_plot );

%% Derivation of the Nullspace. Proof statements from Sec. 3.1
% Calculate forward kinematics. See equ. 1 (\ref{eq:ser_fkin})
T_E = RS.fkineEE(q0);
xE = RS.t2x(T_E);
% Calculate analytic Jacobian of the manipulator
J_x = RS.jacobia(q0); % See equ. 6 (\refl{eq:ser_diff_kin})
% Define selection matrix of Task coordinates and define task Jacobian
Py = eye(5,6); % see text after equ. 6 (\ref{eq:ser_diff_kin})
J_y = Py * J_x;
% Nullspace projector of equ. 5 (\ref{eq:ser_position_ik}).
% Calculate Psi_dq for Psi=0.
Phi_q = RS.constr2grad(q0, T_E(1:3,:), true);
% Select all entries except number 4. This corresponds to the redundant
% task space coordinate (reciprocal Euler angles, [SchapplerTapOrt2019a])
Psi_q = Phi_q([1:3 5 6],:);
N_Psi = eye(RS.NJ) - pinv(Psi_q)*Psi_q;
% Nullspace projector of equ. 7 (\ref{eq:ser_accel_ik})
N = eye(RS.NJ) - pinv(J_y)*J_y;
% Identity of the nullspace projectors
% Statement of equ. 11 (\label{eq:ser_nullspace_projector})
test_N = N-N_Psi;
assert(all(abs(test_N(:))<1e-10), ['nullspace projectors from Psi_dq and ', ...
  'J_y have to be identical']);

%% Calculation of Nullspace Gradients. Proof statements from Sec. 3.2
% General case from Sec. 3.2.1.
% Take condition number as performance criterion h.
% See equ. 8 (\ref{eq:ser_diffquot_general})
h_dq = NaN(1,RS.NJ);
for i = 1:RS.NJ
  dqi = 1e-8;
  dq = zeros(RS.NJ,1);
  dq(i) = dqi;
  J_x_i = RS.jacobia(q0+dq);
  h_dq(i) = (cond(J_x_i) - cond(J_x))/dqi;
end

% Calculate gradient in x coordinates instead of q coordinates.
% Corresponds to equ. 9 (\ref{eq:ser_gradient_opspace}); second part.
h_dx = NaN(1,6); % dh/dx; first term of second part of equation.
for i = 1:6
  % difference quotient with 1e-8 is a compromise for linearization error
  % and numerical precision of the calculation.
  dxi = 1e-8;
  dx = zeros(6,1);
  dx(i) = dxi;
  dq = J_x\dx;
  J_x_i = RS.jacobia(q0+dq);
  h_dx(i) = (cond(J_x_i) - cond(J_x))/dxi;
end
% Use second part of the equation. J_x is dx/dq.
h_dq_fromx = h_dx * J_x;
% Gradient from all task space coordinates is equal to joint space gradient
abserr_hdq = h_dq_fromx - h_dq;
relerr_hdq = abserr_hdq./h_dq;
% This supports equ. 9
assert(all(abs(relerr_hdq)<1e-3) && all(abs(abserr_hdq)<1e-4), ...
  'gradient vector from x or q coordinates has to be identical');

% Nullspace projection of the complete gradient dh/dx.
v_fromdq = N*h_dq(:);
v_fromdx = N*h_dq_fromx(:);
abserr_v = v_fromdq - v_fromdx;
relerr_v = abserr_v ./ v_fromdq;
% syntax check. Since h_dq is identical, the projection has to be identical
% as well
assert(all(abs(relerr_v)<1e-4) && all(abs(abserr_v)<1e-3), ['nullspace ', ...
  'vector from full x or q coordinates has to be identical']);

% Decomposition in single components of x. This corresponds to the single
% terms in equ. 12 (\ref{eq:serrob_nullspace_elim})
v2_allx = NaN(RS.NJ,6); % column-wise nullspace vectors v
h_dq_fromdx_allx = NaN(6,RS.NJ);
for i = 1:6
  % One of the summands in equ. 10 (ref{eq:ser_gradient_algebraic_conversion})
  h_dq_i = h_dx(i)*J_x(i,:);
  h_dq_fromdx_allx(i,:) = h_dq_i; % store for later evaluation
  % one of the expressions in equ. 12
  v2_allx(:,i) = N*h_dq_i(:);
end
% This proofs equ. 12
assert(all(all(abs(v2_allx(:,1:5))<1e-10)), ['gradients w.r.t other ', ...
  'coordinates than phi_z have to be zero after nullspace projection']);
% Re-check for plausibility, if the last gradient is identical to the
% previously calculated. This would find programming errors in this script.
abserr_decomp_x = sum(h_dq_fromdx_allx) - h_dq_fromx;
assert(all(abs(abserr_decomp_x)<1e-10), ['decomposition of the gradient ', ...
  'w.r.t x did not work (syntax error?)']);

% Using the phiz-gradient for projection is equivalent.
% This prooves eq. 14 (\ref{eq:serrob_nullspaceproj_equivalence})
abserr_v_phiz = v2_allx(:,6) - v_fromdq;
relerr_v_phiz = abserr_v_phiz ./ v_fromdq;
assert(all(abs(relerr_v_phiz)<1e-3) && all(abs(abserr_v_phiz)<1e-3), ...
  'nullspace vector from phiz or q coordinates has to be identical');

% Both h_dq are not equal, if only one task space coordinate is changed.
% This shows the second part of equ. 14.
test_h_from_q_vs_phiz = h_dq - h_dq_fromdx_allx(6,:);
assert(any(abs(test_h_from_q_vs_phiz)>1e-2), ...
  'gradient vectors regarding q and phiz have to be different');

fprintf('The script run through. All checks passed. Sec. 3 of paper is valid.\n');