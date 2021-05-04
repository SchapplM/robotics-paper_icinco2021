% Computation of the root locus of an I2-system with P controller.
% This proves the statement in Sec. 5 of the paper that this kind of linear
% control system is not asymptotically stable.

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2021-05
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear

% Define the transfer function
Fs = tf(1, [1 0 0]);

% Get the roots of the closed loop system, which assumes control by a P
% controller
R = rlocus(Fs, [0.1,1,10]);
assert(all(real(R(:))==0), 'All poles of the system have to be on the imag axis');