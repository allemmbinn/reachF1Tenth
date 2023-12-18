function f = vehicle_kextend_7Eq(x,u)
% vehicle_kextend_7EqDT - system dynamics for the discrete-time version of the 
% Extended Kinematic Model
%
% Syntax:  
%    f = vehicle_kextend_7EqDT(x,u,h)
%
% Inputs:
%    x - state vector
%    u - input vector
%    h - time step size
%
% Outputs:
%    f - new state vector
% 
% References:
%    [1] M. Althoff et al. "Reachability analysis of nonlinear systems with 
%        uncertain parameters using conservative linearization", CDC 2008

% Author:        Matthias Althoff
% Written:       11-Dec-2023
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % parameter
    lr = 0.17145;
    lf = 0.15875;

    % differential equations
    f(1,1) = (x(4)*cos(x(3)) - x(5)*sin(x(3)));           % x
    f(2,1) = (x(4)*sin(x(3)) + x(5)*cos(x(3)));           % y
    f(3,1) = (x(6));                                      % phi
    f(4,1) = (u(1));                                      % vx
    f(7,1) = (u(2));                                      % delta
    f(5,1) = (lr/(lf+lr) * (f(7,1)*x(4) + x(7)*f(4,1)));  % vy
    f(6,1) = (1/(lf+lr) * (f(7,1)*x(4) + x(7)*f(4,1)));   % w
                                      

%------------- END OF CODE --------------