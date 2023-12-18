% function completed = example_nonlinear_reach_f1Tenth()

% example_nonlinear_reach_f1Tenth() - example of 
%     nonlinear reachability analysis
%
% Syntax:
%    completed = example_nonlinear_reach_05_autonomousCar()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false.

% Authors:       Allen Emmanuel Binny
% Written:       10-December-2023
% Last update:   
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------
params.tFinal = 3.00;

% intial set
params.R0 = zonotope([[0; 0; 0; 0; 0 ; 0; 0],...
                      diag([0.3, 0.3, 0.001, 0.1, 0.1, 0.001, 0.001])]);

% uncertain inputs
params.U = zonotope([[0; 0],diag([0.1;0.001])]);

% Reachability Settings ---------------------------------------------------

options.timeStep = 0.01;
options.taylorTerms = 3;
options.zonotopeOrder = 100;

options.alg = 'lin';
options.tensorOrder = 2;


% System Dynamics ---------------------------------------------------------
vehicle = nonlinearSys(@vehicle_kextend_7Eq,7,2);

% Specifications ----------------------------------------------------------
% original specs: |vx| <= 7, |delta| <= 0.42 should be fulfilled at all times
vlim = 7;
deltalim = 0.41;
alim = 7.5;
ratlim = deltalim/vlim;

hs1 = halfspace([0 0 0 -1 0 0 0],-vlim);
hs2 = halfspace([0 0 0 1 0 0 0],-vlim);
hs3 = halfspace([0 0 0 -ratlim 0 0 -1],-deltalim);
hs4 = halfspace([0 0 0 ratlim 0 0 1],-deltalim);
spec = specification({hs1,hs2,hs3,hs4},'unsafeSet');

% Reachability Analysis ---------------------------------------------------
% Input Dictionary --------------------------------------------------------
n = 2; % Discritisation of the input
uinp = zeros(2,n*n);
for i=1:n
    for j=1:n
        acc = -alim + 2*alim/(n+1) * (i);          
        del = -deltalim + 2*deltalim/(n+1) * (j); 
        uinp(:,(i-1)*n+j) = [acc;del];
    end
end
% Defining the Control Inputs
tComp = 0;
maxU = [7;0];
inpU = [repmat(maxU,1,100),repmat(maxU,1,100),repmat(maxU,1,100)];
params.u = inpU;
tic
R(n^6+1) = reach(vehicle, params, options);  
% tComp = tComp + toc; 
for i=1:n*n
    for j=1:n*n
        for k=1:n*n
            inpU = [repmat(uinp(:,i),1,100),repmat(uinp(:,j),1,100),repmat(uinp(:,k),1,100)];
            params.u = inpU;
            % tic
            [R((i-1)*n*n + (j-1) * n + k),res] = reach(vehicle, params, options,spec);  
            % tComp = toc + tComp
        end
    end
end
tComp = toc + tComp;
disp(['computation time of reachable set: ',num2str(tComp)]);

% % Simulation --------------------------------------------------------------
% % simulation settings
% simOpt.points = 10;
% % random simulation
% simRes = simulateRandom(vehicle, params, simOpt);
% Visualization -----------------------------------------------------------
figure;  
hold on; box on
projDim = [1 2];

% plot reachable sets
useCORAcolors("CORA:contDynamics")
plot(R,projDim,'DisplayName','Reachable set');
    
% plot initial set
plot(R(1).R0,projDim, 'DisplayName','Initial set');
    
% % plot simulation results      
% plot(simRes,projDim,'DisplayName','Simulations');    

% label plot    
xlabel('x');    
ylabel('y');    
legend()

% example completed
completed = true;
% ------------------------------ END OF CODE ------------------------------