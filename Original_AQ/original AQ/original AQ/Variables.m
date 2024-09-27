% Declaration of the common blocks & Setting problem variables
% The units are in centimeters

variables.pdeg = 1; % polynomial order
variables.type = 3; % type of stenosis

variables.N = 100; % number of elements
variables.L = 6; % length of vessel
xanchor = 0; % left ends
variables.cardiogrid = linspace(xanchor, variables.L, variables.N+1);

% Physical parameters
variables.alpha = 4/3;
variables.nu = 0.04;
E = 2e4; % Young's module
h = 0.06; % thickness of the vessel
variables.sigma = 0.5; % poission ratio
variables.rho = 1.055; % density
G0 = (E*h) / (0.18^2) / (1 - variables.sigma^2);
variables.G0_ = (E*h) / (1 - variables.sigma^2);
FinalTime = 1.0;
% dt = (cardiogrid(2) - cardiogrid(1)) / sqrt(G0/rho/2);
variables.dt = 1e-5;
variables.M = 30; % coefficient in limiter

% Gaussian quadrature abscissae and weights
variables.xg(1) = -0.981560634246719;
variables.xg(2) = -0.904117256370475;
variables.xg(3) = -0.769902674194305;
variables.xg(4) = -0.587317954286617;
variables.xg(5) = -0.367831498998180;
variables.xg(6) = -0.125233408511469;
variables.xg(7) =  0.125233408511469;
variables.xg(8) =  0.367831498998180;
variables.xg(9) =  0.587317954286617;
variables.xg(10) =  0.769902674194305;
variables.xg(11)=  0.904117256370475;
variables.xg(12)=  0.981560634246719;

variables.wg(1) =  0.047175336386512;
variables.wg(2) =  0.106939325995318;
variables.wg(3) =  0.160078328543346;
variables.wg(4) =  0.203167426723066;
variables.wg(5) =  0.233492536538355;
variables.wg(6) =  0.249147045813403;
variables.wg(7) =  0.249147045813403;
variables.wg(8) =  0.233492536538355;
variables.wg(9) =  0.203167426723066;
variables.wg(10)=  0.160078328543346;
variables.wg(11)=  0.106939325995318;
variables.wg(12)=  0.047175336386512;
