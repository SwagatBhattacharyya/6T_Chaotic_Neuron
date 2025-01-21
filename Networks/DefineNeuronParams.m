%% Initialize
clear, clc;

%% Define constants
% Biases
Params.UT   = 26.0e-3;
Params.VFG0 = 74.0234;      
Params.IFGB = 9.0820e-08;  
Params.vtau = 2.5/Params.UT;

% Capacitances
Params.cp2  = 1e-15;
Params.cin  = 10e-15;
Params.cfb  = 10e-15;

Params.cp   = 7.9122e-11;
Params.cmem = 5.7357e-11;
Params.cout = 2.1600e-11;

% Transistor Params
Params.VDD  = 96;
Params.VT0  = 18.8887;
Params.K    = 0.6787;
Params.Ith  = 1.6572e-06;
Params.Beta = 0.0654;

%% Derived constants
Params.C_z = 1/(1/Params.cin + 1/Params.cfb);
Params.C_alpha_square = (Params.cmem+Params.C_z)*(Params.cout+Params.C_z)-Params.C_z^2;

% Voltage-based scaling parameters
Params.Gamma_2 = exp(Params.K*(Params.VDD-Params.VFG0));
Params.E1 = exp(-Params.K*Params.VT0/2);
Params.Gamma_1 = (Params.vtau-Params.VT0)/2;

% Time constants
Params.tau_m   = Params.UT*Params.C_alpha_square/((Params.cout+Params.C_z)*Params.K*Params.Ith);
Params.tau_p   = Params.UT*Params.cp*exp(Params.K*(Params.VT0+Params.VFG0-Params.VDD))/Params.Ith;
Params.tau_out = Params.UT*Params.C_alpha_square/((Params.cmem+Params.C_z)*Params.Ith);
Params.tau_1   = Params.UT*Params.cp2*exp(Params.K*Params.VT0)/Params.Ith;

% Triangle generator
Params.M = 1000;         % Big M for triangle gen
Params.SpikeThresh = 20; % In UT, based off orbit diagram
Params.tau_tr = 10e-6; 
Params.tau_tf = Params.tau_tr/5; % The fall time constant should be ~5x faster.

%% Save
save('NeuronParameters.mat')