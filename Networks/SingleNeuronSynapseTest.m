%% This Script Simulates the response of a Single Neuron and Synapse
% Written by Swagat Bhattacharyya 11/29/24
%% Initialize
clear, clc;
addpath(genpath(pwd));
load('NeuronParameters.mat');

%% Define Network Connectivity Via the Synapse Weight Matrix 
Params.NeuronPopulation = 1;           % Just one neuron
Params.ReservoirWeight_Mat = [0]*1e-9; % No connectivity here...

%% Estimate y0
y0 = Params.VDD*ones(Params.NeuronPopulation,1);
for i = 1:Params.NeuronPopulation
    y0 = [y0;-0.001; 96; 0.004; 0]; % Vmem, Vinv, Vspike, Vr
end

%% Define Simulation Params and Run
% Prepare input signals
InputCurrentLevel = 30e-9;
SimTime = 0.01+3e-9/InputCurrentLevel;
Params.time = 0:1e-5:SimTime;
Params.Input = zeros(size(Params.time));
Params.Input(Params.time >= 1e-3) = InputCurrentLevel;
tspan = [0, SimTime];

% Run the actual sim
[t,y] = ode15s(@(t,y)NetworkODE(t,y,Params),tspan,y0);

%% Plot
figure;
subplot(2,1,1), plot(t*1e3,y(:,1));
ylabel('V_{tri} (U_T)')
subplot(2,1,2), plot(t*1e3,y(:,2));
ylabel('V_{mem} (U_T)')
xlabel('Time (ms)')
ylabel('Voltage (U_T)')