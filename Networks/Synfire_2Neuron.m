%% This Script Simulates the response of a Reservoir of Chaotic Neurons
% Written by Swagat Bhattacharyya 11/29/24
%% Initialize
clear, clc;
addpath(genpath(pwd));
load('NeuronParameters.mat');

%% Define Ensemble Connectivity Via the Synapse Weight Matrix 
SynapseBias = 40e-9;
Params.NeuronPopulation = 2;
Params.ReservoirWeight_Mat = [0 SynapseBias;SynapseBias 0];

%% Compute Input Current By linearly transforming the data
SimTime = 50e-3;
Params.time = 0:1e-5:SimTime;
Params.Input = zeros(Params.NeuronPopulation,size(Params.time,2));
Params.Input(1,(Params.time >= 10e-3) & (Params.time < 11e-3)) = 30e-9;
tspan = [0, SimTime];

%% Define Simulation Params and Run
% Guess initial conditions
y0 = Params.VDD*ones(Params.NeuronPopulation,1);
for i = 1:Params.NeuronPopulation
    y0 = [y0;-0.001; 96; 0.004; 0]; % Vmem, Vinv, Vspike, Vr
end

% Run the actual sim
[t,y] = ode15s(@(t,y)NetworkODE(t,y,Params),tspan,y0);

%% Plot
Xmax = 10;
XOffset = 9;
figure;
subplot(5,1,1:2), hold on;
plot(t*1e3-XOffset,y(:,Params.NeuronPopulation+1),'k');
plot(t*1e3-XOffset,y(:,Params.NeuronPopulation+5),'b--');
ylabel('V_{mem} (U_T)')
xlim([0,Xmax]);
subplot(5,1,3:4), hold on;
plot(t*1e3-XOffset,y(:,1),'k');
plot(t*1e3-XOffset,y(:,2),'b--');
ylabel('V_{tri} (U_T)')
xlim([0,Xmax]);
subplot(5,1,5), plot(Params.time*1e3-XOffset,Params.Input*1e9,'k');
xlim([0,Xmax]);
xlabel('Time (ms)')
ylabel('I_{in} (nA)')
set(gcf, 'Units','centimeters', 'Position',[8 0 14 9])

%% Plot V_mem2 vs V_mem1
idx = t > 12e-3;
t_clip = t(idx);
Vmem1 = y(idx,(Params.NeuronPopulation+1));
Vmem2 = y(idx,(Params.NeuronPopulation+5));
figure;
plot(Vmem1,Vmem2);
xlabel('V_{mem,1}');
ylabel('V_{mem,2}');
axis tight