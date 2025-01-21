%% Initialize
clear, clc;
addpath(genpath(pwd));
load('NeuronParameters.mat');

%% Define Input
time = 0:1e-5:10e-3;
Vmem_ext = zeros(size(time));
Vmem_ext((time>=1e-3) & (time<2e-3)) = 21;
tspan = [0,max(time)]; 
x0 = 0;

%% Simulate
[t,y] = ode15s(@(t,y)Triangle_Standalone(t,y,Params,time,Vmem_ext),time,x0);

%% Plot
figure, hold on;
plot(time*1e3,Vmem_ext);
plot(t*1e3,y);
legend('V_{mem}','V_{tri}')
ylabel('Voltage (U_T)')
xlabel('Time (ms)')