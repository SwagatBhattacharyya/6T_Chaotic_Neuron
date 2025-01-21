%% Initialize
clear, clc;
addpath(genpath(pwd));
load('NeuronParameters.mat');

%% Simulation Parameters
x0 = [-0.001 96 0.004 0].'; % Vmem, Vinv, Vspike, Vr

% Prepare input signals
InputCurrentLevel = 30e-9;
SimTime = 0.01+3e-9/InputCurrentLevel;
time = 0:1e-6:SimTime;
Iin = zeros(size(time));
Iin(time >= 1e-3) = InputCurrentLevel;
tspan = [0, SimTime];
[t,y] = ode15s(@(t,y)Neuron_Standalone(t,y,Params,time,Iin),tspan,x0);

%% Plot
figure;
subplot(5,1,1);
plot(t*1e3,y(:,1));
ylabel('V_{mem} (U_T)')
subplot(5,1,2);
plot(t*1e3,y(:,2));
ylabel('V_{inv} (U_T)')
subplot(5,1,3);
plot(t*1e3,y(:,3));
ylabel('V_{spike} (U_T)')
subplot(5,1,4);
plot(t*1e3,y(:,4));
ylabel('V_{r} (U_T)')
subplot(5,1,5);
plot(time*1e3,Iin*1e9);
ylabel('I_{in} (nA)')
xlabel('Time (ms)')

figure;
plot3(y(:,1),y(:,2),y(:,3));
xlabel('V_{mem} (U_T)')
ylabel('V_{inv} (U_T)')
zlabel('V_{spike} (U_T)')