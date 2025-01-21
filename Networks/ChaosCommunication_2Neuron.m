%% This Script Simulates the response of a Reservoir of Chaotic Neurons
% Written by Swagat Bhattacharyya 11/29/24
%% Initialize
clear, clc;
addpath(genpath(pwd));
load('NeuronParameters.mat');

%% Define Ensemble Connectivity Via the Synapse Weight Matrix 
Params.NeuronPopulation = 2;
SynapseBias = 450e-9;
Params.ReservoirWeight_Mat = [0 SynapseBias;SynapseBias 0];

%% Compute Input Current By linearly transforming the data
fs = 1e6; % Sampling frequency in Hz
SimTime = 50e-3;
t = 0:1/fs:SimTime-1/fs; % Time vector
Params.time = 0:1/fs:SimTime;
Params.Input = zeros(Params.NeuronPopulation,size(Params.time,2));
Params.Input(1,(Params.time >= 1e-3) & (Params.time < 2e-3)) = 30e-9;
tspan = [0, SimTime];

%% Define Simulation Params and Run
% Guess initial conditions
y0 = Params.VDD*ones(Params.NeuronPopulation,1);
for i = 1:Params.NeuronPopulation
    y0 = [y0;-0.001; 96; 0.004; 0]; % Vmem, Vinv, Vspike, Vr
end

% Run the actual sim
[tvec,y] = ode15s(@(t,y)NetworkODE(t,y,Params),t,y0);

%% Do the chaos communication
% Extract mask and unmask
idx = t > 3e-3;
t_clip = t(idx);
Key1 = 0.5e-3*y(idx,(Params.NeuronPopulation+1)).*y(idx,(Params.NeuronPopulation+2)).*y(idx,(Params.NeuronPopulation+3));
Key2 = 0.5e-3*y(idx,(Params.NeuronPopulation+5)).*y(idx,(Params.NeuronPopulation+6)).*y(idx,(Params.NeuronPopulation+7));

% Original signal: sum of two sinusoids
signal = 0.75*sin(2 * pi * 1000 * t) + 0.5*sin(2 * pi * 2000 * t) + sin(2 * pi * 3000 * t);
signal = 0.5e-3*signal(idx).';

% FFT of the original signal
fft_original = fft(signal);
frequencies = (0:length(signal)-1) * fs / length(signal);

% Masking using Key1
masked_signal = signal + Key1;

% FFT of the masked signal
fft_masked = fft(masked_signal);

% Unmasking using Key2
unmasked_signal = masked_signal - Key2;

% FFT of the unmasked signal
fft_unmasked = fft(unmasked_signal);

% Plotting
figure;
Xrange = [0.5,3.5];
% Original signal FFT
stop = floor(length(fft_original)/2);
subplot(3, 1, 1);
set(gcf,'defaultlinelinewidth',1)
plot(frequencies(1:stop)/1000, abs(fft_original(1:stop)),'k');
ylabel('|X_{in}(f)|')
xlim(Xrange);
ylim([0,12]);

% Masked signal FFT
subplot(3, 1, 2);
plot(frequencies(1:stop)/1000, abs(fft_masked(1:stop)),'k');
ylabel('|X_{send}(f)|')
xlim(Xrange);
ylim([0,22]);

% Unmasked signal FFT
subplot(3, 1, 3);
plot(frequencies(1:stop)/1000, abs(fft_unmasked(1:stop)),'k');
ylabel('|X_{out}(f)|')
xlabel('Frequency (kHz)')
xlim(Xrange);
ylim([0,12]);