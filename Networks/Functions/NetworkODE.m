function dydt = NetworkODE(t,y,Params)
% y order: Triangle Inputs, Neuron State
%% Some More Parameters...
EffCoupling = 0.03; % Combination of both the CG-gate and gate-channel coupling

%% Extract and format stuff
dydt = zeros(length(y),1);                   % Preallocation
Vmem = y((Params.NeuronPopulation+1):4:end); % The membrane voltage of each neuron
Vtri = y(1:Params.NeuronPopulation);

%% Compute input current
InputCurrent = Params.ReservoirWeight_Mat*exp(-EffCoupling*Vtri); % Exponent goes between 0.05 and 1;
for i = 1:length(InputCurrent)
    InputCurrent(i) = InputCurrent(i) + interp1(Params.time,Params.Input(i,:),t);
end

%% Solve for triangle generator outputs
for i = 1:Params.NeuronPopulation
    dydt(i) = Triangle(t,y(i),Params,Vmem(i));     % Solve for dydt of the triangle generators
end

%% Solve for the neuron update
for i = 1:Params.NeuronPopulation
    StartIdx = (Params.NeuronPopulation+1)+4*(i-1); % Corresponds to the index of the first state variable of the current neuron
    dydt(StartIdx:(StartIdx+3)) = Neuron(t,y(StartIdx:(StartIdx+3)),Params,InputCurrent(i)); % Solve for dydt of the current neuron
end
end