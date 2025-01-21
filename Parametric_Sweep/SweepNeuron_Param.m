%% Sweep one of the parameters and find the Lyapunov exponent and the deviation of extrema points
%% Initialize
clear, clc;
addpath(genpath(pwd));

%% Define Sweep Type (1=sweep beta_s, 2=sweep beta_m, 3=sweep V_FG0)
SweepType = 1;

%% Define constants
% Biases
Params.UT   = 26.0e-3;
Params.VFG0 = 74.0234;    
Params.IFGB = 9.0820e-08; 

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
Params.Beta_m = 0.0654;
Params.Beta_s = 0.0654;

%% Derived constants
Params.C_z = 1/(1/Params.cin + 1/Params.cfb);
Params.C_alpha_square = (Params.cmem+Params.C_z)*(Params.cout+Params.C_z)-Params.C_z^2;

% Voltage-based scaling parameters
Params.Gamma_2 = exp(Params.K*(Params.VDD-Params.VFG0));
Params.E1 = exp(-Params.K*Params.VT0/2);

% Time constants
Params.tau_m   = Params.UT*Params.C_alpha_square/((Params.cout+Params.C_z)*Params.K*Params.Ith);
Params.tau_1   = Params.UT*Params.cp2*exp(Params.K*Params.VT0)/Params.Ith;
Params.tau_p   = Params.UT*Params.cp*exp(Params.K*(Params.VT0+Params.VFG0-Params.VDD))/Params.Ith;
Params.tau_out = Params.UT*Params.C_alpha_square/((Params.cmem+Params.C_z)*Params.Ith);

% Input parameters
if(SweepType==1)
    InputLevels = logspace(-3,-1,24*12).'; % Beta_s
elseif(SweepType==2)
    InputLevels = logspace(-3,-1,24*12).'; % Beta_m
elseif(SweepType==3)
    InputLevels = linspace(70,80,24*12).'; % V_FG0
end

%% Simulation Parameters
x0 = [-0.001 96 0.004 0].'; % Vmem, Vinv, Vspike
tf = 2e-3;
vtaus = [2.5,1.5,0.5]/Params.UT;
SpikeRate = zeros(length(InputLevels),length(vtaus));
Lyapunov = zeros(length(InputLevels),length(vtaus));
i = length(InputLevels);
j = length(vtaus);
Data(i,j).t = 0;
Data(i,j).y = 0;
PromThresh = 0.1;
FullOrbitData(length(InputLevels),length(vtaus),4).Min = [];
FullOrbitData(length(InputLevels),length(vtaus),4).Max = [];
Min = zeros([length(InputLevels),length(vtaus),4]);
Max = zeros([length(InputLevels),length(vtaus),4]);
PeriodData(length(InputLevels),length(vtaus)).Period = [];
Spread = zeros([length(InputLevels),length(vtaus)]);
for j = 1:length(vtaus)
    % Define Vtau-dependent params
    Params.vtau = vtaus(j);
    Params.Gamma_1 = (Params.vtau-Params.VT0)/2;

    parfor i = 1:length(InputLevels)
        % Prepare input signals
        InputCurrentLevel = 40e-9;
        SimTime = 0.01+3e-9/InputCurrentLevel;
        time = 0:1e-5:SimTime;
        Iin = zeros(size(time));
        Iin(time >= 1e-3) = InputCurrentLevel;
        tspan = [0, SimTime];

        %% Run transient simulation
        tstart = tic;
        if(SweepType==1)
            [t,y] = ode15s(@(t,y)Neuron_Betas(t,y,Params,time,Iin,InputLevels(i)),tspan,x0,odeset('Events',@(t,y) myevent(t,y,tstart)));
        elseif(SweepType==2)
            [t,y] = ode15s(@(t,y)Neuron_Betam(t,y,Params,time,Iin,InputLevels(i)),tspan,x0,odeset('Events',@(t,y) myevent(t,y,tstart)));
        elseif(SweepType==3)
            [t,y] = ode15s(@(t,y)Neuron_VFG0(t,y,Params,time,Iin,InputLevels(i)),tspan,x0,odeset('Events',@(t,y) myevent(t,y,tstart)));
        end
        Data(i,j).t = t;
        Data(i,j).y = y;
        try
            %% Lyapunov
            % Resample
            time = t(t > tf);
            myData = y(t > tf,1);
            [~,idx] = findpeaks(myData,"MinPeakProminence",0.1);
            SpikeRate(i,j) = 1./mean(diff(time(idx)),'omitnan');
            Spread(i,j) = std(mean(diff(time(idx)))./diff(time(idx)));
            PeriodData(i,j).Period = diff(time(idx));
            Fs = SpikeRate(i,j)*100 + 10;
            DataRes = interp1(time,myData,min(time)+(0:2000)/Fs);
            eRange = (1/SpikeRate(i,j))*Fs;

            % Calculate
            [~,eLag,eDim] = phaseSpaceReconstruction(DataRes,'MaxLag',100);
            Lyapunov(i,j) = lyapunovExponent(DataRes,Fs,eLag,eDim,'ExpansionRange',min([ceil(eRange/8),1e3]));
        catch
            Lyapunov(i,j) = NaN;
        end

        % Estimate min/max information
        for k = 1:4
            try
                Peaks = findpeaks(y(t > tf,k),'MinPeakProminence',PromThresh);
                Troughs = -findpeaks(-y(t > tf,k),'MinPeakProminence',PromThresh);
                FullOrbitData(i,j,k).Min = Troughs;
                FullOrbitData(i,j,k).Max = Peaks;
                if ((~isempty(Peaks)) || (~isempty(Troughs)))
                    Min(i,j,k) = mean(Troughs);
                    Max(i,j,k) = mean(Peaks);
                else
                    Min(i,j,k) = mean(y(t > tf,k));
                    Max(i,j,k) = mean(y(t > tf,k));
                    FullOrbitData(i,j,k).Min = Min(i,j,k);
                    FullOrbitData(i,j,k).Max = Max(i,j,k);
                end
            catch
                Min(i,j,k) = NaN;
                Max(i,j,k) = NaN;
            end
        end
    end
end

%% Plot Orbit Diagram
Labels = {'V_{mem} (V)','V_{inv} (V)','V_{spike} (V)','V_r (V)'};
Marker = {'k.','b.','m.'};
MinBoundary = zeros([2,length(InputLevels)]);
MaxBoundary = zeros([2,length(InputLevels)]);
for k = 1
    figure;
    for i = 1:length(InputLevels)
        for j = 1:length(vtaus)
            ydata = Params.UT*sort([FullOrbitData(i,j,k).Min-mean(FullOrbitData(i,j,k).Min);FullOrbitData(i,j,k).Max-mean(FullOrbitData(i,j,k).Max)]);
            if(~isempty(ydata))
                semilogx(InputLevels(i)*ones(size(ydata)),ydata,Marker{j},'MarkerSize',5);
                hold on;
            end
        end
    end

    legend('V_{\taur}=2.5V','V_{\taur}=1.5V','V_{\taur}=0.5V','location','best');
    if(SweepType==1)
        xlabel('\beta_s');
    elseif(SweepType==2)
        xlabel('\beta_m');
    elseif(SweepType==3)
        xlabel('V_{FG0}');
    end
    ylabel(Labels{k});
end

%% Make the Lyapunov Exponent Plot
W = 1;
figure;
semilogx(InputLevels,Lyapunov(:,1),'k');
hold on;
semilogx(InputLevels,Lyapunov(:,2),'b--');
semilogx(InputLevels,Lyapunov(:,3),'m.--');
ylabel('\lambda (s^{-1})')

figure;
semilogx(InputLevels,Lyapunov(:,1)./SpikeRate(:,1),'k');
hold on;
semilogx(InputLevels,Lyapunov(:,2)./SpikeRate(:,2),'b--');
semilogx(InputLevels,Lyapunov(:,3)./SpikeRate(:,3),'m.--');
ylabel('\lambda T_{spike}')
if(SweepType==1)
    xlabel('\beta_s');
elseif(SweepType==2)
    xlabel('\beta_m');
elseif(SweepType==3)
    xlabel('V_{FG0}');
end