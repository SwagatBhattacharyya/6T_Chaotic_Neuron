%% Sweep the input current and find the Lyapunov exponent and the orbit diagram
%% Initialize
clear, clc;
addpath(genpath(pwd));

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
Params.Beta = 0.0654;

%% Derived constants
Params.C_z = 1/(1/Params.cin + 1/Params.cfb);
Params.C_alpha_square = (Params.cmem+Params.C_z)*(Params.cout+Params.C_z)-Params.C_z^2;

% Voltage-based scaling parameters
Params.Gamma_2 = exp(Params.K*(Params.VDD-Params.VFG0));
Params.E1 = exp(-Params.K*Params.VT0/2);

% Time constants
Params.tau_m   = Params.UT*Params.C_alpha_square/((Params.cout+Params.C_z)*Params.K*Params.Ith);
Params.tau_1   = Params.UT*Params.cp2*exp(Params.K*Params.VT0)/Params.Ith;
Params.tau_p   = Params.UT*Params.cp*exp(Params.K*Params.VT0)/Params.Ith;
Params.tau_out = Params.UT*Params.C_alpha_square/((Params.cmem+Params.C_z)*Params.Ith);

%% Simulation Parameters
x0 = [-0.001 96 0.004 0].'; % Vmem, Vinv, Vspike
tf = 2e-3;
InputCurrentLevels = logspace(-12,-4,24*8).';
vtaus = [2.5,1.5,0.5]/Params.UT;
SpikeRate = zeros(length(InputCurrentLevels),length(vtaus));
Lyapunov = zeros(length(InputCurrentLevels),length(vtaus));
i = length(InputCurrentLevels);
j = length(vtaus);
Data(i,j).t = 0;
Data(i,j).y = 0;
PromThresh = 0.1;
FullOrbitData(length(InputCurrentLevels),length(vtaus),4).Min = [];
FullOrbitData(length(InputCurrentLevels),length(vtaus),4).Max = [];
Min = zeros([length(InputCurrentLevels),length(vtaus),4]);
Max = zeros([length(InputCurrentLevels),length(vtaus),4]);
for j = 1:length(vtaus)
    % Define Vtau-dependent params
    Params.vtau = vtaus(j);
    Params.Gamma_1 = (Params.vtau-Params.VT0)/2;

    parfor i = 1:length(InputCurrentLevels)
        % Prepare input signals
        InputCurrentLevel = InputCurrentLevels(i);
        SimTime = 0.01+3e-9/InputCurrentLevel;
        time = 0:1e-5:SimTime;
        Iin = zeros(size(time));
        Iin(time >= 1e-3) = InputCurrentLevel;
        tspan = [0, SimTime];

        %% Run transient simulation
        tstart = tic;
        [t,y] = ode15s(@(t,y)Neuron_IF(t,y,Params,time,Iin),tspan,x0,odeset('Events',@(t,y) myevent(t,y,tstart)));

        Data(i,j).t = t;
        Data(i,j).y = y;
        try
            %% Lyapunov
            % Resample
            time = t(t > tf);
            myData = y(t > tf,1);
            [~,idx] = findpeaks(myData,"MinPeakProminence",0.1);
            SpikeRate(i,j) = 1./mean(diff(time(idx)),'omitnan');
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
            end
            catch
                Min(i,j,k) = NaN;
                Max(i,j,k) = NaN;
            end
        end
    end
end

%% Plot IF Curve
figure;
SpikeRate(SpikeRate == 0) = NaN;
loglog(InputCurrentLevels,SpikeRate(:,1),'k')
hold on;
loglog(InputCurrentLevels,SpikeRate(:,2),'b--')
loglog(InputCurrentLevels,SpikeRate(:,3),'m.--')
xlabel('I_{in} (A)')
ylabel('Mean Spike Rate (Hz)')

figure;
semilogx(InputCurrentLevels,Lyapunov(:,1),'k');
hold on;
semilogx(InputCurrentLevels,Lyapunov(:,2),'b--');
semilogx(InputCurrentLevels,Lyapunov(:,3),'m.--');
xlabel('I_{in} (A)')
ylabel('\lambda (s^{-1})')

figure;
semilogx(InputCurrentLevels,Lyapunov(:,1)./SpikeRate(:,1),'k');
hold on;
semilogx(InputCurrentLevels,Lyapunov(:,2)./SpikeRate(:,2),'b--');
semilogx(InputCurrentLevels,Lyapunov(:,3)./SpikeRate(:,3),'m.--');
xlabel('I_{in} (A)')
ylabel('\lambda T_{spike}') 

%% Plot Orbit Diagram
Labels = {'V_{mem} (V)','V_{inv} (V)','V_{spike} (V)','V_r (V)'};
Marker = {'k.','b.','m.'};
W = 5;
P = 4000;
MinBoundary = zeros([2,length(InputCurrentLevels)]);
MaxBoundary = zeros([2,length(InputCurrentLevels)]);
for k = 1:4
    figure;
    semilogx(InputCurrentLevels,Params.UT*medfilt1(Min(:,1,k),W),'k--');
    hold on;
    semilogx(InputCurrentLevels,Params.UT*medfilt1(Min(:,2,k),W),'b--');
    semilogx(InputCurrentLevels,Params.UT*medfilt1(Min(:,3,k),W),'m--');
    semilogx(InputCurrentLevels,Params.UT*medfilt1(Max(:,1,k),W),'k');
    semilogx(InputCurrentLevels,Params.UT*medfilt1(Max(:,2,k),W),'b');
    semilogx(InputCurrentLevels,Params.UT*medfilt1(Max(:,3,k),W),'m');
    for i = 1:length(InputCurrentLevels)
        for j = 1:length(vtaus)
            ydata = Params.UT*sort([FullOrbitData(i,j,k).Min;FullOrbitData(i,j,k).Max]);
            if(~isempty(ydata))
                Midline = mean(ydata);
                Top = ydata(ydata >= Midline);
                Bot = ydata(ydata < Midline);
                
                Msize = P*std(Top)+0.1;
                if isempty(Msize) || isnan(Msize)
                    Msize = 0.1;
                end
                semilogx(InputCurrentLevels(i)*ones(size(Top)),Top,Marker{j},'MarkerSize',Msize);
                hold on;

                Msize = P*std(Bot)+0.1;
                if isempty(Msize) || isnan(Msize)
                    Msize = 0.1;
                end
                semilogx(InputCurrentLevels(i)*ones(size(Bot)),Bot,Marker{j},'MarkerSize',Msize);
                hold on;
            end
        end
    end

    legend('V_{\taur}=2.5V','V_{\taur}=1.5V','V_{\taur}=0.5V','location','best');
    xlabel('I_{in} (A)');
    ylabel(Labels{k});
end