function dydt = Neuron_Standalone(t,y,Params,time,Iext)
%% Extract
Iin    = interp1(time,Iext,t);
Vmem   = y(1);
Vinv   = y(2);
Vspike = y(3);
V1     = y(4);

%% Compute Differential Part
dVmemdt = (Iin/(Params.K*Params.Ith)-Params.Gamma_1*(Vmem-V1))/Params.tau_m;
dVinvdt = (exp(-Params.K*Params.Beta*(Vmem + Vspike))*(1-exp(Vinv-Params.VDD))-exp(Params.K*Vmem)*(1-exp(-Vinv))/Params.Gamma_2)/Params.tau_p;
dVspikedt = (Params.IFGB/Params.Ith-((log((1+Params.E1*exp(Params.K*Vinv/2))./(1+Params.E1*exp((Params.K*Vinv-Vspike)/2)))).^2))/Params.tau_out;
dV1dt = ((Params.Gamma_1*Params.K/(Params.E1^2))*(Vmem-V1)-(1-exp(-V1))/exp(-Params.K*Vspike))/Params.tau_1;
dydt  = [dVmemdt;dVinvdt;dVspikedt;dV1dt];
end