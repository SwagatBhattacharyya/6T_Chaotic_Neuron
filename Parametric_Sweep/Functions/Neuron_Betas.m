function dydt = Neuron_Betas(t,y,Params,time,Iext,Beta_s)
%% Extract
Iin    = interp1(time,Iext,t);
Vmem   = y(1);
Vinv   = y(2);
Vspike = y(3);
Vr     = y(4);

%% Derived Constants
Im1_0 = Params.Ith*exp(-Params.K*Params.VT0);
Im3_0 = Params.Ith*exp(-Params.K*Params.VT0);
Im4_0 = Params.Ith*exp((Params.K*(Params.VDD-Params.VT0)));

%% Compute Algebraic Part
V2 = Params.VFG0 + Params.Beta_m*Vmem + Beta_s*Vspike;
E1 = exp((Params.K*(Vinv-Params.VT0))/2);

%% Compute Transistor Currents
Im1 = Im1_0*exp(Params.K*Vspike)*(1-exp(-Vr));
Im2 = Params.Ith*Params.K*Params.Gamma_1*(Vmem-Vr);
Im3 = Im3_0*exp(Params.K*Vmem)*(1-exp(-Vinv));
Im4 = Im4_0*exp(-Params.K*V2)*(1-exp(Vinv-Params.VDD));
Im5 = Params.Ith*((log(1+E1)).^2 - (log(1+E1*exp(-Vspike/2))).^2);
Im6 = Params.IFGB;

%% Compute Differential Part
dVmemdt   = (Iin-Im2)*(Params.cout+Params.C_z)/Params.C_alpha_square;
dVinvdt   = (Im4-Im3)/Params.cp;
dVspikedt = (Im6-Im5)*(Params.cmem+Params.C_z)/Params.C_alpha_square;
dV1dt     = (Im2-Im1)/Params.cp2;
dydt      = [dVmemdt;dVinvdt;dVspikedt;dV1dt]/Params.UT;
end