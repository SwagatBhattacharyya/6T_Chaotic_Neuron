function dydt = Triangle_Standalone(t,y,Params,time,Vmem_ext)
% Triangle Generator ODE (for generating Rall Alpha)
Vmem    = interp1(time,Vmem_ext,t);
dydt = sigmoid(Vmem-Params.SpikeThresh,Params.M)*sigmoid(Params.VDD-y,Params.M)/Params.tau_tr-...
    sigmoid(Params.SpikeThresh-Vmem,Params.M)*sigmoid(y,Params.M)/Params.tau_tf;
end

function yout = sigmoid(xin,M)
yout = 1./(1+exp(-M*xin));
end