function dydt = Triangle(t,y,Params,Vmem)
% Triangle Generator ODE (for generating Rall Alpha)
dydt = sigmoid(Vmem-Params.SpikeThresh,Params.M)*sigmoid(Params.VDD-y,Params.M)/Params.tau_tr-...
    sigmoid(Params.SpikeThresh-Vmem,Params.M)*sigmoid(y,Params.M)/Params.tau_tf;
end

function yout = sigmoid(xin,M)
yout = 1./(1+exp(-M*xin));
end