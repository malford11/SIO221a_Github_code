function [f_vector] = frqvector(data,dt_min)
%FRQVECTOR creates frequency vector in cycles per day
N = length(data);
dt=(dt_min/60)/24; %in cycles per day
T=dt*N;
dfrq=1/T;
fn=1/2/dt;
f_vector=-fn:dfrq:(fn-dfrq);
end

