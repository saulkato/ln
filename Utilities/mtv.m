function timevec=mtv(trace,dt)
%timevec=mtv(trace,dt)
%
%maketimevec from t=0, dt steps, to trace length)
%Saul Kato 10/3/12
%

if nargin<2
    dt=0.05;
end

timevec=dt*(0:length(trace)-1);