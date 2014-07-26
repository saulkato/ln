function out=deriv(traces)
% out=DERIV(traces)
%
% take time derivative of an array of traces
% does no smoothing, so it will increase noise
%
% use with fastsmooth.m to smooth traces
%
% Saul Kato
% version 1.0
% first version 12/9/09
% latest version 3/9/10

    n=size(traces,1);
    newtraces=0*traces;
    
    for i=1:size(traces,2)

        newtraces(1,i)=traces(2,i)-traces(1,i);
        newtraces(n,i)=traces(n,i)-traces(n-1,i);
        newtraces(2:n-1,i)=(traces(3:n,i)-traces(1:(n-2),i))/2;

    end

    out=newtraces;
end