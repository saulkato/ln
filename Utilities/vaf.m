function v = vaf(x,y)
% compute variance accounted for between two signals
% v = vaf(x,y)
%
% uses built-in matlab function STD
%
% first variable is the true signal
% 20131204 modified to accept multi-d signals as matrix of column vectors

residual=x-y;
v=100*(1-sum(std(residual).^2)/sum(std(x).^2));