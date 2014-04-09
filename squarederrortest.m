%this is a test of using an anonymous error function rather than a hardcoded square
%error function.
%it is only about 1% slower and saves having to write a new error function
%for each new function definition.



cd('/Users/skato/Desktop/Dropbox/ZimmerLab/BargmannProjects/Christine/AIAunc103/N2-nv');
load('lnmresults.mat');
pretruncatetime=2; dt=0.05;
kern=lnmstruct.h_est_causal(1:end-pretruncatetime/dt);

starting_paramguess=[2 2 20 -5 -15];
inp=dt*(0:(length(kern)-1))';
options=optimset('Display','off');  %'iter' 'final' 'none'


tic
for i=1:100
[pm,fval]=fminsearch(@lti3odifffitwangp2,starting_paramguess,options,inp,kern);
end
toc

tic
for i=1:100
[pm,fval]=fminsearch(@squarederror,starting_paramguess,options,@lti3odifffnwangp2,inp,kern);
end
toc
