function [pm shift]=nest(interm,outp,nltype,Starting)
%[pm shift]=nest(interm,outp)
%estimate a nonlinear instantaneous i/o relationship
%
%currently only works for pow3
%f(x)=a(x-x0)^p
%params=[a p x0];
%x0 is fixed to be min(interm);
%a and p are estimated using pow2fit and pow2fn
%
%Saul Kato
%120710

%if interm and outp are different lengths, truncate interm
%pinterm=interm(1:length(outp));

if nargin<3
    nltype='pow';
end

options=optimset('Display','off');  %'iter' 'final' 'none'

%L-N curvefit, 
if strcmp(nltype,'sr')  %sigmoid soft rectifier
    Starting=[1 2 0 1];
    %curve fit
    pm=fminsearch(@softrectifierfit,Starting,options,interm(1:length(outp))',outp);
elseif strcmp(nltype,'pow')  %power fn  
    disp('pow is deprecated.');
    Starting=[1 1 1];
    %curve fit
    shift=min(min(interm));
    interm_cut=interm(1:length(outp),:);
    pm=fminsearch(@powerfit,Starting,options,interm_cut'-shift,outp(:));  %+.0001;
    pm=[pm shift];
    
elseif strcmp(nltype,'pows')  %power fn  
    shift=min(min(interm));
    Starting=[1 1 1];
    %curve fit
    interm_cut=interm(1:length(outp),:);
    pm=fminsearch(@powerfitnoshift,Starting,options,interm_cut'-shift,outp(:));  %+.0001;
    pm=[pm shift];
    
elseif strcmp(nltype,'powns')  %power fn without shift
    Starting=[1 1 3];
    %curve fit
    interm_cut=interm(1:length(outp),:);
    pm=fminsearch(@powerfit,Starting,options,interm_cut',outp(:));
elseif strcmp(nltype,'powfix') %power fn with fixed power
    Starting=[1 1];
    %curve fit
    shift=min(min(interm));
    pm=fminsearch(@powerfit,Starting,options,(interm(1:length(outp))'-shift),outp); %+.0001;
elseif strcmp(nltype,'pow3')
    shift=min(min(interm));
    Starting=[1 1];  
    pm=fminsearch(@pow2fit,Starting,options,(interm'-shift),outp);
elseif strcmp(nltype,'powh')
    if nargin<4
        Starting=[1 0 2.3 min(min(interm))];  
    end
    interm_cut=interm(1:length(outp),:);
    pm=fminsearch(@powh,Starting,options,interm_cut,outp);    
elseif strcmp(nltype,'powhfix')
    if nargin<4
        Starting=[20 0 .13];  
    end
    interm_cut=interm(1:length(outp),:);
    pm=fminsearch(@powhfix,Starting,options,interm_cut',outp);
elseif strcmp(nltype,'powhfixnooff')
    if nargin<4
        Starting=[1 min(min(interm))];  
    end
    interm_cut=interm(1:length(outp),:);
    pm=fminsearch(@powhfixnooff,Starting,options,interm_cut',outp);
    
    
  %  figure;
  %  plot(interm,outp,'b');
  %  hold on;
  %  z=linspace(min(interm),max(interm),100);
  %  plot(z,pow2fn(z-min(interm),pm),'r');
  %  figure;
    
else  %'poly5'  %degree-5 polynomial
    warning('OFF','MATLAB:polyfit:RepeatedPointsOrRescale');
    pm = polyfit(interm(1:length(outp)),outp,5);
    warning('ON','MATLAB:polyfit:RepeatedPointsOrRescale');
end

if ~exist('shift') shift=0; end
