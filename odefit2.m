function odemodel=odefit2(lnmstruct,options)
%fit an odemodel to an lns struct
%Saul Kato
%20130110

if nargin<2
    options=0;
end

if ischar(lnmstruct)
    if strcmp(lnmstruct,'getmodels');
        odemodel={'3o5gc','3o5wang','3o5wang2','3o5wangp2_lneval2';};
        return;
    end
end

if isfield(options,'subplotflag')
    subplotflag=options.subplotflag;
else
    subplotflag=0;
end

if isfield(options,'plots')
    plots=options.plots;
else
    plots=1;
end

if isfield(options,'starting')
    starting=options.starting;
else
    starting=[.5 .55 10 -.025 -.0175];
end

if isfield(options,'modeltype')
    modeltype=options.modeltype;
else
    modeltype='3o5wangp2_lneval2'; %'3o5gc';
end

if isfield(options,'dt')
    dt=options.dt;
else
    dt=0.1;
end

if isfield(options,'pretruncatetime')
    pretruncatetime=options.pretruncatetime;
else
    pretruncatetime=0;
end

%[ltikernel,ltikerneltype,ltikernel_pm]=ltikernelfit(lns.h_est_causal(1:end-pretruncatetime/dt)/sum(lns.h_est_causal(1:end-pretruncatetime/dt)),modeltype,dt,starting);
[ltikernel,ltikerneltype,ltikernel_pm]=ltikernelfit2(lnmstruct.h_est_causal(1:end-pretruncatetime/dt),modeltype,dt,starting,lnmstruct);

odemodel.ltikernel=ltikernel;
odemodel.ltikerneltype=ltikerneltype;
odemodel.ltikernel_pm=ltikernel_pm;
extraptv=0:dt:2*lnmstruct.iv(end)*dt; 
odemodel.extraptv=extraptv;
if strcmp(modeltype,'3o5wang')
    odemodel.ltikernelextrap=lti3odifffnwang(ltikernel_pm,extraptv); 
elseif strcmp(modeltype,'3o5wang2')
    odemodel.ltikernelextrap=lti3odifffnwang2(ltikernel_pm,extraptv); 
elseif strcmp(modeltype,'3o5wangp2') || strcmp(modeltype,'3o5wangp2_lneval2')
    odemodel.ltikernelextrap=lti3odifffnwangp2(ltikernel_pm,extraptv); 
else
    odemodel.ltikernelextrap=lti3odifffn5exact(ltikernel_pm,extraptv); 
end

odemodel.ltikernelextrapdc=lti3ofn5exact(ltikernel_pm,extraptv);

odemodel.starting=starting;

if plots
    if ~(subplotflag) figure('Position',[0 0 800 600]); end;
          
    xlim([dt*-lnmstruct.h_len2 extraptv(end)]);
    ylim([min(lnmstruct.h_est) max(lnmstruct.h_est)]*1.2);
    vline(0);
    hline(0);
    hold on;
    plot(lnmstruct.iv*dt,lnmstruct.h_est,'k');
    plot(extraptv,odemodel.ltikernelextrap,'b');
    plot(extraptv,-abs(odemodel.ltikernelextrapdc/max(abs(odemodel.ltikernelextrapdc))*max(abs(odemodel.ltikernelextrap))),'g');
    plot(mtv(ltikernel,dt),lti3odifffn5exact(starting,mtv(ltikernel,dt)),'Color',color('r'));
end

end %function

