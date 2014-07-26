function odemodel=odefit2(lnmstruct,options)
%fit an odemodel to an lns struct
%Saul Kato
%20130110


%odemodel={'3o5gc','3o5wang','3o5wang2','3o5wangp2_lneval2';};

odeModelArray(1).index=1;
odeModelArray(1).name='3o5wang2';
odeModelArray(1).numParams=5;
odeModelArray(1).paramLabels={'a','b','c','d','e'};
odeModelArray(1).func=@lti3odifffnwang2;
odeModelArray(1).funcDeconvolved=@lti3ofn5exact;
odeModelArray(1).paramMax=[10,40,40,100,1];
odeModelArray(1).paramMin=[0.001,0.001,0.001,.00001,-1];
odeModelArray(1).standardLossFuncFlag=true;
odeModelArray(1).fitseed=[10 1.051 37 .026 -.019];

odeModelArray(2).index=2;
odeModelArray(2).name='3o5wangp2';
odeModelArray(2).numParams=5;
odeModelArray(2).paramLabels={'Trise','TfallF','TfallS','kas/kaf','kas+kaf'};
odeModelArray(2).func=@lti3odifffnwangp2;
odeModelArray(2).funcDeconvolved=@lti3ofn5exact;
odeModelArray(2).paramMax=[10,40,40,100,1];
odeModelArray(2).paramMin=[0.001,0.001,0.001,.00001,-1];
odeModelArray(2).standardLossFuncFlag=true;
odeModelArray(2).fitseed=[10 1.051 37 .026 -.019];


odeModelArray(3).index=3;
odeModelArray(3).name='4o2parallel';
odeModelArray(3).numParams=6;
odeModelArray(3).paramLabels={'TriseF','TriseS','TfallF','TfallS','kas/kaf','kas+kaf'};
odeModelArray(3).func=@lti4o2parallelgc;
odeModelArray(3).funcDeconvolved=@lti4o2parallel;
odeModelArray(3).paramMax=[30,30,30,40,100,1];
odeModelArray(3).paramMin=[0.001,0.001,0.001,0.001,.00001,-1];
odeModelArray(3).standardLossFuncFlag=true;
odeModelArray(3).fitseed=[20 10 1.051 37 .026 -.019];


%put odeModel names into a cell array
odeModels={};
for i=1:length(odeModelArray)
    odeModels=[odeModels odeModelArray(i).name];
end

if nargin<2
    options=0;
end

if ischar(lnmstruct)
    if strcmp(lnmstruct,'getmodels');
        odemodel=odeModels;
        return;
    elseif strcmp(lnmstruct,'getmodel');
        if ischar(options)            
            index = find(strcmp(odeModels, options));
        else
            index = options;
        end
        if ~isempty(index)
            odemodel=odeModelArray(index);
        else
            disp('no model found.');
            beep;pause(.1);beep;
        end
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

%select odemodel
thisModelIndex=find(strcmp(odeModels,modeltype));
thisModel=odeModelArray(thisModelIndex);

if isfield(options,'dt')
    dt=options.dt;
else
    dt=0.05;
end

if isfield(options,'pretruncatetime')
    pretruncatetime=options.pretruncatetime;
else
    pretruncatetime=0;
end

%[ltikernel,ltikerneltype,ltikernel_pm]=ltikernelfit(lns.h_est_causal(1:end-pretruncatetime/dt)/sum(lns.h_est_causal(1:end-pretruncatetime/dt)),modeltype,dt,starting);
[ltikernel,ltikerneltype,ltikernel_pm]=ltikernelfit2(lnmstruct.h_est_causal(1:end-pretruncatetime/dt),thisModel.func,dt,starting,lnmstruct,thisModel.standardLossFuncFlag);

odemodel.ltikernel=ltikernel;
odemodel.ltikerneltype=ltikerneltype;
odemodel.ltikernel_pm=ltikernel_pm;
extraptv=0:dt:2*lnmstruct.iv(end)*dt; 
odemodel.extraptv=extraptv;

% if strcmp(modeltype,'3o5wang')
%     odemodel.ltikernelextrap=lti3odifffnwang(ltikernel_pm,extraptv); 
% elseif strcmp(modeltype,'3o5wang2')
%     odemodel.ltikernelextrap=lti3odifffnwang2(ltikernel_pm,extraptv); 
% elseif strcmp(modeltype,'3o5wangp2') || strcmp(modeltype,'3o5wangp2_lneval2')
%     odemodel.ltikernelextrap=lti3odifffnwangp2(ltikernel_pm,extraptv); 
% else
%     odemodel.ltikernelextrap=lti3odifffn5exact(ltikernel_pm,extraptv); 
% end

odemodel.ltikernelextrap=odeModelArray(thisModelIndex).func(ltikernel_pm,extraptv);
odemodel.ltikernelextrapdc=odeModelArray(thisModelIndex).funcDeconvolved(ltikernel_pm,extraptv);

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
    plot(mtv(ltikernel,dt),odeModelArray(thisModelIndex).funcDeconvolved(starting,mtv(ltikernel,dt)),'Color',color('r'));
end

end %function

