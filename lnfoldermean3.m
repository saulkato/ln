function lnmstruct=lnfoldermean2(foldername,options,quickflag)   

    if nargin<1
        foldername=pwd;
    end
    
    if nargin<2
        options=0;
    end
    
    if nargin<3
        quickflag=[];
    end
    
    if ischar(options)   %look for quickflags
        quickflag=options;
        options=0;
    end
     
    if isfield(options,'writefolder')
        writefolder=options.writefolder;
    else
        writefolder=pwd;
    end
    
    if isfield(options,'geno')
        geno=options.geno;
    else
        geno='';
    end
    
    if isfield(options,'trialestflag')
        trialestflag=options.trialestflag;
    else
        trialestflag=1;
    end
    
    if isfield(options,'kerneltime')
        kerneltime=options.kerneltime;
    else
        kerneltime=26;
    end
    
    if isfield(options,'plotrow')
        plotrow=options.plotrow;
    else
        plotrow=1;
    end
    
    if isfield(options,'plotnr')
        nr=options.plotnr;
    else
        nr=2;
    end    
    
    if isfield(options,'plots')
        plots=options.plots;
    else
        plots=1;
    end
    
    if isfield(options,'subplotflag')
        subplotflag=options.subplotflag;
    else
        subplotflag=0;
    end
    
    if isfield(options,'plotnxl')
        nxlim=options.plotnxl;
    else
        nxlim=[-.1 .1];
    end
    
    if isfield(options,'plotnyl')
        nylim=options.plotnyl;
    else
        nylim=[-1 2];
    end
    
    if isfield(options,'iter')
        iter=options.iter;
    else
        iter=2;
    end

    if isfield(options,'crossvalfilt')
        crossvalfilt=options.crossvalfilt;
    else
        crossvalfilt=0;
    end
    
    if isfield(options,'odepretruncatetime')
        odepretruncatetime=options.odepretruncatetime;
    else
        odepretruncatetime=2;
    end
    
    if isfield(options,'prekerneltime')
        prekerneltime=options.prekerneltime;
    else
        prekerneltime=2;
    end
    
    if isfield(options,'nltype')
        nltype=options.nltype;
    else
        nltype='poly5'; %'pows';
    end
    
    if isfield(options,'svd_len')
        svd_len=options.svd_len;
    else
        svd_len=50; %150
    end
    
    if isfield(options,'meanquantizeflag')
        meanquantizeflag=options.meanquantizeflag;
    else
        meanquantizeflag=false;  
    end
    
    
    if isfield(options,'odefitflag')
        odefitflag=options.odefitflag;
    else
        odefitflag=false;
    end
    
    if isfield(options,'odeplotflag')
        odeplotflag=options.odeplotflag;
    else
        odeplotflag=false;
    end
    
    if isfield(options,'odemodeltype')
        odemodeltype=options.odemodeltype;
    else
        odemodeltype= '3o5wangp2'; %'3o5gc';  %AWC 
    end   
    
    if isfield(options,'odefitstarting')
        odefitstarting=options.odefitstarting;
    elseif exist('odeparams.mat','file')
        load('odeparams.mat');
        odefitstarting=odeseedparams;
    else
        odefitstarting=[.25 .55 5 -.25/2 -.035/2];  %AWC 
    end
    
    if isfield(options,'saveflag')
        saveflag=options.saveflag;
    else
        saveflag=1;
    end

    if isfield(options,'bleachflag')
        bleachflag=options.bleachflag;
    else
        bleachflag=4;
    end
    
    if isfield(options,'comparisonkernel')
        comparisonkernelflag=1;
    else
        comparisonkernelflag=0;
    end
    
    if isfield(options,'comparisonlegend')
        comparisonlegendflag=1;
    else
        comparisonlegendflag=0;
    end
    
%run lnest on mean i/o

[dye,dyetime,filename]=loaddyefolder(foldername);
[xpos,ypos,traces,time,filename,metadata,baseline]=loadlogfolder(foldername);

dt=metadata.dt;

if isfield(metadata,'endbleachrange')
    endbleachrange=metadata.endbleachrange;
else
    disp('no endbleachrange in metafile.  guessing [end-300:end-1].');
    endbleachrange=[size(traces,1)-300 size(traces,1)-1  ];
end
   

[dye_c,time_c]=timecorrect(dye,dyetime,dt);
if isfield(metadata,'starttimes')
    starttimes=metadata.starttimes;
    disp('taking starttimes from metafile');
else
    starttimes=findtriggers(quantizedye(dye_c))
end
rangeflag=0;rng=[metadata.startframe metadata.endframe];
timeregisterflag=1;timecorrectflag=1;cooperativityflag=0;normflag=0;snrmaxflag=0;dedecayflag=0;fminfac=0.9;

%balance trace lengths
if size(dye_c,1)>size(traces,1)
    disp('lnfoldermean: trace and dye lengths differ. cropping.');
    dye_c(end,:)=[];
    dyetime(end,:)=[];
elseif size(dye_c,1)<size(traces,1)
    disp('lnfoldermean: trace and dye lengths differ. cropping.');
    traces(end,:)=[];
end

[traces,flagstr,tau_bleach,~,~,time_corr]=processtraces(traces,dyetime,rng,timecorrectflag,timeregisterflag,bleachflag,rangeflag,cooperativityflag,normflag,starttimes,snrmaxflag,dedecayflag,endbleachrange,fminfac,metadata.dt);

if timeregisterflag
    dye_tc=timeregister(dye_c,starttimes,150);
else
    dye_tc=dye_c;
end

meaninp=mean(dye_tc,2);
meaninpq=quantizedye(mean(dye_tc,2));
meanqinp=mean(quantizedye(dye_tc),2);
meanoutp=mean(traces,2);

%esimate trial mean kernels
disp('lnfoldermean: estimating trial-mean kernel');
h_len=kerneltime/dt;  %number of samples in kernel
h_len2=round(prekerneltime/dt);
detrendtype='constant';
derivflag=0;

%lnmstruct=lnest8(meanqinp,meanoutp,h_len,svd_len,rng,1,0,detrendtype,nltype,h_len2,iter,derivflag); %1/12/2014

if meanquantizeflag
    lnmstruct=lnest9(meaninpq,meanoutp,h_len,svd_len,rng,1,0,detrendtype,nltype,h_len2,iter,derivflag);
    lnmstruct.meaninp=meaninpq;
else
    lnmstruct=lnest9(meaninp,meanoutp,h_len,svd_len,rng,1,0,detrendtype,nltype,h_len2,iter,derivflag);
    lnmstruct.meaninp=meaninp;
end

lnmstruct.meanoutp=meanoutp;



if odefitflag
    options_ode.modeltype=odemodeltype;
   % [ltikernel,ltikerneltype,ltikernel_pm]=ltikernelfit(lnmstruct.h_est_causal(1:end-2/dt),modeltype,dt,odefitstarting);
    options_ode.starting=odefitstarting;
    options_ode.dt=0.05;
    options_ode.pretruncatetime=odepretruncatetime;
    options_ode.subplotflag=1;
    options_ode.plots=1;    
    lnmstruct.ode=odefit2(lnmstruct,options_ode);
    lnmstruct.ode.ltikerneltv=mtv(lnmstruct.ode.ltikernel,options_ode.dt);
    
    crossvalopts.flag=1;
    crossvalopts.hhat=lnmstruct.ode.ltikernel;
    h_len_trunc=(kerneltime-options_ode.pretruncatetime)/dt;
    
    lnmstruct.ode.lns=lnest9(meaninp,meanoutp,h_len_trunc,svd_len,rng,1,0,detrendtype,nltype,0,iter,derivflag,crossvalopts);
end




if plots
    nc=4;
    if ~subplotflag 
        fig1handle=figure('Position',[0 0 1200 800]);
    end
    
    %plot L
    subplot(nr,nc,nc*(plotrow-1)+(1:2));
    ylabel(geno);
    xlim(dt*[-h_len2 h_len]);
    ylim([min(lnmstruct.h_est) max(lnmstruct.h_est)]*1.2);
    hline(0);
    hold on;
    plot(lnmstruct.iv*dt,lnmstruct.h_est);
    %plot(lnmstruct.iv*dt,lnmstruct.h_est_pass1,'--','Color','g');  
    if comparisonkernelflag
        refLNstruct=load(options.comparisonkernel);
        refkernel=refLNstruct.lnmstruct.h_est;
        plot(lnmstruct.iv(1:length(refkernel))*dt,refkernel/max(abs(refkernel))*max(abs(lnmstruct.h_est)),'r');
        if comparisonlegendflag
            legend(options.comparisonlegend);             
        else
            legend({'this expt','reference kernel'});             
        end
       
    elseif strcmp(quickflag,'refash')  %plot against ash reference kernel
        refLNstruct=load('/Users/skato/Desktop/Dropbox/WormData2/PaperData/newbleach/newmeankernels.mat');
        refkernel=refLNstruct.ln_sm2{2}.h_est;
        plot(lnmstruct.iv(1:length(refkernel))*dt,refkernel/max(abs(refkernel))*max(abs(lnmstruct.h_est)),'r');
        legend({'this expt','N2 ASH reference'});
    elseif strcmp(quickflag,'refawc')  %plot against ash reference kernel
        refLNstruct=load('/Users/skato/Desktop/Dropbox/MyPaper/newfigs/draft13/ALLlns5.mat');
        refkernel=refLNstruct.lns{1}.h_est;
        plot(lnmstruct.iv(1:length(refkernel))*dt,refkernel/max(abs(refkernel))*max(abs(lnmstruct.h_est)),'r');
        legend({'this expt','N2 AWC reference'});
    end
    
    
    if odeplotflag
        
        plot(lnmstruct.ode.extraptv,lnmstruct.ode.ltikernelextrap,'g');
        
    end
    
    
    vline(0);
    xlabel('lag (s)');
    grid on;
    set(gca,'XMinorTick','on');
    
    subplot(nr,nc,nc*(plotrow-1)+3);
    %plot scatter
    plot(lnmstruct.interm,lnmstruct.trueout,'.','MarkerSize',1,'Color',[0.5 0.5 0.5]);
    hold on;
    irv=linspace(lnmstruct.intermrng(1),lnmstruct.intermrng(2));
    plot(irv,nfn(lnmstruct.npm,irv,nltype),'b');
    xlim(1.1*[min(lnmstruct.interm) max(lnmstruct.interm)]);
    ylim(1.1*[min(lnmstruct.trueout) max(lnmstruct.trueout)]);
    
    vline(0); 
    xlabel('in');ylabel('out');
    
    
    %plot reference N      
    if comparisonkernelflag
        refnpm=refLNstruct.lnmstruct.npm;
        refnmax=nfn(refnpm,refLNstruct.lnmstruct.intermrng(2),nltype);
        refnmin=nfn(refnpm,refLNstruct.lnmstruct.intermrng(1),nltype); 
        refirv=linspace(refLNstruct.lnmstruct.intermrng(1),refLNstruct.lnmstruct.intermrng(2),100);
        plot(irv,(nfn(refnpm,refirv,nltype)-refnmin)/(refnmax-refnmin)*(nfn(lnmstruct.npm,irv(end),nltype)-nfn(lnmstruct.npm,irv(1),nltype))+nfn(lnmstruct.npm,irv(1),nltype),'r');

    elseif strcmp(quickflag,'refash')  %plot against ash reference kernel
        nash=load('/Users/skato/Desktop/Dropbox/WormData3/YifanNaClandDHCA/N2Gly/lnmresults.mat');
        refashnpm=nash.lnmstruct.npm;
%        plot(irv,nfn(refashnpm,irv,nltype)/max(nfn(refashnpm,irv,nltype))*max(nfn(lnmstruct.npm,irv,nltype)),'r--');
        refnmax=nfn(refashnpm,nash.irv(end),nltype);
        refnmin=nfn(refashnpm,nash.irv(1),nltype);
        plot(irv,(nfn(refashnpm,nash.irv,nltype)-refnmin)/(refnmax-refnmin)*(nfn(lnmstruct.npm,irv(end),nltype)-nfn(lnmstruct.npm,irv(1),nltype))+nfn(lnmstruct.npm,irv(1),nltype),'r');
    elseif strcmp(quickflag,'refawc')  %plot against awc reference kernel
        refawcnpm=refLNstruct.lns{1}.npm;
        %size(refawcnpm)
        refnmax=nfn([refawcnpm refLNstruct.lns{1}.intermrng(1)],refLNstruct.lns{1}.intermrng(2),nltype)
        refnmin=nfn([refawcnpm refLNstruct.lns{1}.intermrng(1)],refLNstruct.lns{1}.intermrng(1),nltype)
        refirv=linspace(refLNstruct.lns{1}.intermrng(1),refLNstruct.lns{1}.intermrng(2),100);
        %plot(irv,nfn(refawcnpm,refirv,nltype),'r--');
        nfn(lnmstruct.npm,irv(end),nltype)
        plot(irv,(nfn([refawcnpm refLNstruct.lns{1}.intermrng(1)],refirv,nltype)-refnmin)/(refnmax-refnmin)*(nfn(lnmstruct.npm,irv(end),nltype)-nfn(lnmstruct.npm,irv(1),nltype))+nfn(lnmstruct.npm,irv(1),nltype),'r');
    end

    if odefitflag
       [odemodel_simout,odevaf1,odevaf2,model_interm,sse2]=lneval2(lnmstruct.inp,lnmstruct.ode.ltikernelextrap*lnmstruct.hgain,lnmstruct.npm,lnmstruct.nltype,lnmstruct.trueout,0,1,0,0.5);
    end
    
    %write out information
    subplot(nr,nc,nc*(plotrow-1)+4);
    ax=axis;
    
    if odefitflag
         txt={['n=' num2str(size(traces,2))];...
         ['mean VAF (LN)=',num2str(lnmstruct.vaf2,3),'%'];...
         ['mean VAF (L only)=',num2str(lnmstruct.vaf1,3),'%'];...
         ['L length=',num2str(h_len),' samples'];['SVD comps=',num2str(svd_len)];...   
         ['dt=',num2str(dt),' s'];['frame range= [',num2str(rng(1)),'-',num2str(rng(2)),']'];...
         ['N func=',nltype];['detrend=',detrendtype];...
         ['flags:',flagstr];...
         ['mean VAF (LN-ODE)=',num2str(odevaf2,3),'%'];...
         ['ODE params: ',num2str(lnmstruct.ode.ltikernel_pm,4)];...
         };
    else
         txt={['n=' num2str(size(traces,2))];...
         ['mean VAF (LN)=',num2str(lnmstruct.vaf2,3),'%'];...
         ['mean VAF (L only)=',num2str(lnmstruct.vaf1,3),'%'];...
         ['L length=',num2str(h_len),' samples'];['SVD comps=',num2str(svd_len)];...   
         ['dt=',num2str(dt),' s'];['frame range= [',num2str(rng(1)),'-',num2str(rng(2)),']'];...
         ['N func=',nltype];['detrend=',detrendtype];...
         ['flags:',flagstr];...
         };
     
    text(ax(1),ax(4),txt,'VerticalAlignment','top');
    axis off;
    box off;
           
%     if odefitflag
%         subplot(nr,nc,nc*(plotrow-1)+(5:6));
%         extraptv=0:dt:2*lnmstruct.iv(end)*dt;        
%         xlim([dt*-h_len2 extraptv(end)]);
%         ylim([min(lnmstruct.h_est) max(lnmstruct.h_est)]*1.2);
%         vline(0);
%         hline(0);
%         hold on;
% 
%         plot(lnmstruct.iv*dt,lnmstruct.h_est,'k');
%         plot(extraptv,lti3odifffn5(lnmstruct.ode.ltikernel_pm,extraptv),'b');
%         plot(mtv(lnmstruct.ode.ltikernel,dt),lti3odifffn5(odefitstarting,mtv(lnmstruct.ode.ltikernel,dt)),'Color',color('lr'));
%     end
end 

if ~subplotflag
    
    %plot mean input and output traces
    subplot(nr,nc,nc*(plotrow-1)+(5:8));
    otv=mtv(lnmstruct.meanoutp,dt);
    plot(mtv(lnmstruct.meaninp,dt),lnmstruct.meaninp/2-1/2,'r');
    hold on;
    plot(otv,lnmstruct.meanoutp,'Color',color('gray'));
    
    %plot(otv(rng(1):rng(end)),lnmstruct.simtrace-lnmstruct.shift,'b');

    plot(otv(rng(1):rng(end)),lnmstruct.simtrace+mean(lnmstruct.meanoutp),'b');

    
    if odeplotflag
        
        [odemodel_simout,odevaf1,odevaf2,model_interm,sse2]=lneval2(lnmstruct.inp,lnmstruct.ode.ltikernelextrap*lnmstruct.hgain,lnmstruct.npm,lnmstruct.nltype,lnmstruct.trueout,0,1,0,0.5);

        plot(otv(rng(1):rng(end)),odemodel_simout+mean(lnmstruct.meanoutp),'g');
        
    end

    %plot reference output        
    if strcmp(quickflag,'refash')  %plot against ash reference kernel
            plot(mtv(nash.meanoutp,dt),nash.meanoutp,'r');
            if comparisonlegendflag
                legend(options.comparisonlegend);
            else
                legend({'this expt','N2 ASH reference'});
            end
    end

    ylim([-0.5,  ceil(2*max(lnmstruct.meanoutp))/2]);
    xlim([0 otv(end)]);
    ylabel('\DeltaF/F');xlabel('time (s)');
    title('trial-averaged input and output, and simulation');
    proplot;
    set(gca,'XMinorTick','on');
    legend({'input','true output','model output','LN-ODE model output'});

    mtit(strrep(foldername,'_','\_'));
end

pdfflag=1;
if pdfflag==1
    figtitle2=['LN_mean_' num2str(prekerneltime) '+' num2str(kerneltime) 's_' num2str(svd_len) 'svd' ];
    tit_f1_plain=[filename{1}(1:end-11) filename{1}(end-8:end-11)];
        rangestr=['[' num2str(rng(1)) '-' num2str(rng(2)) ']'];
    save2pdf([writefolder '/' figtitle2 '_' tit_f1_plain flagstr rangestr '.pdf'],gcf,600);
end

drawnow;


if trialestflag
    %estimate individual trial kernels
    disp('lnfoldermean: estimating individual trial kernels');
    for n=1:size(traces,2)
        disp('lnfoldermean: estimating individual trial kernel');
        lnmstruct.tracelns{n}=lnest8(dye_tc(:,n),traces(:,n),h_len,svd_len,rng,1,0,detrendtype,nltype,h_len2,iter,derivflag);
        %run mean kernel on individual trials
        crossvalopts.flag=1;
        crossvalopts.hhat=lnmstruct.h_est;  
        disp('lnfoldermean: estimating mean kernel used on indiv. trial');
        lnmstruct.tracelnsM{n}=lnest8(quantizedye(dye_tc(:,n)),traces(:,n),h_len,svd_len,rng,1,0,detrendtype,nltype,h_len2,iter,derivflag,crossvalopts);

        if length(crossvalfilt)>1
            crossvalopts.hhat=options.crossvalfilt;
            disp('lnfoldermean: estimating crossval. kernel used on indiv. trial');
            lnmstruct.tracelnsCV{n}=lnest8(quantizedye(dye_tc(:,n)),traces(:,n),h_len,svd_len,rng,1,0,detrendtype,nltype,h_len2,iter,derivflag,crossvalopts);
        end
    end 
    
    
    
    
end

if odefitflag
    if trialestflag
        %do ode est on each trial
        for n=1:size(traces,2)
            lnmstruct.tracelns{n}.ode=odefit2(lnmstruct.tracelns{n},options_ode);
            
            %get the performance of the trialODE filter
            crossvalopts.flag=1;
            crossvalopts.hhat=lnmstruct.tracelns{n}.ode.ltikernel;    
            lnmstruct.tracelnsO{n}=lnest8(quantizedye(dye_tc(:,n)),traces(:,n),h_len_trunc,svd_len,rng,1,0,detrendtype,nltype,0,iter,derivflag,crossvalopts);

            %get the performance of the meanODE filter
            crossvalopts.hhat=lnmstruct.ode.ltikernel;
            lnmstruct.tracelnsMO{n}=lnest8(quantizedye(dye_tc(:,n)),traces(:,n),h_len_trunc,svd_len,rng,1,0,detrendtype,nltype,0,iter,derivflag,crossvalopts);

            %
%             if crossvalfilt
%                lnmstruct.tracelnsCV{n}=lnest8(quantizedye(dye_tc(:,n)),traces(:,n),crossvalfilt,svd_len,rng,1,0,detrendtype,nltype,0,iter,derivflag,crossvalopts);
%             end
            
        end
    end
end
   

%plot individual trials
if trialestflag
    
    if plots
        if odefitflag
            nc=5;
        else
            nc=3;
        end
        if ~subplotflag 
            figure('Position',[0 0 800 1200]);
        end

        nr=size(lnmstruct.tracelns,2);
        for i=1:nr
            plotrow=i;

            subplot(nr,nc,nc*(plotrow-1)+(1:2));

            hold on;
            plot(lnmstruct.tracelns{i}.iv*dt,lnmstruct.tracelns{i}.h_est);

            plotpass1flag=0;
            if plotpass1flag
                plot(lnmstruct.tracelns{i}.iv*dt,lnmstruct.tracelns{i}.h_est_pass1,'--','Color','g');
            end

            if comparisonkernelflag
                refkernelindiv=refLNstruct.lnmstruct.tracelns{i}.h_est;
                plot(lnmstruct.iv(1:length(refkernelindiv))*dt,refkernelindiv/max(abs(refkernelindiv))*max(abs(lnmstruct.tracelns{i}.h_est)),'r');
                if comparisonlegendflag
                    legend(options.comparisonlegend);
                else
                    legend({'this expt','reference kernel'});  
                end
            elseif strcmp(quickflag,'refash')  %plot against ash reference kernel
                refLNstruct=load('/Users/skato/Desktop/Dropbox/WormData2/PaperData/newbleach/newmeankernels.mat');
                refkernel=refLNstruct.ln_sm2{2}.h_est;
                plot(lnmstruct.iv(1:length(refkernel))*dt,refkernel/max(abs(refkernel))*max(abs(lnmstruct.tracelns{i}.h_est)),'r');
                %legend({'this expt','N2 ASH'});
            elseif strcmp(quickflag,'refawc')  %plot against awc reference kernel
                refLNstruct=load('/Users/skato/Desktop/Dropbox/MyPaper/newfigs/draft13/ALLlns5.mat');
                refkernel=refLNstruct.lns{1}.h_est;
                plot(lnmstruct.iv(1:length(refkernel))*dt,refkernel/max(abs(refkernel))*max(abs(lnmstruct.tracelns{i}.h_est)),'r');
                %legend({'this expt','N2 AWC'});
            end
            
            ylabel(geno);
            xlim(dt*[-h_len2 h_len]);
            ylim([min(lnmstruct.tracelns{i}.h_est) max(lnmstruct.tracelns{i}.h_est)]*1.1);

            vline(0);
            hline(0);

            if i==nr xlabel('lag (s)'); end
                
%             subplot(nr,nc,nc*(plotrow-1)+3);
%             %plot scatter fixedscale
%             plot(lnmstruct.tracelns{i}.interm,lnmstruct.tracelns{i}.trueout,'.','MarkerSize',1);
%             hold on;
%             irv=linspace(lnmstruct.tracelns{i}.intermrng(1),lnmstruct.tracelns{i}.intermrng(2));
%             plot(irv,nfn(lnmstruct.tracelns{i}.npm,irv,lnmstruct.tracelns{i}.nltype),'r');
%             xlim(nxlim);
%             ylim(nylim);

            subplot(nr,nc,nc*(plotrow-1)+3);
            %plot scatter 
            plot(lnmstruct.tracelns{i}.interm,lnmstruct.tracelns{i}.trueout,'.','MarkerSize',1,'Color',[0.5 0.5 0.5]);
            hold on;
            plot(irv,nfn(lnmstruct.tracelns{i}.npm,irv,lnmstruct.tracelns{i}.nltype),'b');
            hline(0);
            vline(0);
                xlim(1.1*[min(lnmstruct.tracelns{i}.interm) max(lnmstruct.tracelns{i}.interm)])
            ylim(1.1*[min(lnmstruct.tracelns{i}.trueout) max(lnmstruct.tracelns{i}.trueout)]);
            
            if comparisonkernelflag
                refnpm=refLNstruct.lnmstruct.tracelns{i}.npm;
                refnmax=nfn(refnpm,refLNstruct.lnmstruct.tracelns{i}.intermrng(2),nltype);
                refnmin=nfn(refnpm,refLNstruct.lnmstruct.tracelns{i}.intermrng(1),nltype); 
                refirv=linspace(refLNstruct.lnmstruct.tracelns{i}.intermrng(1),refLNstruct.lnmstruct.tracelns{i}.intermrng(2),100);
                plot(irv,(nfn(refnpm,refirv,nltype)-refnmin)/(refnmax-refnmin)*(nfn(lnmstruct.tracelns{i}.npm,irv(end),nltype)-nfn(lnmstruct.tracelns{i}.npm,irv(1),nltype))+nfn(lnmstruct.tracelns{i}.npm,irv(1),nltype),'r');
            elseif strcmp(quickflag,'refash')  %plot against ash reference kernel
                nash=load('/Users/skato/Desktop/Dropbox/WormData3/YifanNaClandDHCA/N2Gly/lnmresults.mat');
                refashnpm=nash.lnmstruct.npm;
        %        plot(irv,nfn(refashnpm,irv,nltype)/max(nfn(refashnpm,irv,nltype))*max(nfn(lnmstruct.npm,irv,nltype)),'r--');
                refnmax=nfn(refashnpm,nash.irv(end),nltype);
                refnmin=nfn(refashnpm,nash.irv(1),nltype);
                plot(irv,(nfn(refashnpm,nash.irv,nltype)-refnmin)/(refnmax-refnmin)*(nfn(lnmstruct.npm,irv(end),nltype)-nfn(lnmstruct.npm,irv(1),nltype))+nfn(lnmstruct.npm,irv(1),nltype),'r');
            elseif strcmp(quickflag,'refawc')  %plot against awc reference kernel
                refawcnpm=refLNstruct.lns{1}.npm;
                %size(refawcnpm)
                refnmax=nfn([refawcnpm refLNstruct.lns{1}.intermrng(1)],refLNstruct.lns{1}.intermrng(2),nltype)
                refnmin=nfn([refawcnpm refLNstruct.lns{1}.intermrng(1)],refLNstruct.lns{1}.intermrng(1),nltype)
                refirv=linspace(refLNstruct.lns{1}.intermrng(1),refLNstruct.lns{1}.intermrng(2),100);
                %plot(irv,nfn(refawcnpm,refirv,nltype),'r--');
                nfn(lnmstruct.npm,irv(end),nltype)
                plot(irv,(nfn([refawcnpm refLNstruct.lns{1}.intermrng(1)],refirv,nltype)-refnmin)/(refnmax-refnmin)*(nfn(lnmstruct.tracelns{i}.npm,irv(end),nltype)-nfn(lnmstruct.tracelns{i}.npm,irv(1),nltype))+nfn(lnmstruct.tracelns{i}.npm,irv(1),nltype),'r');
            end

      
            if odefitflag
                subplot(nr,nc,nc*(plotrow-1)+(4:5));
                extraptv=0:dt:2*lnmstruct.tracelns{i}.iv(end)*dt;        
                xlim([dt*-h_len2 extraptv(end)]);
                ylim([min(lnmstruct.tracelns{i}.h_est) max(lnmstruct.tracelns{i}.h_est)]*1.2);
                vline(0);
                hline(0);
                hold on;
                %need to fix these regerences
                plot(lnmstruct.tracelns{i}.iv*dt,lnmstruct.tracelns{i}.h_est,'k');
                plot(extraptv,lti3odifffn5(lnmstruct.ode.ltikernel_pm,extraptv),'b');
                plot(mtv(lnmstruct.ode.ltikernel,dt),lti3odifffn5(odefitstarting,mtv(lnmstruct.ode.ltikernel,dt)),'Color',color('lr'));
            end


            ax=axis;
            text(ax(2),ax(4),[num2str(lnmstruct.tracelns{i}.vaf2,3) '% '],'VerticalAlignment','top','HorizontalAlignment','right');

            end
        end
    
   mtit(strrep(foldername,'_','\_'));

   if pdfflag==1
    figtitle3=['LN_mean_indiv_' num2str(prekerneltime) '+' num2str(kerneltime) 's_' num2str(svd_len) 'svd' ];
    tit_f1_plain=[filename{1}(1:end-11) filename{1}(end-8:end-11)];
        rangestr=['[' num2str(rng(1)) '-' num2str(rng(2)) ']'];
    save2pdf([writefolder '/' figtitle3 '_' tit_f1_plain flagstr rangestr '.pdf'],gcf,600);
end

end

drawnow;

lnmstruct.options=options;
lnmstruct.rng=rng;


if saveflag
    save([foldername '/lnmresults.mat']);
end

end