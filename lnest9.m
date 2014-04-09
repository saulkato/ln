function lnstruct=lnest9(input,output,hlen,svdlength,range,estmethod,plots,detrendtype,nltype,hlen2,iter,derivflag,crossvalopts)
%lnstruct=lnest9(input,output,hlen,svdlength,range,estmethod,plots,detrendtype,nltype,hlen2,iter)
%lnstruct fields: h_est npm gain offset shift vaf1 vaf2 interm simtrace hhat2 inp_offset outp_offset
%
%Estimate impulse responses from a matrix of input/output record
%using cross correlation technique
%input= input vector or array (time is 1st dim)
%output= output vector or array
%svdlength= number of components to use for SVD smoothing (only used if
%            estmethod=0)
%range=  2-element vector indicating what subset of i/o
%        record to use for estimation  default: all range
%hlen=length of impulse response, in samples, to estimate
%estmethod= method of estimation
%    0: cross-correlation deconvolution
%    1: direct regression
%
%npm=estimated parameters of nonlinear function
%nltype='poly5' polynomial of degree 5
%       'sr'  soft rectifier (integrated sigmoid)
%       'pow' power law function
%
%Saul Kato
%created 100115
%last update 100817  lnest2 expanded to support different choice of N fn
%version 3 100927    lnest3 added normalization of gain
%version 4 101110    lnest4 offset, gain for different nfns
%          101210   added asymmetric 2-sided estimation
%version 5 111610  estimate from multiple traces
%version 7 020811 trying to debug notch, plus factor lest.m, plus add re-estimation
%version 8 051212 cleanup, struct output
%
%multilnest42->lnest8->nest,lneval2,lest,powerfninverse,nfn
%hhat2 is redundant with hhat
%disp(['lnest8: nltype=' nltype])
%flags


if (nargin<13) crossvalopts.flag=0; end
if (nargin<12) derivflag=0; end
if (nargin<11) iter=1; end   %L interation
if (nargin<10) || ~exist('hlen2')
    hlen2=0;
    type='causal';
else
    type='2-sided';
end
if (nargin<9) nltype='powns'; end
if (nargin<8) detrendtype='constant'; end
if (nargin<7) plots=0; end
if (nargin<6) estmethod=1; end %full regression method
if (nargin<5)
    range=[1 size(input,1)];
elseif range==0
    range=[1 size(input,1)];
end

if (nargin<4) svdlength=0; end

if isvector(input)   %reformat input and output into column vectors
    if size(input,1)==1
        input=input';
        output=output';
    end
end

N=size(input,2);  %number of trials

    %trim input and output
    inp=input(range(1):range(2),:);
    outp=output(range(1):range(2),:);

    meancenter=1;
    %mean center
    if meancenter==1 
        inp_offset=mean(inp);
        outp_offset=mean(outp);
        if strcmp(detrendtype,'none') 
        else
            inp=detrend(inp,detrendtype);
            outp=detrend(outp,detrendtype);
        end
    else
        inp_offset=0;
        outp_offset=0;
    end

    %take derivative of output before running estimation
    if derivflag
        outp_orig=outp;
        outp=deriv(outp);
    end

%% LEST section

nltype

    [hhat hhat_causal vafout1 vafout1_causal interm interm_causal hgain]=kernelest(inp,outp,hlen,svdlength,estmethod,hlen2);

%% Now do NEST

% figure;plot(interm,'r');hold on;plot(outp,'b');
% figure;plot(hhat);
% figure;plot(interm(hlen:end),outp(hlen:end),'k.');
% figure;

    %N curvefit
    [npm shift]=nest(interm(hlen:end-hlen2),outp(hlen:end-hlen2),nltype);
    
    [npm_causal shift_causal]=nest(interm_causal,outp,nltype);

%% Compute VAFs

    %VAF of full L-N system
    [simtrace,vaf1,vaf2,interm]=lneval2(inp,hhat*hgain,npm,nltype,outp,0,0,hlen2,0.5);

   
    
%% ITERATION METHOD
%re-estimate linear kernel, only for 2-sided filter now.
%
    hhat_pass1=hhat;
    npm_iter(1,:)=npm;
    hhat_iter(:,1)=hhat;
    hgain_iter(1)=hgain;

    if iter>1 && crossvalopts.flag==0
        for i=2:iter

            inverse_interm=powerfninverse(outp,npm_iter(i-1,:));

            [hhat_iter(:,i), ~,~,~, interm2, ~, hgain_iter(i)]=kernelest(inp,inverse_interm,hlen,svdlength,estmethod,hlen2);
            [npm_iter(i,:) shift_iter(i)]=nest(interm2(hlen:end-hlen2),outp(hlen:end-hlen2),nltype);

        end

        if plots==1
            figure('Position',[0 0 800 800]);
            subplot(2,2,1);
            plot(hhat_iter);

            subplot(2,2,2);
            plot(interm2-shift,outp,'.r');
            hold on; 
            plot(interm2-shift,nfn(npm_iter(i+1,:),interm2-shift,nltype),'k');
            plot(inverse_interm-shift,outp,'.g');
        end


        
        [notused,vaf2L,vaf2,notused3]=lneval2(inp,hhat*hgain,npm,nltype,outp,0,0,hlen2,0.5);
        [notused,vaf2Li,vaf2i,notused3]=lneval2(inp,hhat_iter(:,i)*hgain_iter(i),npm_iter(i,:),nltype,outp,0,0,hlen2,0.5);
        
 
        if vaf2i>vaf2
            disp(['lnest9: vaf improved from ' num2str(vaf2) ' to ' num2str(vaf2i)]);
            hhat=hhat_iter(:,i);
            hgain=hgain_iter(i);
            shift=shift_iter(i);
            npm=npm_iter(i,:);

        else
            disp(['lnest9: vaf fell from ' num2str(vaf2) ' to ' num2str(vaf2i) '. using initial estimate.']);

        end
    end


%test plotting for iter method
%figure;
%plot(hhat2,'g');hold on;plot(hhat,'r');legend({'re-estimated','orig estimate'});

%% calc offset and gain

    offset=nfn(npm,0,nltype);
    gain=nfn(npm,.1,nltype); %set gain to polynomial value at 1;

%% plots and struct output
    if plots
        figure;
        subplot(4,1,1);
        plot(inp,'b');

        subplot(4,1,2);
        plot(outp,'r');
        hold on;
        %plot(interm(1:length(outp)),'b');
        plot(nfn(npm,interm(1:length(outp)),nltype),'g');  

        subplot(4,1,3);
        plot(-hlen2+1:hlen,hhat,'.-');
        hline(0);
        vline(0);

        subplot(4,1,4);
        scatter(interm(1:length(outp)),outp,8,1:length(outp),'filled');
        hold on;
        ndomain=linspace(min(interm),max(interm));
        plot(ndomain,nfn(npm,ndomain,nltype),'LineWidth',3,'Color','r');
        %legend([num2str(vafout) ; num2str(vaf2)]);
    end

    lnstruct.nw=size(hhat,2);
    lnstruct.h_est=hhat;
    lnstruct.hgain=hgain;
    lnstruct.h_len=hlen;
    lnstruct.h_len2=hlen2;
    lnstruct.h_iter=hhat_iter;
    lnstruct.h_est_causal=hhat(hlen2+1:end);
    lnstruct.h_est_pass1=hhat_pass1;
    lnstruct.npm=npm;
    lnstruct.nltype=nltype;
    lnstruct.gain=gain;
    lnstruct.offset=offset;
    lnstruct.shift=shift;
    lnstruct.vaf1=vaf1;
    lnstruct.vaf2=vaf2;
    lnstruct.interm=interm;
    lnstruct.intermrng=[min(interm) max(interm)];
    lnstruct.simtrace=simtrace;
    lnstruct.inp_offset=inp_offset;
    lnstruct.outp_offset=outp_offset;
    lnstruct.trueout=outp;
    lnstruct.inp=inp;
    lnstruct.iv=(-hlen2:hlen-1);

end

function [hhat, hhat_causal, vafL, vafL_causal, simtrace, simtrace_causal, hgain]=kernelest(inp,outp,hlen,svdlength,estmethod,hlen2)

    if ~exist('hlen2')
        hlen2=0;
        type='causal';
    else
        type='2-sided';
    end
    
    timelen=size(inp,1);
    N=size(inp,2);  %number of trials

    if estmethod==0 %cross correlation method

        %compute autocorrelation of input
        [inp_corr,ac_lags]=xcorr(inp,hlen,'biased');

        %compute Hessian
        Hess=zeros(hlen);
        for i=1:hlen
            for j=1:i
                Hess(i,j)=inp_corr(hlen+1+i-j);
                Hess(j,i)=Hess(i,j);
            end
        end

        %SVD smoothing
        if svdlength==0
            Hinv=inv(Hess);    
        else
           %{
            [U,S,VT]=svd(Hess);
            Hinv=VT*inv(S)*U'; 
            %}
            [m,n]=size(Hess);
            [U,S,V]=svd(Hess);
            %r=rank(S);
            r=svdlength;
            SR=S(1:r,1:r);
            SRc=[SR^-1 zeros(r,m-r);zeros(n-r,r) zeros(n-r,m-r)];
            Hinv=V*SRc*U.';
        end

        %cross correlation
        [xc,xc_lags]=xcorr(outp,inp,hlen,'biased');

        %estimate h
        hhat=Hinv*xc((hlen+1):(2*hlen));

    elseif strcmp(type,'2-sided')

       % use direct inversion
       Uneg=zeros(timelen,hlen);   %causal part of U  
       Unegn=zeros(timelen,hlen,N);

       for n=1:N
           Uneg(:,1)=Uneg(:,1)+inp(:,n); %first column is inp
           Unegn(:,1,n)=inp(:,n);

           for i=2:hlen    

               Uneg(:,i)=Uneg(:,i)+circshift(inp(:,n),i-1);  %delayed copies of input
               Unegn(:,i,n)=circshift(inp(:,n),i-1);

%                  Uneg(1:(i-1),i)=0;
%                  Unegn(1:(i-1),i,n)=0;
               
               %extrapolate input
                Uneg(1:(i-1),i)=inp(1,n);
                Unegn(1:(i-1),i,n)=inp(1,n);
               
           end

       end

       Upos=zeros(timelen,hlen2);   %acausal part of U
       Uposn=zeros(timelen,hlen2,N);

       for n=1:N
           for i=1:hlen2    
               thiscol=round(hlen2-i+1);  %advanced copies of input
               Upos(:,thiscol)=Upos(:,thiscol)+circshift(inp(:,n),-i);
               Upos((end-i+1):end,thiscol)=0;

               Uposn(:,thiscol,n)=circshift(inp(:,n),-i);
               Uposn((end-i+1):end,thiscol,n)=0;

           end
       end

       U=[Upos , Uneg];  %concatenate Upos and Uneg

       UTn=zeros(hlen+hlen2,timelen,N);
       for n=1:N
          UTn(:,:,n)=[Uposn(:,:,n) , Unegn(:,:,n)]';
       end
       Hess=U'*U;

       [m,n]=size(Hess);
       [VV,S,V]=svd(Hess);
       %r=rank(S);
       if svdlength==0
           r=size(S,1);
       else 
            r=svdlength;
       end
       SR=S(1:r,1:r);
       SRc=[SR^-1 zeros(r,m-r);zeros(n-r,r) zeros(n-r,m-r)];
       Hinv=V*SRc*VV.';

       UTz=zeros(hlen+hlen2,1);
       for n=1:N
          UTz=UTz+UTn(:,:,n)*outp(:,n);
       end

       hhat=Hinv*UTz;

    else  %1-sided causal
       % use direct inversion
       % does not currently accept array inputs, only vector

       %con
       U=zeros(timelen,hlen);
       U(:,1)=inp;
       %
       for i=2:hlen
           U(:,i)=circshift(inp,i-1);
           U(1:(i-1),i)=0;
       end

       hhat=inv(U'*U)*U'*outp;

    end

    %swap in cross validation filter if flag is set
    if isfield('crossvalopts','flag')
        if crossvalopts.flag==1
            hhat=crossvalopts.hhat;
            hlen=length(hhat)-hlen2;
        end
    end

    [hhat hgain]=normalize(hhat,0); %area-normalize

    hhat_causal=hhat(end-hlen+1:end);

    %VAF of just L stage with gain adjustment
    [simtrace,vafL,~,~]=lneval2(inp,hhat*hgain,0,[],outp,0,1,hlen2,0.5);
    [simtrace_causal,vafL_causal,~,~]=lneval2(inp,hhat_causal*hgain,0,[],outp,0,0,0,0.5);
  
end