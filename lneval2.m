 function [model_outp,vaf1,vaf2,model_interm,sse2]=lneval2(inp,h,pm,nltype,outp,plots,centeroutputflag,h_len2,padvalue)
%[model_outp,vaf1,vaf2,model_interm,sse2]=lneval2(inp,h,pm,nltype,outp,plots,centeroutputflag,h_len2)
%run L-N model on input array to generate simulated output and measure
%against real output.
%
%use pm=0 to skip use of N stage
%
%LNEVAL2 2/14/2011 added two-sided kernel evaluation
%
%uses nfn.m
%Saul Kato

if (nargin<9)  padvalue=NaN; end
if (nargin<8)  h_len2=0; end  %0 indicates pure causal kernel
if (nargin<7)  centeroutputflag=1; end
if (nargin<6)  plots=0;  end
if (nargin<5)  outp=0;   end
if (nargin<4)  nltype='pow';  end


h_len=length(h);
h_len1=h_len-h_len2;

model_outp=zeros(size(inp));
model_interm=zeros(size(inp));

for i=1:size(inp,2)

   %figure;plot(interm);hold on; 
  
   if isnan(padvalue) 
       thispadvalue=inp(1,i); 
   else
       thispadvalue=padvalue; 
   end
   interm=conv([thispadvalue*ones(h_len1,1);inp(:,i);inp(end,i)*ones(h_len2,1)],h); %convolution on padded input
   interm=interm(h_len+1:h_len+size(inp,1));
  
   %plot(interm,'--r');
   model_interm(:,i)=interm(1:length(outp));  %truncate interm
   %model_interm(:,i)=interm;  
    
    if pm==0   %No N stage
        model_outp(:,i)=model_interm(:,i);
    else
        model_outp(:,i)=nfn(pm,model_interm(:,i),nltype);
    end
    
    if centeroutputflag==1
        outp_centered_trunc(:,i)=detrend(outp(h_len1:end-h_len2,i),'constant');
    else
        outp_centered_trunc(:,i)=outp(h_len1:end-h_len2,i);
    end
    
    vaf1(i)=vaf(outp_centered_trunc(:,i),model_interm(h_len1:end-h_len2,i));
    vaf2(i)=vaf(outp_centered_trunc(:,i),model_outp(h_len1:end-h_len2,i));
    %vaf2(i)=vaf(outp_centered(:,i),model_outp(h_len:end,i));
    sse2(i)=sse(outp_centered_trunc(:,i),model_outp(h_len1:end-h_len2,i));
end

%plots=0;

if plots
    figure;
    plot(inp,'Color',[0.8 0.8 0.8]);
    hold on;
    plot(outp_centered_trunc,'r');
    plot(model_outp(h_len:end,:),'g');
    plot(model_interm(h_len:end,:),'b');
    intitle('LNEVAL subroutine');
    legend({'input','real output','model outp','model interm. signal'},'FontSize',8);
end