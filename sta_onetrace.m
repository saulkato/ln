function stastruct=sta_onetrace(itrace,otrace,thresh,plotflag,dt,pulselength,saveCopyDir,filename,startframe,endframe)
%stimulus  triggered averaging
%one trial

%uses mseq pulselength from meta.mat *important to have filled out*


if nargin<6
    saveCopyDir='/Users/skato/Desktop/Dropbox/MATLAB/Chron/Output/Yifan';
end
    
if nargin<5
       mdt=2;
else
    mdt=pulselength;
end

    
if nargin<4
    plotflag=1;
end

if nargin<3
    thresh=0.05;
end


halfwindow=200;  %size of 1/2 time window in frames (from -window to +window)

time_seg=((-halfwindow):(halfwindow-1))*dt;


%find transitions
[numswings,numupswings,numdownswings,allswings,upswings,downswings]=swings(itrace,thresh);

as=find(allswings); 
us=find(upswings); 
ds=find(downswings);

%find all pulse onsets

currentstate=zeros(size(itrace));
upstates=upswings;
downstates=downswings;

if as(1)==us(1)  %check first swing is up

    
    for i=1:(numupswings-1)
        upstates( us(i): round(mdt/dt):(ds(i)-round(mdt/dt/2)) )=1;
        downstates( ds(i): round(mdt/dt):(us(i+1)-round(mdt/dt/2)) )=1;

        
    end

else %first swing is down
    
    for i=1:(numdownswings-1)
        downstates( ds(i): round(mdt/dt):(us(i)-round(mdt/dt/2)) )=1;
        upstates( us(i): round(mdt/dt):(ds(i+1)-round(mdt/dt/2)) )=1;
    end
    
end

numupstates=sum(upstates);
numdownstates=sum(downstates);

inp_up_seg=zeros(2*halfwindow,numupstates);
outp_up_seg=zeros(2*halfwindow,numupstates);

inp_down_seg=zeros(2*halfwindow,numdownstates);
outp_down_seg=zeros(2*halfwindow,numdownstates);

j=1;k=1;

for i=halfwindow:(length(otrace)-halfwindow+1)
    if upstates(i)==1     

        if j<=size(inp_up_seg,2)
            inp_up_seg(:,j)=itrace((i-halfwindow):(i+halfwindow-1));
            outp_up_seg(:,j)=otrace((i-halfwindow):(i+halfwindow-1))-otrace(i);

            %subtract value at t=now
            outp_down_seg(:,k)=outp_down_seg(:,k)-outp_down_seg(halfwindow,k);        
        end
        j=j+1;
    end
    if downstates(i)==1     

        if k<=size(outp_down_seg,2) && k<=size(outp_up_seg,2)

            inp_down_seg(:,k)=itrace(1+i-halfwindow:i+halfwindow);
            outp_down_seg(:,k)=otrace(1+i-halfwindow:i+halfwindow)-otrace(i);

             %subtract value at t=now
            outp_down_seg(:,k)=outp_down_seg(:,k)-outp_up_seg(halfwindow,k);

        end
        k=k+1;
    end

end

%pick seg subset
inp_up_mean=mean(inp_up_seg,2);
outp_up_mean=mean(outp_up_seg,2);
inp_down_mean=mean(inp_down_seg,2);
outp_down_mean=mean(outp_down_seg,2);

outp_mean=mean(otrace(startframe:endframe));

%create synthetic output
outp_synth=conv(outp_up_mean(halfwindow:end),upstates)+conv(outp_down_mean(halfwindow:end),downstates);
outp_synth=outp_synth-mean(outp_synth(startframe:endframe))+outp_mean;


vafout=vaf(otrace(startframe:endframe)-outp_mean,outp_synth(startframe:endframe)-outp_mean);

%PLOTS
nr=4;
if plotflag==1

    %traces and triggers and synthetic output
    figure('Position',[0 0 1200 800]);
    subplot(nr,2,1:2);
    plot(itrace/2-1/2,'r');
    hold on;
    plot(otrace,'b');

%     for i=1:length(us)
%         vline(us(i));
%     end
%         
%     for i=1:length(ds)
%         vline(ds(i),'k');
%     end

    plot(-1.5+upstates/2,'c.-','MarkerSize',8);
    plot(-2.5+downstates/2,'m.-','MarkerSize',8);
    
    plot(outp_synth,'g');
    
    legend({'input','output','on-triggers','off-triggers','synthetic outp'},'Location','SouthEast');
    textur(['vaf=' num2str(vafout,3) '%']);
    vline(startframe);
    vline(endframe);
    subplot(nr,2,3);

    hold on;
%         for i=1:numupstates
%             plot(time_seg',inp_up_seg(:,i),'Color',[1 .6 .6]);
%         end
    plot(time_seg',inp_up_mean,'LineWidth',1.5,'Color','r');
        %hline(50);
    xlabel('time (s)');
    title('ON-triggered averaged input');
    set(gca,'XMinorTick','on');


    subplot(nr,2,5);
%         for i=1:numdownstates
%             plot(time_seg',inp_down_seg(:,i),'Color',[1 .6 .6]);
%             hold on;
%         end
    plot(time_seg',inp_down_mean,'LineWidth',1.5,'Color','r');
    %hline(50);
    xlabel('time (s)');
    title('OFF-triggered averaged input');
    set(gca,'XMinorTick','on');


    %output traces
    subplot(nr,2,4);

    hold on;
    
%     for i=1:numupstates
%         %normfactor=max([abs(max(outp_up_seg(:,i))) abs(min(outp_up_seg(:,i)))]);
%         normfactor=1;
%         plot(time_seg',(outp_up_seg(:,i))/normfactor,'Color',[.6 .6 1]);
%     end

    
    plot(time_seg',outp_up_mean,'LineWidth',1.5,'Color','b');
    hline(0);
    xlabel('time (s)');
    ylim([min([outp_up_mean; outp_down_mean]) max([outp_up_mean; outp_down_mean]) ]);
    title('ON-triggered averaged output');
    %xlim([-15 15]);
    set(gca,'XMinorTick','on');

    subplot(nr,2,6);
    hold on;
%     for i=1:numdownstates
%         %normfactor=max([abs(max(outp_up_seg(:,i))) abs(min(outp_up_seg(:,i)))]);
%         normfactor=1;
%         plot(time_seg',outp_down_seg(:,i)/normfactor,'Color',[.6 .6 1]);
%     end

    plot(time_seg',outp_down_mean,'LineWidth',1.5,'Color','b');
    hline(0);
    xlabel('time (s)');
    %ylim([-50 50]);
    title('OFF-triggered averaged output');
    ylim([min([outp_up_mean; outp_down_mean]) max([outp_up_mean; outp_down_mean]) ]);
    %xlim([-15 15]);
    set(gca,'XMinorTick','on');


    subplot(nr,2,7)
    hold on;
    plot(time_seg',(inp_up_mean-inp_down_mean)/2,'LineWidth',1.5,'Color','r');
    %hline(0);
    xlabel('time (s)');
    title('ON and OFF merged input');
    ylim([min([inp_up_mean; inp_down_mean]) max([inp_up_mean; inp_down_mean]) ]);
    set(gca,'XMinorTick','on');

    
    subplot(nr,2,8)
    hold on;
    plot(time_seg',(outp_up_mean-outp_down_mean)/2,'LineWidth',1.5,'Color','b');
    hline(0);
    xlabel('time (s)');
    title('ON and OFF merged output');
    ylim([min([outp_up_mean; outp_down_mean]) max([outp_up_mean; outp_down_mean]) ]);    
    set(gca,'XMinorTick','on');

        
end

mtit(strrep(filename{1}(1:end-4),'_','\_'));


save2pdf([pwd filesep 'stta-' filename{1}(1:end-4) '.pdf'],gcf,600);

if exist('saveCopyDir','var') && ~isempty(saveCopyDir)
    save2pdf([saveCopyDir filesep 'stta-' filename{1}(1:end-4) '.pdf'],gcf,600);
end

stastruct.outp_synth=outp_synth;


