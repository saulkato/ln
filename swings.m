function [numswings,numupswings,numdownswings,dtrace,uptrace,downtrace]=swings(trace,thresh)
%[numswings,numupswings,numdownswings,swings,upswings,downswings]=SWINGS(trace,thresh)
%
% count number of swings in a signal trace that cross a middle band
% centered around the midpoint level of the signal
%
% a swing must entirely cross the band in order to be counted
%
% band is defined by thresh, should be between 0 and 1.
% thresh=0.01 means a tiny band
% thresh=0.9 means a large band
%
% numswings     number of total swings
% numupswings   number of upswings
% numdownswings number of downswings
% swings        occupancy {1,0} vector of swings (aligns with trace)
% upswings      occupancy {1,0} vector of upswings (aligns with trace)
% downswings    occupancy {1,0} vector of downswings (aligns with trace)
%
% calling swings with no arguments plots a demonstration
%
% Saul Kato
% created 4/22/10
% modified 9/12/20 to add up and down counting and to return occupancy vector of
% swing points

%demo
if nargin<1  %test trace
    trace=2*sin(.1*[0:100]).*sin(.5*[0:100]);
    figure;
    plot(trace);
end

if nargin<2
    thresh=0.2;
end

maxtrace=max(trace);
mintrace=min(trace);
midpoint=(maxtrace+mintrace)/2;
halfrange=(maxtrace-mintrace)/2;

if trace(1)>0
    rtrace(1)=1;
else
    rtrace(1)=-1;
end

for i=2:length(trace)

    if trace(i)>thresh*halfrange+midpoint
        rtrace(i)=1;
    elseif trace(i)<-thresh*halfrange+midpoint
        rtrace(i)=-1;
    else
        rtrace(i)=0;
    end
   
end

if nargin<1
    hold on
    plot(rtrace,'r.');
    hline(thresh*halfrange+midpoint,'g');
    hline(-thresh*halfrange+midpoint,'g');    
end

for i=2:length(rtrace)
    if rtrace(i)==0 
        rtrace(i)=rtrace(i-1);
    end 
end

if nargin<1
    hold on
    plot(rtrace,'c-');
end

dtrace=zeros(size(trace));
uptrace=zeros(size(trace));
downtrace=zeros(size(trace));

for i=2:length(trace)
    if rtrace(i)==rtrace(i-1)
        dtrace(i)=0;
        uptrace(i)=0;
        downtrace(i)=0;
    else
        dtrace(i)=1;  
        if rtrace(i)>rtrace(i-1)
            uptrace(i)=1;
            downtrace(i)=0;
        else
            downtrace(i)=1;
            uptrace(i)=0;
        end
    end
end

%demo
if nargin<1  
    hold on;
    plot(dtrace,'r');
    vline(find(uptrace),[.8 .2 .2]);
    vline(find(downtrace),[.2 .2 .8]);
    hline(midpoint);
end

numswings=sum(dtrace);
numupswings=sum(uptrace);
numdownswings=sum(downtrace);
