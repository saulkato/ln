function [xpos,ypos,traces,time,filename,metadata,baseline,deletelist,area]=loadlogfolder(folder,bleachflag,timecorrectflag,column,pixcalc,fcutoff,cellflag,percentflag,range,structout)
% [xpos,ypos,traces,time,filename,metadata,baseline]=loadlogfolder(folder,bleachflag,timecorrectflag)
%
% leave folder argument empty to choose directory interactively
% leave outputs empty to create variables in workspace

% -folder (optional)
% -bleachflag (perform bleach correction if 1)
%
%Saul Kato
%created 042510
%modified 080310
%mod 092710 added time correction
%mod 100410 added metafile load
%mod 120203 added brightness cutoff
%mod 120221 added area write out
%mod 120301 added trace range

areacolumn=6;

if nargin<10 || isempty(structout)
    structout=0;
end

if nargin<9 || isempty(range)
    range=0;
end

if nargin<8 || isempty(percentflag)
    percentflag=0;
    percent=1;
end

if percentflag==1
    percent=100;
else
    percent=1;
end

if nargin<7 || isempty(cellflag)
    cellflag=0;
end

if nargin<6 || isempty(fcutoff)
    fcutoff=0;
end

if nargin<5 || isempty(pixcalc)
    pixcalc='auto';
end

if nargin<4 || isempty(column)
    column=3;  %5 is fluoresence
end

if nargin<3 || isempty(timecorrectflag)
    timecorrectflag=0;
end

if nargin<2 || isempty(bleachflag)
    bleachflag=0;
end

if nargin<1
     folder = uigetdir('Select directory of log files');
end

cd(folder);

fnames=dir('*.log')
fnamestxt=dir('*.txt');

if isempty(fnames) || size(fnames,1)==size(fnamestxt,1)
    disp('loadlogfolder: no log files in this folder. trying .txt files');
    fnames=dir('*.txt');
    pixcalc=5; 
    if isempty(fnames)
        disp('loadlogfolder: no .txt files found either.  Giving up.');
        xpos=[];ypos=[];traces=[];time=[];filename=[];metadata=[];baseline=[];deletelist=[];area=[];
        return;
    end
end

i=1;
    
while i<=size(fnames,1)
    if size(fnames,2)==0
        disp('loadlogfolder: NO non-dye LOG FILES FOUND. trying .txt files');
        fnames=dir('*.txt');
        pixcalc=5;   
        if isempty(fnames)
            disp('loadlogfolder: no .txt files found either.  Giving up.');
            xpos=[];ypos=[];traces=[];time=[];filename=[];metadata=[];baseline=[];deletelist=[];area=[];
            return;
        end       

    end
    if strcmp(fnames(i).name(end-5:end-4),'ye') || ...
       strcmp(fnames(i).name(end-5:end-4),'YE') % ||  ...
     %  strcmp(fnames(i).name(end-5:end-4),'SK')
        fnames(i)=[];
        i=i-1;
    end
    i=i+1;
end

if range~=0
    fnames=fnames(range);
end

n=size(fnames,1);

%load metafile if it exists
if ~isempty(dir([folder '/meta.mat']))
     metadata=load([folder '/meta.mat']);
else
     disp('loadlogfolder: no meta.mat file in this folder.');
     m=load(fnames(1).name);
     metadata.startframe=1;
     metadata.endframe=size(m,1)-300;
end

if cellflag==0
    for i=1:n
        disp(fnames(i).name);
        m=load(fnames(i).name);
        mm{i}=m;

        if pixcalc(1)==1 || (strcmp(pixcalc(1),'a'))  && size(m,2)>3 %total area, with backsub
            time(:,i)=m(:,2)/1000;
            column=3;
            baseline(i)=mean(m(10:40,column));
            traces(:,i)=percent*m(:,column)/baseline(i);  %set df/f baseline
            area(:,i)=m(:,areacolumn);
            
        elseif pixcalc==2  %no backsub

            time(:,i)=m(:,2)/1000;
            column=5; areacolumn=6;
            pix=m(:,column)./m(:,areacolumn);      
            baseline(i)=mean(pix(10:40));

            traces(:,i)=percent*pix/baseline(i);
            area(:,i)=m(:,areacolumn);
        
        elseif pixcalc==3  %avg w backsub

            time(:,i)=m(:,2)/1000;
            column=5; areacolumn=6; avgbackgroundcolumn=7;  %5/6-7
            for j=1:size(m,1)
                if m(j,areacolumn)
                    pix(j)=m(j,column)./m(j,areacolumn)-m(j,avgbackgroundcolumn);    
                else
                    pix(j)=0;
                end
            end
            baseline(i)=mean(pix(10:40));
            traces(:,i)=percent*pix/baseline(i);
            area(:,i)=m(:,areacolumn);
            
        elseif pixcalc==4  %manuel logs

            time(:,i)=m(:,2)/1000;
            column=3; areacolumn=6; avgbackgroundcolumn=7;  %3/6-7
            for j=1:size(m,1)
                if m(j,areacolumn)
                    pix(j)=m(j,column); %./m(j,areacolumn); %-m(j,avgbackgroundcolumn);    
                else
                    pix(j)=0;
                end
            end
            baseline(i)=mean(pix(10:40));
            traces(:,i)=pix; %percent*pix/baseline(i);
            area(:,i)=m(:,areacolumn);
            
        elseif pixcalc==5 
            if size(m,2)==7  %new yifan andrew .txt logs
                disp('whatup')
                dt=0.05;  %should make this a metadata reference
                time(:,i)=m(:,1)*dt;
                column=2;  avgbackgroundcolumn=4;
                baseline(i)=mean(m(10:40,column)-m(10:40,avgbackgroundcolumn));
                traces(:,i)=(m(:,column)-m(:,avgbackgroundcolumn))/baseline(i); 
                area(:,i)=zeros(size(m(:,1)));

            else   %andrew .txt logs (NeuronTrackingSingle.txt)
                disp('.txt log does not have a time stamp and there is no meta file.  Guessing dt=0.05.');
                dt=0.05;  %should make this a metadata reference
                time(:,i)=m(:,1)*dt;
                column=2; areacolumn=9; avgbackgroundcolumn=4;
                baseline(i)=mean(m(10:40,column));
                traces(:,i)=(m(:,column)-m(:,avgbackgroundcolumn))/baseline(i); 
                if size(m,2) >= areacolumn
                    area(:,i)=m(:,areacolumn);
                else
                    area(:,i)=zeros(size(m(:,1)));
                end
            end
            
        elseif (strcmp(pixcalc(1),'a') && size(m,2)==3) %new stkanalyze

            m(end,:)=[];
            time(:,i)=((0:length(m(:,1))-1)*0.1)';
            column=2;
            baseline(i)=mean(m(10:40,column));
            traces(:,i)=percent*m(:,column)/baseline(i);
            area(:,i)=100;
        end    
        
        if ~isfield(metadata,'dt')
            metadata.dt=round(100*(time(2,1)-time(1,1)))/100;  %guess exposure time from first time interval of first trace  
            dt=metadata.dt;
        end
        
        if size(m,2)>7
            xpos(:,i)=m(:,8);
            ypos(:,i)=m(:,9);
        else
            xpos=zeros(size(m(:,1)));
            ypos=zeros(size(m(:,1)));
        end
        filename(i)={fnames(i).name};
        
        if timecorrectflag==1
            [traces(:,i),notused]=timecorrect(traces(:,i),time,dt);
            [area(:,i),notused]=timecorrect(area(:,i),time,dt);
        end
        
        if bleachflag==1
            traces(:,i)=bleachcorrect(traces(:,i));
        end          
    end
    
    %remove traces that are too bright
    deletelist=[];
    if ~fcutoff==0
        for i=1:n  %write out deletelist for loadlogfolder
            if baseline(i)>fcutoff
                deletelist=[deletelist i];
            end
        end       
        
        i=1;
        while i<=n
            if baseline(i)>fcutoff      
                traces(:,i)=[];
                trialname{i}=[];
                time(:,i)=[];
                area(:,i)=[];
                mm{i}=[];
                if size(m,2)>7  
                    xpos(:,i)=[];
                    ypos(:,i)=[]; 
                end
                baseline(i)=[];
                i=i-1;
                n=n-1;
            end
            i=i+1;
        end
    end
    
    
else  %cell version
    
    for i=1:n
        m=load(fnames(i).name);
        time{i}=m(:,2)/1000;  %time vector
        area{i}=m(:,areacolumn);
        total_baseline{i}=mean(m(10:40,column));
        traces{i}=m(:,column)/total_baseline{i};  %set df/f baseline
        baseline{i}=mean(m(10:40,column)./m(10:40,areacolumn));                
        if size(m,2)>7
            xpos{i}=m(:,8);
            ypos{i}=m(:,9);
        else
            xpos=zeros(size(m(:,1)));;
            ypos=zeros(size(m(:,1)));;
        end

        filename{i}=fnames(i).name;
        
        if timecorrectflag==1
            [traces{i},notused]=timecorrect(traces{i},time,dt);
            [area{i},notused]=timecorrect(area{i},time,dt);
        end
        
        if bleachflag==1
            traces{i}=bleachcorrect(traces{i});
        end
    end
        
    %remove traces that are too bright
    if ~fcutoff==0
        for i=1:n
            if baseline{i}>fcutoff
                traces{i}=[];
                filename{i}=[];
                area{i}=[];
                time{i}=[];
                mm{i}=[];
                disp('blah')
                if size(m,2)>7  xpos{i}=[];ypos{i}==[]; end
                baseline{i}=[];
                i=i-1;
                n=n-1;
            end
        end
    end
    
end

if structout
    out.xpos=xpos;
    out.ypos=ypos;
    out.time=time;
    out.area=area;
    out.fluor=traces;
    out.metadata=metadata;
    if exist('baseline')
        out.baseline=baseline;
    end
    out.mm=mm;
    clear xpos;
    xpos=out;
end

if nargout==0
    
    assignin('base','xpos',xpos);
    assignin('base','ypos',ypos);
    assignin('base','time',time);
    assignin('base','area',area);
    assignin('base','fluor',traces);
    assignin('base','filename',filename);
    assignin('base','metadata',metadata);
    if exist('baseline')
        assignin('base','baseline',baseline);
    end
    assignin('base','mm',mm);

end
