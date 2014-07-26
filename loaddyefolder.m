function [dye,time,filename]=loaddyefolder(folder,deletelist)
%[dye,time,filename]=loaddyefolder(folder)
% leave folder argument empty to choose directory interactively
% leave outputs empty to create variables in workspace
%
%Saul Kato
%created 042510
%modified 050210, 101010
%020212 added deletelist to work with fcutoffs in loadlogfolder

if nargin<2
    deletelist=[];
end

if nargin<1
     folder = uigetdir('Select directory of log files');
end

cd(folder);

fnames=dir('*.log');

if size(fnames,1)==0
    disp('loaddyefolder.m warning: no dye traces found in the dir. trying .txt files');
    fnames=dir('*.txt');
    correctlength=5001;  %hack
    if isempty(fnames)
        disp('loaddyefolder.m: no .txt files found either.  quitting.');
        xpos=[];ypos=[];traces=[];time=[];filename=[];metadata=[];baseline=[];deletelist=[];area=[];
        return;
    end
else
    folder
    correctlength=getconsensusloglength(folder);
    i=1;
    while i<=size(fnames,1)
        if (~strcmp(fnames(i).name(end-5:end-4),'ye')) && ...
           (~strcmp(fnames(i).name(end-5:end-4),'YE'))
            fnames(i)=[];
            i=i-1;
        end
        i=i+1;
    end
end


n=size(fnames,1);

arrayversion=1; cellversion=0;

if arrayversion==1
    for i=1:n
        fnames(i).name;
        m=load(fnames(i).name);
        filename(i)={fnames(i).name};  
        m(correctlength:end,:)=[];
        
        if size(m,2)==2  %this is an imageJ dye log file         
            dye(:,i)=(m(:,2)-min(m(:,2)));
            dye(:,i)=dye(:,i)/max(dye(:,i));           
            time(:,i)=zeros(size(dye(:,i)));  %should eventually pull tiff timestamp data  
            
        elseif size(m,2)==3
            disp('loaddyefolder.m: three columns in dye file');
            time(:,i)=m(:,1);
            col=size(m,2)-1;
            m(:,col)=m(:,col)-min(m(:,col));
            m(:,col)=m(:,col)/max(m(:,col));  
            dye(:,i)=m(:,col);
            
        elseif size(m,2)==9  %this is an Andrew .txt file
            time(:,i)=m(:,1)/20;
            col=6;  %dye column see NeuronTrackingSingle.txt
            m(:,col)=m(:,col)-min(m(:,col));
            m(:,col)=m(:,col)/max(m(:,col));  
            dye(:,i)=m(:,col);
            
        else %this is a metamorph log file with 4 columns          
            time(:,i)=m(:,2);
            col=size(m,2)-1;
            m(:,col)=m(:,col)-min(m(:,col));
            m(:,col)=m(:,col)/max(m(:,col));  
            dye(:,i)=m(:,col);      
        end
    end
end

if cellversion==1
    for i=1:n
        fnames(i).name;
        m=load(fnames(i).name);
        filename{i}=fnames(i).name;
        m(end,:)=[];
        if size(m,2)==2  %this is an imageJ dye log file
            dye(:,i)=(m(:,2)-min(m(:,2)));
            dye(:,i)=dye(:,i)/max(dye(:,i));           
            time{i}=0;  %should eventually pull tiff timestamp data  
        else
            time{i}=m(:,2);
            col=size(m,2)-1;
            m(:,col)=m(:,col)-min(m(:,col));
            m(:,col)=m(:,col)/max(m(:,col));  
            dye{i}=m(:,col); 
        end
    end
    
end

%check for missing time column
if max(max(time))==0
    disp('loaddyefolder: time column is missing. creating synthetic time column.');
    dt=0.1;
    size(time)
    for i=1:n
        time(:,i)=(0:(size(dye,1)-1))*dt;
    end
end

%if deletelist exists delete those traces
dye(:,deletelist)=[];
time(:,deletelist)=[];
filename(deletelist)=[];

if nargout==0
    assignin('base','dye',dye);
    assignin('base','timedye',time);
    assignin('base','filenamedye',filename);
end
