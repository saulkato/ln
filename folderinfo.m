function folderinfo(folder,doubletimecheckflag)
%folderinfo
%report info about folder
%   -number of log files
%   -sizes of log files
%
% created 4/29/10
% Saul Kato

if nargin<2
    doubletimecheckflag=0;
end

if nargin<1
%     folder = uigetdir('Select directory of log files');
    folder=pwd;
end

cd(folder);
fnames=dir('*.log');
n=size(fnames,1);
fprintf('folder: %s\n',folder');
fprintf('number of log files: %d\n',n);
for i=1:n
    m=load(fnames(i).name);
    fprintf('%s: %d lines, %d columns\n',fnames(i).name,size(m,1),size(m,2));
    
    %find double time entries
    if doubletimecheckflag
        timevec=m(:,2);
        tdiff=timevec-circshift(timevec,1);
        entries=find(tdiff==0);
        if ~isempty(entries)
            disp(['     DOUBLE TIME entry at ' num2str(length(entries)) ' lines.']);
            disp(entries);
        end
    end
end

end
