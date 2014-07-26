function stafolder(folder)
%Stimulus Triggered Analysis on a folder
%

if nargin<1
    folder=pwd;
end

cd(folder);
[starttimes iorstruct]=checksyncfolder(folder,1,0);

metadata=load('meta.mat');

plotflag=1;
thresh=.1;
pulselength=metadata.mseqpulselength;
sta_onetrace(iorstruct.meanqinp,iorstruct.meanoutp,thresh,plotflag,iorstruct.dt,pulselength,[],iorstruct.filename,metadata.startframe,metadata.endframe);

end