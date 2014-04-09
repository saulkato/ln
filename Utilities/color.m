function colortriple=color(colornamestring,colorrange,mycolormap)
% colortriple=color(colornamestring,colorrange)
%
%RGB values for commonly used colors
%or a value from a colormap range

colortriple=[0 0 0];

if nargin==1
    
    if strcmp(colornamestring,'lb') || strcmp(colornamestring,'lightblue')
        colortriple=[0.7 0.7 1];
    elseif strcmp(colornamestring,'llb') || strcmp(colornamestring,'lightlightblue')
        colortriple=[0.9 0.9 1];
    elseif strcmp(colornamestring,'mlb') || strcmp(colornamestring,'mediumlightblue')
        colortriple=[0.5 0.5 1];
    elseif strcmp(colornamestring,'lr') || strcmp(colornamestring,'lightred')
        colortriple=[1 0.7 0.7];
    elseif strcmp(colornamestring,'lg') || strcmp(colornamestring,'lightgreen')
        colortriple=[0.7 1 0.7];
    elseif strcmp(colornamestring,'dg') || strcmp(colornamestring,'darkgreen')
        colortriple=[0 0.7 0];
    elseif strcmp(colornamestring,'gray') || strcmp(colornamestring,'grey')
        colortriple=[0.5 0.5 0.5];
    elseif strcmp(colornamestring,'r') || strcmp(colornamestring,'red')
        colortriple=[1 0 0];
    elseif strcmp(colornamestring,'g') || strcmp(colornamestring,'green')
        colortriple=[0 1 0];
    elseif strcmp(colornamestring,'k') || strcmp(colornamestring,'black')
        colortriple=[0 0 0];
    elseif strcmp(colornamestring,'b') || strcmp(colornamestring,'blue')
        colortriple=[0 0 1];
    end
    
elseif nargin==2
    
    %colormap('jet')
    cmap=jet(colorrange); 
    if colorrange==3  %overrwite bad n=3 color mapping for jet
        cmap(3,:)=[1 0 0];
    end
    colortriple=cmap(colornamestring,:);
    
elseif nargin==3  %'rg'
        cmap=[(1:(-1/(colorrange-1)):0)', (0:(1/(colorrange-1)):1)', zeros(colorrange,1)];
        colortriple=cmap(colornamestring,:);
else
    colortriple=[0.5 0.5 0.5];
end