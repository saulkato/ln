function vline(x,colr,pattern)
%vline(x,colr)
%draw vertical line at each value of x and color
%
%Saul Kato
%created 4/13/10
%modified 4/22/10 to support vector input, a list of vertical positions

if nargin<3 
     pattern='-';
end

if nargin<2
    colr=[0.5 0.5 0.5];
end

if nargin<1
    x=0;
end

yl=ylim;

for i=1:length(x)
line([x(i) x(i)],[yl(1) yl(2)],'Color',colr,'LineStyle',pattern);

end
