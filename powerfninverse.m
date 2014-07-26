function out=powerfninverse(input,pm)
%
%inverse of POWERFN
%
%Saul Kato
%2/09/11
%
%if pm is omitted default values are
% pm(1)=p=2.3
% pm(2)=x0=0

if length(pm)<4
    pm(4)=0;
end

if nargin<2
    pm(1)=1;
    pm(2)=1;
    pm(3)=2.3;
end

if length(pm)<3
    p=2.3;
else
    p=pm(3);
end

out=real(((input-pm(2))/pm(1)).^(1/p))+pm(4);

end