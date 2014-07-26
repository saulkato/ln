function outvec=nfn(pm,invec,type,shift)
%outvec=nfn(pm,invec,type,shift)
%
%shift optional
%evaluate a nonlinearity function
%type: 'poly5','poly','sr','pow','powfix'
%
%uses polyval,softrectifierfn,powerfn
%
%Saul Kato
%101110

if nargin<4
    shift=min(invec); %-0.0001;
    %shift=0;
end

if strcmp(type,'poly5')
    outvec=polyval(pm,invec);
elseif strcmp(type,'poly')
    outvec=polyval(pm,invec);
elseif strcmp(type,'sr')
    outvec=softrectifierfn(invec,pm);
elseif strcmp(type,'pow')
    disp('pow is deprecated.'); 
    outvec=powerfn(invec-shift,pm);    
elseif strcmp(type,'pows')
    outvec=powerfn(invec,pm);
elseif strcmp(type,'powns')
    outvec=powerfn(invec,pm);
elseif strcmp(type,'powfix')
    outvec=powerfn(invec-shift,pm);
elseif strcmp(type,'hardrect')
    outvec=hardrectifierfn(invec,pm);
end