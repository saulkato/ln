function f=getfolders(folder)

d=dir(folder);
k=1;
for i=1:length(d)
    if d(i).isdir && ~(strcmp(d(i).name,'.') ||  strcmp(d(i).name,'..'))
        
        f{k}=d(i).name;
        k=k+1;
        
    end
    
end