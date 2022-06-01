function [nI] = gray2rgb256(I,map)

if(size(map,1)==255),map=[map;map(end,:)];end
if(size(I,3)>1)
    nI=zeros(size(I,1),size(I,2),size(I,3),3);
    for ii=1:size(I,3)

aux=gray2rgb256(I(:,:,ii),map);
nI(:,:,ii,:)=permute(aux,[1 2 4 3]);
    end
else
nI=zeros(size(I,1),size(I,2),3);
df=size(I,1)*size(I,2);
df2=2*df;
I=double(I);
vals=unique(I);
if(numel(vals)==1)
    nI(:,:,1)= map(1,1);
    nI(:,:,2)= map(1,2);
    nI(:,:,3)= map(1,3);
else
    
    %vv=unique(I);
    %rg=linspace(min(vv),max(vv),size(map,1));
    vals(isnan(vals))=[];%remove nans
    for ii=vals'
        
        qui=find(I==ii);
        
        nI(qui)=map(ii+1,1);
        nI(qui+df)=map(ii+1,2);
        nI(qui+df2)=map(ii+1,3);
        
        
    end
    
end
end
end