function [L3,L4,Re] = segment3n4(ima)
Re=zeros(size(ima(:,:,1)));
for ii=3:4
L=NormArray(ima(:,:,ii));
L=bwlabel(L>.5);
if(ii==4)
   % rely on segmentation in ii=3 and expect only areas in the vecinity of one of the two ends
  [a,b]=find(L3==1);
  P=polyfit(b,a,1);
%   figure,imagesc(L3);hold on;plot([1:size(L3,2)],P(1)*[1:size(L3,2)]+P(2))
[allx,ally] = puntsEnmig([1 size(L3,2)],round(P(1)*[1 size(L3,2)]+P(2)));
out=1;
vals=zeros(length(allx),1);
for jj=1:length(allx)
    try
    gg=L3(ally(jj),allx(jj));
    vals(jj)=gg;
    end
end
mask=zeros(size(L3));
aux=find(vals==1,1,'first');
mask(ally(aux),allx(aux))=1;
aux=find(vals==1,1,'last');
mask(ally(aux),allx(aux))=1;
mask=imdilate(mask,getnhood(strel('disk',10,0)));
L=L.*mask;
end
A=cell2mat(struct2cell(regionprops(L,'area')));
in=find(A==max(A));
L(L~=in)=0;
L(L==in)=1;

Re(L==1)=ii-2;
if(ii==3)
    L3=L;
else
    L4=L;
end
end

% Re=Re-imerode(Re,ones(3,3));

end