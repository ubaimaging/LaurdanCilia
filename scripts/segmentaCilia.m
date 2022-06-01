function [mk]=segmentaCilia(I,gg)
% I1=I(:,:,:,1);% cilia
% I2=I(:,:,:,2);% base
% gg=3;
filtr=[3,3];
mk=zeros(size(I,1),size(I,2));
prob=mk;prob(round(size(I,1)/2),round(size(I,1))/2)=1;prob=bwdist(prob);
prob=max(prob(:))-prob;prob=NormArray(prob);
for ii=1:2
    aux=medfilt2(I(:,:,:,ii),filtr);
    switch gg
        case 1
            L = watershedSize(aux,200,1200,0);
        case 2
            aux=NormArray(aux.*prob);
            L = watershedSize(aux,200,1200,0);
        case 3
            aux=double(aux>graythresh(aux));
            L=bwlabel(aux);
        case 4
            aux=NormArray(aux.*prob);
            aux=double(aux>graythresh(aux));
            L=bwlabel(aux);
        case 5
            aux=NormArray(aux.*prob);
            aux(aux<.5)=.5;
            aux=NormArray(aux);
 aux=double(aux>graythresh(aux));
            L=bwlabel(aux);            
    end
    
    
    A=cell2mat(struct2cell(regionprops(L,'area')));
    [~,q]=max(A);
    if(~isempty(q))
        L(L~=q)=0;
    end
    L(L>0)=1;
    mk(L==1)=ii;
end

end