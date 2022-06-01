function [IPF] = watershedSize(IM0,minsize,maxsize,verb)
if(nargin<4),verb=0;end
steps=100;
IM0=double(IM0);
if(sum(round(IM0(:)))==sum((IM0(:))))
    rng=linspace(min(IM0(:)),max(IM0(:)),min([max(IM0(:))-min(IM0(:))+1,steps]));
else
    rng=linspace(min(IM0(:)),max(IM0(:)),steps);
end
[a0,b0,c0]=size(IM0);
IP=zeros(a0,b0,c0);
contador=0;a=0;b=0;c=0;mc=1;
if(verb==1),ensenya( ' -Selecting cells. 00%');end
tots=numel(rng);
for ii=1:tots
    if(verb==1),percentatge(ii,tots); end   
    L=bwlabeln(IM0>rng(end-ii+1));
    
    numcells=max(L(:));
    for jj=1:numcells
        tamany=numel(find(L==jj));
        if(tamany<maxsize)&&(tamany>minsize)
            valors=unique(IP(L==jj));valors(valors==0)=[];
            if(numel(valors)<=1)
                contador=contador+1;IP(L==jj)=contador;
            else
                clear a b c;aux=find(L==jj);
                [a,b,c]=ind2sub([a0,b0,c0],aux);
                for kk=1:numel(a)
                    try
                    valor=IP(a(kk),b(kk),c(kk));
                    catch
                        gg=0;
                    end
                    if(valor==0)
                        entorn=IP(max(a(kk)-mc,1):min(a(kk)+mc,a0),max(c(kk)-mc,1):min(c(kk)+mc,c0),max(c(kk)-mc,1):min(c(kk)+mc,c0));
                        vals=unique(entorn);vals(vals==0)=[];
                        if(numel(vals)==1)
                            IP(a(kk),b(kk),c(kk))=vals;
                        end
                    end
                end
            end
        end
    end
end

aux=unique(IP);
aux(aux==0)=[];
IPF=zeros(size(IP));
for ii=1:length(aux)
    IPF(IP==aux(ii))=ii;
end
end