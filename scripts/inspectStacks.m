function [] = inspectStacks(vol,DX,DZ,col,fil,tt)
% fil col is optional for display
sz=get(0,'screensize');
sz(2)=41;sz(4)=sz(4)-116;
vol=double(vol);


N=size(vol,3);
if(nargin==4)
   fil=ceil(N/col); 
end
if(nargin<4)
[fil,col]=squareDistrib2(N,sz(4)/sz(3));
end


imratio=size(vol,1)*fil/size(vol,2)/col;
figratio=sz(4)/sz(3);

if(figratio/imratio<1)
    sz(3)=sz(3)*figratio/imratio;
else
    sz(4)=sz(4)*imratio/figratio;
end

figure;set(gcf,'position',sz);
if(size(vol,4)==1)
rg=[min(min(min(vol))) max(max(max(vol)))];
if(rg(1)==rg(2)),rg(2)=rg(1)+1;end
else
    if(max(vol(:))>1)
        rg=[0,256];
    else
    rg=[0,1];
    end
end
mgsz=2;

mg=[mgsz/sz(4) mgsz/sz(3) mgsz/sz(4) mgsz/sz(3)];
axx=zeros(1,N);
for ii=1:N
    
    aux=subNM(fil,col,ii,mg);
    axx(ii)=aux;
Im=vol(:,:,ii,:);
if(size(vol,4)==3),Im=permute(Im,[1 2 4 3]);if(rg(2)==256),Im=uint8(Im);end;end
    if(((nargin>=2)&&(ii==N))&&(~isempty(DX)))
        Im=enxufaLlegenda3(Im,DX,rg(2),'se',0,0);
    end
    imagesc(Im,rg);
    axis image;
    set(gca,'xtick',[],'ytick',[]);
    txt=[' frame#' niceDigits(ii,2)];
    if(nargin>=3)&&(~isempty(DZ))
        txt=[' Stack#' niceDigits(ii,2) ', +' niceNums(DZ*(ii-1)) '{\mu}m'];
    end
    salta=0;
    if((nargin<6)||isempty(tt))
    h=text(1,1,txt);
    else
       if(numel(tt)==N)
          h=text(1,1,tt{ii});
       else
           salta=1;
       end
    end
    if(salta==0)
   set(h,'fontname','verdana','fontsize',14,'verticalalignment','top','horizontalalignment','left','color','w'); 
    end
set(gca,'box','on');    
end
linkaxes(axx);
set(gcf,'color','w');

colormap(superjet(255));




end
