function [] = figure2comps0(ima,IT,laumap,nom,RF,mk)
if(nargin<6)||(isempty(mk)),mk=ones(size(IT));end
mkedge=mk-imerode(mk,ones(3,3));
mk(mk>1)=1;
norma=sum(ima(:,:,1:2),3);
L1=ima(:,:,1)./norma;% norelax
L2=ima(:,:,2)./norma;% relax
aux1=L1;aux1(mk==0)=[];
aux2=L2;aux2(mk==0)=[];
aux=cat(3,aux1,aux2);
%     figureHists(aux,2,[nom 'c'],RF);

figure(68);clf,set(gcf,'position',[91  392  1332  370]);
mg=[.01 .03 .01 .02];
%     fracim=2*(L2-.5);fracim(fracim<-1)=-1;fracim(fracim>1)=1;
fracim=L2;fracim(fracim<0)=0;fracim(fracim>1)=1;
Imbw=repmat(NormArray(IT),[1,1,3]);%Imbw=uint8(256*Imbw);
subNM(1,3,1,mg);imagesc(Imbw);set(gca,'xtick',[],'ytick',[]);
axis image;
%     subNM(1,3,2,mg);imagesc(fracim,[-1,1]);set(gca,'xtick',[],'ytick',[]);axis image;h=colorbar;
%     set(h,'ticks',[-1 0 1],'ticklabels',{'no-relajado','50/50','relajado'})

fracimcol=gray2rgb(NormArray(fracim),laumap);

subNM(1,3,2,mg);imagesc(fracimcol+100*mkedge,[0,1]);set(gca,'xtick',[],'ytick',[]);axis image;
set(gca,'xtick',[],'ytick',[]);axis image;h=colorbar;
%     set(h,'ticks',[-1 0 1],'ticklabels',{'no-relajado','50/50','relajado'})
ti=[0:.1:1];for ii=1:length(ti),tl{ii}=niceNums(ti(ii),1);end
set(h,'ticks',ti,'ticklabels',tl);
niceplot;colormap(laumap);
aux=IT;aux=aux/quantile(aux(:),.999);aux(aux>1)=1;
subNM(1,3,3,mg);
if(numel(unique(mk))==1)
    imagesc(aux.*fracimcol);
else
fracimcol=fracimcol.*mk+Imbw.*(1-mk);
imagesc(fracimcol,[0,1]);
end
set(gca,'xtick',[],'ytick',[]);
axis image;
niceplot;
saveas(68,[RF filesep nom '.jpg']);

%     clf,set(gcf,'position',[91.0000  298.6000  602.8000  463.4000]);alp=0;colormap(parula);
%     mg=[.01 .03 .01 .02];
%     subNM(2,2,1,mg);imagesc(L1.*mk,[0,1]);set(gca,'xtick',[],'ytick',[]);axis image;colorbar;niceplot
%     xl=xlim();yl=ylim();h=text(xl(2),yl(2),names{1});set(h,'fontname','calibri','fontsize',16,'verticalalignment','bot','horizontalalignment','right','color','w','margin',.1,'backgroundcolor',[1 1 1 alp]);
%     subNM(2,2,3,mg);imagesc(L2.*mk,[0,1]);set(gca,'xtick',[],'ytick',[]);axis image;colorbar;niceplot
%     xl=xlim();yl=ylim();h=text(xl(2),yl(2),names{2});set(h,'fontname','calibri','fontsize',16,'verticalalignment','bot','horizontalalignment','right','color','w','margin',.1,'backgroundcolor',[1 1 1 alp]);
%     subNM(2,2,2,mg);imagesc(Imbw.*L1.*mk);set(gca,'xtick',[],'ytick',[]);axis image
%     xl=xlim();yl=ylim();h=text(xl(2),yl(2),names{1});set(h,'fontname','calibri','fontsize',16,'verticalalignment','bot','horizontalalignment','right','color','w','margin',.1,'backgroundcolor',[1 1 1 alp]);
%     subNM(2,2,4,mg);imagesc(Imbw.*L2.*mk);set(gca,'xtick',[],'ytick',[]);axis image
%     xl=xlim();yl=ylim();h=text(xl(2),yl(2),names{2});set(h,'fontname','calibri','fontsize',16,'verticalalignment','bot','horizontalalignment','right','color','w','margin',.1,'backgroundcolor',[1 1 1 alp]);
%     saveas(68,[RF '/' nom 'b.png']);
end