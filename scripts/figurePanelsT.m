function [panels] = figurePanels(ima,IT,refpos,PP,NF,colcode,fov,nom,names,type,nor)
% NF is the normalising factor in units of IT
% type=0 is for lifetime, 1 for spectral
% nor=0 for fractional images, nor=1 for photon-weighted fractional and nor=2 for photon-weighted with shifted LUT
opt=0;% 0 depict color of maximum, 1 do weighted average color (for panel with combined colors)
if((nargin<11)||(isempty(nor))),nor=1;end% default normalise fractional images by intensity
axft=10;% axis fontsize
sz=get(0,'screensize');
mg=[.01,.04,.01,.01]*2;addist=1;ed=[-1:.01:2];
despColorbar=(mg(1)+mg(3))*.2;
despHistogram=(mg(1)+mg(3))-despColorbar;
Nc=size(refpos,1);
if(addist==1),yal='left';else,yal='right';end
figure(67);clf;set(gcf,'position',sz);
if(isempty(NF)),NF=max(IT(:));NF=quantile(IT(:),.999);end
if(NF>1),NF=round(NF);end
IT=double(IT)/NF;IT(IT>1)=1;

rawim = gray2rgb256(round(256*IT),gray(257));

ax=subNM(2,1+ceil(sqrt(Nc)),1,mg);imagesc(rawim);axis image;set(gca,'xtick',[],'ytick',[]);
xl=xlim;yl=ylim();

alp=.5;
if(~isempty(nom))
h=text(xl(2),yl(2),nom);
set(h,'interpreter','none','fontname','calibri','fontsize',16,'verticalalignment','bot','horizontalalignment','right','color','k','margin',.1,'backgroundcolor',[1 1 1 alp]);
end

po1=ax.InnerPosition;
%     po2=ax.OuterPosition;
%         realpo3=po(3)*sz(4)*3/sz(3)/2;

if((Nc==6)||(Nc==7)),po1(2)=po1(2)+.0122;po1(4)=po1(4)*.937;end
axes('position',[po1(1)+po1(3) po1(2) despColorbar po1(4)]);
cmap=gray(257);% declare the colormap
cmapim=permute(cmap(end:-1:1,:),[1,3,2]);% make an image for the colorbar
imagesc(cmapim);
% set(gca,'xtick',[],'ytick',[1 257],'yticklabel',[NF 0],'yaxislocation','right','fontname','calibri','fontsize',axft,'box','on');
  set(gca,'xtick',[],'ytick',[ ],'box','on');
  aux1=xlim;aux2=ylim;
    h=text(mean(aux1),aux2(1)+((1-.985)*diff(aux2)),num2str(NF));set(h,'vertical','bot','horizontal','center','fontname','calibri','fontsize',axft,'margin',.1);
 h=text(mean(aux1),aux2(2),'0');set(h,'vertical','top','horizontal','center','fontname','calibri','fontsize',axft,'margin',.1);


subNM(2,1+ceil(sqrt(Nc)),ceil(sqrt(Nc))+2,mg);
% phasorImage(PP,fov,[],0,0,0);

%    phasorImage(PP,fov,[],-2,linspace(.001,.05,3));
   phasorImage(PP,fov,[],-2,[10.^linspace(-3,-1.3,3)]);
if(type==1)
    universalCircle(fov,'k',-8.3);
else
    universalCircle(fov,'k');
    h=line(fov(1,:),[fov(2,1) fov(2,1)]);set(h,'color','k');
% h=line(fov(1,:),[fov(2,2) fov(2,2)]);set(h,'color','k');
% h=line([fov(1,1) fov(1,1)],fov(2,:));set(h,'color','k');
% h=line([fov(1,2) fov(1,2)],fov(2,:));set(h,'color','k');
end


for kk=1:Nc
    h=plot(refpos(kk,2),refpos(kk,1),'o');
    set(h,'MarkerFaceColor',superjet(1,colcode(kk)),'color','k');

    try
        h=text(refpos(kk,2),refpos(kk,1),[' ' names{kk}]);
    catch
        h=text(refpos(kk,2),refpos(kk,1),['{\it f_' num2str(kk) '}']);
    end
    set(h,'hori','left','vert','bot','fontsize',11,'fontname','calibri','BackgroundColor',[1 1 1 .5],'margin',.1);
    %         if(kk==2),set(h,'vert','bot');end
end
aux=refpos(1:Nc,:);%[~,ord]=sort(aux(:,2));aux=refpos(ord,:);
aux(end+1,:)=aux(1,:);
% for kk=1:Nc,h=line([aux(kk,2) aux(kk+1,2)],[aux(kk,1) aux(kk+1,1)]);set(h,'color','r');end
dadexp=cell(Nc,1);
for kk=1:Nc
    imk=ima(:,:,kk);
    dadexp{kk}=imk;
    imk(imk<0)=0;imk(imk>1)=1;

    if(addist==1)
        aux=ima(:,:,kk);
        if(nor>=1)
                res=100;
     if(max(IT(:))<10),IT=IT*res;end
            IT=round(res*IT/max(IT(:)));
            gg=ceil(IT(:));ggS=aux(:);
            Sx=zeros(sum(gg),1);
            cc=1;
            for ii=1:length(gg)
                Sx(cc:cc+gg(ii)-1)=ggS(ii)*ones(gg(ii),1);
                cc=cc+gg(ii);
            end
            aux=Sx;
        end
        if(nor==0)
            aux(aux<0)=[];aux(aux>1)=[];
            aux(isnan(aux))=[];
        end
        distrib=histc(aux(:),ed)';
        thisNF=quantile(aux(:),.95);% over unity
        distrib=medfilt1(distrib,3);distrib=distrib/(max(distrib)*1.1);
    end
    if(nor==2)
        fac=IT;fac=1*fac/round(max(fac(:))*thisNF);
    end
    if(nor==1)
        fac=IT;fac=1*fac/max(fac(:));
    end
    if(nor==0)
        fac=1;
    end
    imk=imk.*fac;
    imk(1,1)=1;imk(1,2)=0;imk(imk>1)=1;

    cmap=superjet(385,['k' colcode(kk) 'w']);% declare the colormap
    cmap=cmap(1:257,:);
    cmapim=permute(cmap(end:-1:1,:),[1,3,2]);% make an image for the colorbar
    if(nor==2)
cmapim=cat(1,cmapim(1,:,:).*ones(round(size(cmapim,1)/thisNF)-size(cmapim,1),1,3),cmapim);
    end
    imk=gray2rgb256(round(255*imk),cmap);
    %         imwrite(imk,[RF '/COMP' num2str(kk) '_' S(ii).name(1:end-4) '.jpg']);


    ex=0;if(kk>ceil(sqrt(Nc))),ex=1;end
    ax=subNM(2,1+ceil(sqrt(Nc)),1+kk+ex,mg);po=ax.Position;
    set(ax,'position',[po(1)-mg(1) po(2:4)]);
    imagesc(imk,[0 1]);axis image;set(gca,'xtick',[],'ytick',[]);
    xl=xlim;yl=ylim();
    try
        h=text(xl(2),yl(2),names{kk});
    catch
        h=text(xl(2),yl(2),[' {\it f_' num2str(kk) '} ']);
    end
    alp=.5;
    set(h,'fontname','calibri','fontsize',16,'verticalalignment','bot','horizontalalignment','right','color','k','margin',.1,'backgroundcolor',[1 1 1 alp]);
    po1=ax.InnerPosition;
    %     po2=ax.OuterPosition;
    %         realpo3=po(3)*sz(4)*3/sz(3)/2;

    if((Nc==6)||(Nc==7)),po1(2)=po1(2)+.0122;po1(4)=po1(4)*.937;end
    axes('position',[po1(1)+po1(3) po1(2) despColorbar po1(4)]);
    imagesc(cmapim);
    set(gca,'xtick',[],'ytick',[1 256],'yticklabel',[ ],'YAxisLocation',yal,'fontname','calibri','fontsize',14,'box','on');
    if(addist==1)
        axes('position',[po1(1)+po1(3)+despColorbar po1(2) despHistogram po1(4)]);

        h=patch([zeros(1,length(distrib)),distrib(end:-1:1)],[ed,ed(end:-1:1)],'k');
        set(h,'edgecolor','none','facealpha',.7);xlim([0 1]);ylim([0,1]);
        niceplot();set(gca,'yticklabel',[],'xticklabel',[]);
        if(nor>=1),ytl=num2str(NF);end
        if(nor==0),ytl='1';end
        h=text(.5,0,'0');set(h,'vertical','top','horizontal','center','fontname','calibri','fontsize',axft,'margin',.1);
        h=text(.5,.985,ytl);set(h,'vertical','bot','horizontal','center','fontname','calibri','fontsize',axft,'margin',.1);
        if(nor==2)
            h=line([0 1],[thisNF thisNF]);set(h,'color','k')
        h=text(.5,thisNF-(1-.985),num2str(round(thisNF*NF)));set(h,'vertical','bot','horizontal','center','fontname','calibri','fontsize',axft,'margin',.1);
        end
    end
end
% compose single colored image
imcols=zeros(size(ima,1),size(ima,2),3);
if(opt==0)
    [~,q]=max(ima,[],3);
    for ii=1:size(ima,3)
        col=superjet(1,colcode(ii));
        in=find(q==ii);
        for cc=1:3
            ch=imcols(:,:,cc);
            ch(in)=col(cc);
            imcols(:,:,cc)=ch;
        end
    end
end
if(opt==1)
    if(1)
        ima2=ima;ima2(ima2<0)=0;ima2(ima2>1)=1;
    else
        ima2=ima-min(ima,[],3);ima2=ima2./max(ima2,[],3);
    end
    ima2=ima2./sum(ima2,3);
    for cc=1:3
        ch=zeros(size(ima,1),size(ima,2));
        for ii=1:size(ima,3)
            col=superjet(1,colcode(ii));
            ch=ch+ima2(:,:,ii)*col(cc);
        end
        imcols(:,:,cc)=ch;
    end
end
FAC=NormArray(IT);FAC=FAC/quantile(FAC(:),.99);FAC(FAC>1)=1;
imcols=imgaussfilt(imcols,1).*FAC; %figure,imagesc(imcols)
if(mod(Nc,2)~=0)% if extra space
    ax=subNM(2,1+ceil(sqrt(Nc)),2*(1+ceil(sqrt(Nc))),mg);
else
    figure(456);set(gcf,'color','w','position',[488   342   467   420]);set(gca,'position',[0 0 1 1]);
end
imagesc(imcols);axis image;set(gca,'xtick',[],'ytick',[]);
try
    h=text(xl(2),yl(2),'Dominating Component');
    set(h,'fontname','calibri','fontsize',16,'verticalalignment','bot','horizontalalignment','right','color','k','margin',.1,'backgroundcolor',[1 1 1 alp]);
end
if(mod(Nc,2)==0)% if extra space
    [Ff,~]=getframe(456);
    close(456);
end

if(Nc==4)||(Nc==3),set(gcf,'position',[10          10        1142         792]);end
if(Nc==5),set(gcf,'position',[ 47          96        1223         601]);end
% gg=get(gcf,'position');sz=get(0,'screensize');
% gg(1:2)=[1 1];gg(4)=gg(4)*sz(3)/gg(3);gg(3)=sz(3);
% fac1=gg(4)/gg(3);
fac2=(size(FAC,1))/(size(FAC,2));
if(fac2<.68),set(gcf,'position',[1           1        1536         520]);end
if(fac2>1.4),set(gcf,'position',[1           1        1536         864]);end
% gg(4)=round(gg(4)*fac2);
% set(gcf,'position',gg);
    %     saveas(67,[RF '/' nom '_' num2str(Nc) 'c.png']);%close(67);
[F,~]=getframe(gcf);
panels=F;
if(mod(Nc,2)==0)% if extra space
    mm=ceil(size(panels,1)*.52);
    Ff=imresize(Ff,[round(size(panels,1)*.4) size(Ff,2)*round(size(panels,1)*.4)/size(Ff,1)]);
    Ff=cat(1,255*ones(mm,size(Ff,2),3),Ff);
    Ff=cat(1,Ff,255*ones(size(panels,1)-size(Ff,1),size(Ff,2),3));
    panels=cat(2,panels,Ff);
end

end