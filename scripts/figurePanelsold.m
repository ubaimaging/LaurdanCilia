function [] = figurePanels(ima,IT,refpos,PP,Nc,NF,pseudocm,fov,nom,RF,names)
sz=get(0,'screensize');
mg=[.01,.04,.01,.01];addist=1;ed=[-1:.01:2];
despColorbar=.08;
despHistogram=.14;
if(addist==1),yal='left';else,yal='right';end
     figure(67);clf;set(gcf,'position',sz);
         rawim = gray2rgb256(round(256*IT/NF),gray(257));
%    imwrite(rawim,[RF '/RAWim_' S(ii).name(1:end-4) '.jpg']); %
    
           subNM(2,1+ceil(sqrt(Nc)),1,mg);imagesc(rawim);axis image;set(gca,'xtick',[],'ytick',[]);
           xl=xlim;yl=ylim();
            h=text(xl(2),yl(2),nom);
        alp=0;
        set(h,'interpreter','none','fontname','calibri','fontsize',16,'verticalalignment','bot','horizontalalignment','right','color','w','margin',.1,'backgroundcolor',[1 1 1 alp]);
       
    subNM(2,1+ceil(sqrt(Nc)),ceil(sqrt(Nc))+2,mg);phasorImage(PP,fov,[],0,0,0);universalCircle(fov,'k',-1);
    

    for kk=1:Nc
        plot(refpos(kk,2),refpos(kk,1),'r.','markersize',22);
%         h=text(refpos(kk,2),refpos(kk,1),['{\it f_' num2str(kk) '}']);
h=text(refpos(kk,2),refpos(kk,1),[' ' names{kk}]);

        set(h,'hori','left','vert','midd','fontsize',11,'fontname','calibri');
        if(kk==2),set(h,'vert','bot');end
    end
    aux=refpos(1:Nc,:);%[~,ord]=sort(aux(:,2));aux=refpos(ord,:);
    aux(end+1,:)=aux(1,:);
    for kk=1:Nc,h=line([aux(kk,2) aux(kk+1,2)],[aux(kk,1) aux(kk+1,1)]);set(h,'color','r');end
    dadexp=cell(Nc,1);
    for kk=1:Nc
        imk=ima(:,:,kk);
        dadexp{kk}=imk;
       imk(imk<0)=0;imk(imk>1)=1;
%         fac=IT;fac=fac/max(max(fac));
        fac=1;
        imk=imk.*fac;
        imk(1,1)=1;imk(1,2)=0;  

        cmap=superjet(257,pseudocm{kk});% declare the colormap
        cmapim=permute(cmap(end:-1:1,:),[1,3,2]);% make an image for the colorbar
        imk=gray2rgb256(round(255*imk),cmap); 
%         imwrite(imk,[RF '/COMP' num2str(kk) '_' S(ii).name(1:end-4) '.jpg']);
        
        
        ex=0;if(kk>ceil(sqrt(Nc))),ex=1;end
        ax=subNM(2,1+ceil(sqrt(Nc)),1+kk+ex,mg);po=ax.Position;
        set(ax,'position',[po(1)-.05 po(2:4)]);
        imagesc(imk,[0 1]);axis image;set(gca,'xtick',[],'ytick',[]);
        xl=xlim;yl=ylim();
%         h=text(xl(2),yl(2),[' {\it f_' num2str(kk) '} ']);
        h=text(xl(2),yl(2),names{kk});
        alp=0;
        set(h,'fontname','calibri','fontsize',16,'verticalalignment','bot','horizontalalignment','right','color','w','margin',.1,'backgroundcolor',[1 1 1 alp]);
        po=ax.Position;
        realpo3=po(3)*sz(4)*3/sz(3)/2;
        axes('position',[po(1)+realpo3*(1+despColorbar) po(2) realpo3*.05 po(4)]);
        imagesc(cmapim);
        set(gca,'xtick',[],'ytick',[1 256],'yticklabel',[ ],'YAxisLocation',yal,'fontname','verdana','fontsize',14,'box','on');
        if(addist==1)
            axes('position',[po(1)+realpo3*(1+despHistogram) po(2) realpo3*.15 po(4)]);
            
            aux=ima(:,:,kk);aux(aux<0)=[];aux(aux>1)=[];
            distrib=histc(aux(:),ed)';
            
            distrib=medfilt1(distrib,3);distrib=distrib/(max(distrib)*1.1);
            
            h=patch([zeros(1,length(distrib)),distrib(end:-1:1)],[ed,ed(end:-1:1)],'k');
            set(h,'edgecolor','none','facealpha',.7);xlim([0 1]);ylim([0,1]);
            niceplot();set(gca,'yticklabel',[],'xticklabel',[]);
        end
    end
    
    saveas(67,[RF '/' nom '_COMPS' num2str(Nc) '.png']);%close(67);
 
    
end