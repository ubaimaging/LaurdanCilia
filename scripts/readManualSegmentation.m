function [mk] = readManualSegmentation(ruta,DX,I1,I2,rm,RF)

thFac=.75;% threshold is a factor between the lowest intensity on the cilia ridge and the lowest background around the cilia
I0=imread([ruta]);

if(size(I0,1)~=256)% lack of consistency in experiments, if images have different sizes they need to be equalised
I0=imresize(I0,[256 256],'nearest');
end
nom=ruta(find(ruta==filesep,1,'last')+1:end-4);
radi=ceil(500/(DX*1e9));% swap to nm, numerator is the radius in nm
temp1=getnhood(strel('disk',radi,0));% extended area around cilia to get bg
radi=ceil(250/(DX*1e9));% swap to nm, numerator is the radius in nm
temp2=getnhood(strel('disk',radi,0));% area around cilia
radi=ceil(300/(DX*1e9));% swap to nm, numerator is the radius in nm
temp3=getnhood(strel('disk',radi,0));% 

% 2sigma 68%
% 3sigma 95%
% 4sigma 99.7%
% sigmaFac=3; % full width at X max (insert fraction over unity at which to measure width)
mk=[];
aux=double(I0(:,:,1))-double(I0(:,:,2));
aux=aux';
aux=aux(rm+1:end-rm,rm+1:end-rm);
L2=aux;
L=L2*0;
L(L2>0)=1;% cilia
L(L2<0)=2;% base
fill=0;
if(numel(find(L==2))>0)
if(bweuler(L==2)==0),fill=1;end
end
% figure, imagesc(L);% figure, imagesc(I0)

% refcol1=[255 0 255];% magenta fr cilia
% refcol2=[0 255 255];% cian for base
% L=zeros(size(I0,1),size(I0,2));
% L=L(rm+1:end-rm,rm+1:end-rm);
% I0cut=I0;I0cut=I0cut(rm+1:end-rm,rm+1:end-rm,:);
% un=intersect(find(I0cut(:,:,1)==refcol1(1)),find(I0cut(:,:,2)==refcol1(2)));un=intersect(un,find(I0cut(:,:,3)==refcol1(3)));
% L(un)=1;
% un=intersect(find(I0cut(:,:,1)==refcol2(1)),find(I0cut(:,:,2)==refcol2(2)));un=intersect(un,find(I0cut(:,:,3)==refcol2(3)));
% L(un)=2;
% L=L';% figure,imagesc(L)
contin=1;
% if(max(L(:))>2),disp(['ACHTUNG: found more than 2 regions! ' ruta(find(ruta==filesep,1,'last'):end) ')']);contin=0;end
% if(max(L(:))<2),disp(['ACHTUNG: found less than 2 regions! ' ruta(find(ruta==filesep,1,'last'):end) ')']);contin=0;end
 if(max(L(:))>2),contin=0;ensenya(['Expecting two regions in png mask image, found ' num2str(max(L(:))) '!'],'r');end
 if(max(L(:))<2),contin=0;ensenya(['Expecting two regions in png mask image, found ' num2str(max(L(:))) '!'],'r');end
if(contin==1)
    La=cell2mat(struct2cell(regionprops(L,'area')));
    ind=[1,2];
   
    if(length(La)==1)
error('Expecting two regions in png mask image, only found one!');
    else
    if(La(1)<La(2)),ind=[2,1];end
    end
    
    mk0=zeros(size(L,1),size(L,2));
    tam=3;% radi de la finestra
    dst=(2/(DX*1e6));% distancei per traçar perfil dintensitat (introduir valor en microns, dst es en pix)
    pinta=4;
    if(pinta==3),figure(100);clf;set(gcf,'position',[489.0000  507.4000  560.0000  255.6000]);end
    for ii=1:2
        
        mk=mk0;
        
        mk(L==ind(ii))=1;
        %     figure,imagesc(I1);figure,imagesc(imgaussfilt(I1,3))
        if(ii==1),ref=I1;end
        if(ii==2)
            if(fill==0)
            mk12=imdilate(mk,temp3);
            else
             aux1=bwlabel(1-mk);
             aux1(aux1==aux1(1,1))=0;
             mk(aux1>0)=1;
             mk12=mk;
            end
        else
            I11=imgaussfilt(ref,.05/(DX*1e6));
            %     mkr=mk;
            I11=NormArray(I11);
            aux=I11;
            aux(mk==0)=0;
            %     figure,imagesc(aux);
            top=max(max(aux));
            aux(mk==0)=999999;
            bot=min(min(aux));
            dil1=imdilate(mk,temp1);%     figure,imagesc(dil1);
            dil2=imdilate(mk,temp2);%     figure,imagesc(dil2);
            ref2=dil1-dil2;
            bak=min(I11(ref2>0));
            if(bak>bot),bak=bot;end
            thresh=bak+(bot-bak)*thFac;% threshold
            mk3=(I11>thresh).*dil1;
            Lcil=bwlabel(mk3);
            Lcila=cell2mat(struct2cell(regionprops(Lcil,'area')));
            [~,gg]=max(Lcila);
            Lcil(Lcil~=gg)=0;Lcil(Lcil==gg)=1;
            mk3=Lcil;
            %         mk2=bwlabel(I11>thresh);
            %         vals=unique(mk2.*mk);vals(vals==0)=[];
            %         mk3=zeros(size(mk2));
            %         for kk=1:length(vals)
            %             mk3(mk2==vals(kk))=1;%      figure,imagesc(mk3.*I1)
            %         end
        end
        %     if(pinta==1),figure(33);clf,imagesc(mk);axis image;end
        %     [a,b]=find(L==ind(ii));
        %     warning('off','all');
        %     %
        %
        %     % fit exponencial al decay
        %     model1 = fittype('a*exp(-((x-b)/c)^2)+d');
        %     options1 = fitoptions(model1);
        %     options1.Display = 'off';
        %     options1.Robust = 'bisquare';
        %
        %     %
        %     for jj=1:length(a)
        %         disp(jj);
        %         crop=mkr(a(jj)-tam:a(jj)+tam,b(jj)-tam:b(jj)+tam);
        %         [aa,bb]=find(crop==1);
        %         aux=[aa-mean(aa) bb-mean(bb)];
        %         aux=sum(aux.*aux,2);
        %         [~,q]=max(aux);
        %         vec=[aa(q)-mean(aa) bb(q)-mean(bb)];
        %         vec=[vec(2) -vec(1)];
        %         vecn=vec/norm(vec);
        %         vec=round(vecn*dst);
        %         ori=[a(jj) b(jj)]+vec;
        %         fin=[a(jj) b(jj)]-vec;
        %        coor=puntsMig(ori,fin);
        %        dd=diff(coor,[],1);
        %        dd=sqrt(sum(dd.*dd,2));dd=[0;cumsum(dd)];
        %         if(pinta==1),hold on,plot([b(jj)-vec(2) b(jj)+vec(2)],[a(jj)-vec(1) a(jj)+vec(1)],'r');end
        %         if(pinta==2),figure(33);clf;imagesc(crop);axis image;hold on,plot([tam+1-vec(2) tam+1+vec(2)],[tam+1-vec(1) tam+1+vec(1)],'r');pause();end
        %         tra=zeros(length(coor),1);
        %         for kk=1:length(coor)
        %            tra(kk)=I1(coor(kk,1),coor(kk,2));
        %         end
        %             %                        a         b        c        d
        %             aux=mean(dd);
        %     options1.startpoint = [1 aux aux/2 0];
        %     options1.Lower = [0 aux*.5 aux/100 0];
        %     options1.Upper = [2*max(tra) aux*1.5 aux max(tra)];
        %     [fresult1,gof2,out2]  = fit(dd,tra,model1,options1);
        %     plusminus=sigmaFac*fresult1.c/2;
        %
        %         if(pinta==3)
        %             figure(100);clf;subNM(1,2,1,.05);imagesc(I1);set(gca,'xtick',[],'ytick',[]);axis image;
        %             hold on,h=plot([b(jj)-vec(2) b(jj)+vec(2)],[a(jj)-vec(1) a(jj)+vec(1)],'r');set(h,'linewidth',2);
        %             subNM(1,2,2,[.05 .3 .05 .3]);h=plot(dd,tra);set(h,'linewidth',2);hold on;
        %             h=plot(fresult1);set(h,'linewidth',2);
        %             h=plot([fresult1.b-plusminus fresult1.b-plusminus],ylim);set(h,'linewidth',2,'color',[0.9290 0.6940 0.1250]);
        %             h=plot([fresult1.b+plusminus fresult1.b+plusminus],ylim);set(h,'linewidth',2,'color',[0.9290 0.6940 0.1250]);
        %             h=plot([fresult1.b+plusminus fresult1.b-plusminus],[fresult1.d fresult1.d]);set(h,'linewidth',2,'color',[0.9290 0.6940 0.1250]);
        %             legend off;xlabel('pix');ylabel('Counts');
        %             h=text(fresult1.b,fresult1.d,[niceNums(plusminus*2*DX*1e9,0) 'nm']);set(h,'vertical','bot','horiz','center');
        %             xlim([dd(1) dd(end)]);
        %             [gg,~]=getframe(gcf);
        %             imwrite(gg,['slices' filesep niceDigits(jj,3) '.png']);
        %         end
        %
        %         gg=round(ori-vecn.*(fresult1.b-plusminus));
        %         mk(gg(1),gg(2))=1;
        %         gg=round(ori-vecn.*(fresult1.b+plusminus));
        %         mk(gg(1),gg(2))=1;
        % for kk=1:length(coor)
        %
        % end
        %     end
        if(ii==1),mk11=mk3;end
        %         if(ii==2),mk12=mk3;end
    end
    mk=mk11;
    mk(mk12==1)=2;
    
    if(pinta==4)
        ref=[1 0 1;0 1 0;0 0 0];
        II=cat(3,NormArray(I1),NormArray(I2)*1.2);
        II(:,:,3)=II(:,:,1);
        
        aux=mean(II,3);I0=aux;
        for cc=1:3
            ch=II(:,:,cc);
            ch(mk==0)=aux(mk==0);
            II(:,:,cc)=ch;
        end
        [II,MAP] = rgb2ind(II,256);
        
        aux=MAP-ref(1,:);[~,q1]=min(sum(aux.*aux,2));
        aux=MAP-ref(2,:);[~,q2]=min(sum(aux.*aux,2));
        aux=MAP-ref(3,:);[~,q3]=min(sum(aux.*aux,2));
        mk2=mk;
        mk2(mk==0)=q3-1;
        mk2(mk==1)=q1-1;
        mk2(mk==2)=q2-1;
        
        imwrite(uint8(256*NormArray(I0)),gray(256),[RF filesep nom '.gif'],'Loopcount',inf);
        imwrite(II,MAP,[RF filesep nom '.gif'],'WriteMode','append');
        imwrite(uint8(mk2),MAP,[RF filesep nom '.gif'],'WriteMode','append');
        % mk3=mk2/56;mk3=ceil(mk3);
        % imwrite(uint8(mk3),superjet(3,'rgb'),['seg2/' nom '.gif'],'WriteMode','append');
    end
end
end