% measure phasor positions 1st and 2nd h
clear all;close all;clc;
% one folder for each tag

fols{1}='C:\Users\austi\Documents\1_CiliaLaurdan\soloFiles';
fols{2}='C:\Users\austi\Documents\1_CiliaLaurdan\soloFiles';
fols{3}='C:\Users\austi\Documents\1_CiliaLaurdan\soloFiles';
fols{4}='C:\Users\austi\Documents\1_CiliaLaurdan\soloFiles';
fols{5}='C:\Users\austi\Documents\1_CiliaLaurdan\soloFiles';
tags={'21C DPPC','37C DOPC','cpJF549','JF646','ACDAN'};

% tags={'solo ACDAN'};

Ithresh=25;% pixels with more than these photons
eb=5;ti=0;rm=floor(eb/2)*(ti);
NH=3;%
NChan=28;%%% achtung, only grabbing first 28 ch (bc we are pooling with old data), we may want to increase this
pintaCursor=1;
spectralrange=[416,450:50:700];
if(spectralrange(end)<416+(NChan*(728-416)/32)),qui=length(spectralrange)+1;else;qui=length(spectralrange);end
spectralrange(qui)=416+(NChan*(728-416)/32);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phs=[512,512];
fov=[-1 1;-1 1];
refpos=zeros(NH*length(fols),2);
type=0;
for kk=1:length(fols)
    %

    S=dir([fols{kk} filesep '*' tags{kk} '*.lsm']);
    N=length(S);
    if(N==0)
        S=dir([fols{kk} filesep '*' tags{kk} '*.R64']);
        N=length(S);
        if(N>0)
            type=2;
        end
    else
        type=1;
    end
    PPs=cell(N,2);
    tams=256;tams=numel(rm+1:tams-rm);
    ST=zeros(tams,tams,N*2);GT=ST;IT=ST(:,:,1:N);
    vec=zeros(NChan,1);
    ensenya(['#####################################'],'d');
    ensenya(['Files for tag ' tags{kk} ':'],'b');
    px=0;
    for ii=1:N
        disp(S(ii).name);
        if(type==1)
            [Data,LSMinfo] = lsmread([fols{kk} filesep S(ii).name]);
            DX=LSMinfo.voxSizeX;
            Ims=permute(Data,[5,4,2,1,3]);%Y,X,C
            Ims=Ims(:,:,1:NChan);
            [S1,G1,I] = fpt(Ims,1);
            [S2,G2,~] = fpt(Ims,2);
            [S3,G3,~] = fpt(Ims,3);            
        else
            [FLIM,FileNames,PathName]=refread([fols{kk} filesep S(ii).name]);
            I=FLIM{1}(:,:,1);
            S1=FLIM{1}(:,:,3).*sin(FLIM{1}(:,:,2)*3.1416/180);
            G1=FLIM{1}(:,:,3).*cos(FLIM{1}(:,:,2)*3.1416/180);
            S2=FLIM{1}(:,:,5).*sin(FLIM{1}(:,:,4)*3.1416/180);% keep also higher harmonic
            G2=FLIM{1}(:,:,5).*cos(FLIM{1}(:,:,4)*3.1416/180);% keep also higher harmonic
            if(NH>2),error('not implemented for flim');end
        end
          [S1,G1]=phasorSmooth(S1,G1,eb,ti);S1=S1(rm+1:end-rm,rm+1:end-rm);G1=G1(rm+1:end-rm,rm+1:end-rm);
          [S2,G2]=phasorSmooth(S2,G2,eb,ti);S2=S2(rm+1:end-rm,rm+1:end-rm);G2=G2(rm+1:end-rm,rm+1:end-rm);
          [S3,G3]=phasorSmooth(S3,G3,eb,ti);S3=S3(rm+1:end-rm,rm+1:end-rm);G3=G3(rm+1:end-rm,rm+1:end-rm);          
          I=I(rm+1:end-rm,rm+1:end-rm);
        if(Ithresh==0)
            [S1] = phasorFillNaNs(S1);
            [G1] = phasorFillNaNs(G1);
            [S2] = phasorFillNaNs(S2);
            [G2] = phasorFillNaNs(G2);
            [S3] = phasorFillNaNs(S3);
            [G3] = phasorFillNaNs(G3);            
        else
            [S1,G1]=phasorFilterIntensity(S1,G1,I,Ithresh);
            [S2,G2]=phasorFilterIntensity(S2,G2,I,Ithresh);
            [S3,G3]=phasorFilterIntensity(S3,G3,I,Ithresh);            
        end
%         figure,phasorPlot2(S1,G1,phs,2,fov);
        ST(:,:,ii)=S1;
        GT(:,:,ii)=G1;
        ST(:,:,ii+N)=S2;
        GT(:,:,ii+N)=G2;
        ST(:,:,ii+2*N)=S3;
        GT(:,:,ii+2*N)=G3;        
        in=find(~isnan(S1));
        px=px+length(in);

        if(type==1)
            for jj=1:length(vec)
                aux=Ims(:,:,jj);aux=aux(in);
                vec(jj)=vec(jj)+sum(aux);
            end
        end
    end
    % apanyoos pel laurdan
    % vec(8:end)=0;
    % vec(1:7)=0;
    % vec(15:end)=0;
    yesno=1;% show phasor data or not
    figure(22),clf;df=.3;set(gcf,'position',[128.2000  390.6000  190*(NH+1)  193.6000]);
    subNM(1,3,1,[.06 .13 .05 .16]);bar(NormArray(vec));title(['Spectral curve']);niceplot;xlim([0,30]);
    h=text(0.05,.01,{[num2str(N) 'files'],[num2str(sum(px)) 'pix'],[num2str(sum(vec)) 'ph']});set(h,'fontsize',8,'horizontal','left','vertical','bot','background',[1 1 1 .5]);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for hh=1:NH
auxS=ST(:,:,[1:N]+(N*(hh-1)));auxG=GT(:,:,[1:N]+(N*(hh-1)));
auxS(isnan(auxS))=[];auxG(isnan(auxG))=[];
    ax1=subNM(1,1+NH,1+hh,[.08 .04 .01 .04]);
    if(yesno==1)
        if(type==1)
            phasorPlot2(auxS,auxG,phs*.25,0,fov);
            universalCircle(fov,'k',spectralrange);xlim([-1,1]);ylim([-1,1]);
        else
            phasorPlot2(auxS,auxG,phs,0,fov);
            universalCircle(fov,'k');xlim([0,1]);ylim([0,.5]);
        end
    else
        if(type==1)
            universalCircle(fov,'k',-1);xlim([-1,1]);ylim([-1,1]);
        else
            universalCircle(fov,'k');xlim([0,1]);ylim([0,.5]);
        end
        niceplot;
    end
 title(['harmonic' num2str(hh)]);
    if(type==1)
        [s1,g1]=fpt(vec,hh);
    else
        ed=0:.005:.5;hh=histc(auxS(:),ed);%figure,bar(ed,hh)
        s1=ed(hh==max(hh));
        ed=0:.005:1;hh=histc(auxG(:),ed);%figure,bar(ed,hh)
        g1=ed(hh==max(hh));
    end
    hold on;
    if(pintaCursor==1)
        h=plot(g1,s1,'ro');set(h,'markersize',22);
        h=plot(g1,s1,'r+');set(h,'markersize',22);
        h=text(g1+df,s1,{['g=' niceNums(g1,4)],['s=' niceNums(s1,4)]});
        set(h,'horizontal','left','vertical','middle');
    end
    ensenya(['Extracted mean coordinates [s,g] for pix with more than ' num2str(Ithresh) 'counts:']);
    disp(['H' num2str(hh) ': ' niceNums(s1,4) ' ' niceNums(g1,4)]);

    refpos(kk+length(fols)*(hh-1),1)=s1;
    refpos(kk+length(fols)*(hh-1),2)=g1;    
  end
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     auxS=ST(:,:,1:N);auxG=GT(:,:,1:N);
%   
%     auxS(isnan(auxS))=[];auxG(isnan(auxG))=[];
%     ax1=subNM(1,3,2);
%     if(yesno==1)
%         if(type==1)
%             phasorPlot2(auxS,auxG,phs*.25,0,fov);
%             universalCircle(fov,'k',spectralrange);xlim([-1,1]);ylim([-1,1]);
%         else
%             phasorPlot2(auxS,auxG,phs,0,fov);
%             universalCircle(fov,'k');xlim([0,1]);ylim([0,.5]);
%         end
%     else
%         if(type==1)
%             universalCircle(fov,'k',-1);xlim([-1,1]);ylim([-1,1]);
%         else
%             universalCircle(fov,'k');xlim([0,1]);ylim([0,.5]);
%         end
%         niceplot;
%     end
%     title('harmonic1');
%     if(type==1)
%         [s1,g1]=fpt(vec,1);
%     else
%         ed=0:.005:.5;hh=histc(auxS(:),ed);%figure,bar(ed,hh)
%         s1=ed(hh==max(hh));
%         ed=0:.005:1;hh=histc(auxG(:),ed);%figure,bar(ed,hh)
%         g1=ed(hh==max(hh));
%     end
%     hold on;
%     if(pintaCursor==1)
%         h=plot(g1,s1,'ro');set(h,'markersize',22);
%         h=plot(g1,s1,'r+');set(h,'markersize',22);
%         h=text(g1+df,s1,{['g=' niceNums(g1,4)],['s=' niceNums(s1,4)]});
%         set(h,'horizontal','left','vertical','middle');
%     end
%     ensenya(['Extracted mean coordinates [s,g] for pix with more than ' num2str(Ithresh) 'counts:']);
%     disp(['H1: ' niceNums(s1,4) ' ' niceNums(g1,4)]);
%     refpos(kk,1)=s1;
%     refpos(kk,2)=g1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     auxS=ST(:,:,N+1:2*N);
%     auxG=GT(:,:,N+1:2*N);
%     auxS(isnan(auxS))=[];
%     auxG(isnan(auxG))=[];
%     ax2=subNM(1,3,3);
%     if(yesno==1)
%         if(type==1)
%             phasorPlot2(auxS,auxG,phs*.25,0,fov);
%             universalCircle(fov,'k',spectralrange);xlim([-1,1]);ylim([-1,1]);
%         else
%             phasorPlot2(auxS,auxG,phs,0,fov);
%             universalCircle(fov,'k');xlim([0,1]);ylim([0,.5]);
%         end
%     else
%         if(type==1)
%             universalCircle(fov,'k',-1);xlim([-1,1]);ylim([-1,1]);
%         else
%             universalCircle(fov,'k');xlim([0,1]);ylim([0,.5]);
%         end
%         niceplot;
%     end
%     title('harmonic2');
%     if(type==1)
%         [s2,g2]=fpt(vec,2);
%     else
%         ed=0:.005:.5;hh=histc(auxS(:),ed);%figure,bar(ed,hh)
%         s2=ed(hh==max(hh));
%         ed=0:.005:1;hh=histc(auxG(:),ed);%figure,bar(ed,hh)
%         g2=ed(hh==max(hh));
%     end
%     hold on;
%     if(pintaCursor==1)
%         h=plot(g2,s2,'ro');set(h,'markersize',22);
%         h=plot(g2,s2,'r+');set(h,'markersize',22);
%         h=text(g2+df,s2,{['g=' niceNums(g2,4)],['s=' niceNums(s2,4)]});
%         set(h,'horizontal','left','vertical','middle');
%     end
%     disp(['H2: ' niceNums(s2,4) ' ' niceNums(g2,4)]);
%     refpos(kk+length(fols),1)=s2;
%     refpos(kk+length(fols),2)=g2;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if(exist('rawPositions','dir')==0),mkdir('rawPositions');end
    saveas(gcf,['rawPositions' filesep tags{kk} '_new2.png']);
end

for hh=1:NH
disp('');
if(hh==1),fra=['refpos=['];else,fra=['refpos=[refpos;'];end
rg=[1:length(fols)]+(length(fols)*(hh-1));
for ii=rg
    fra=[fra niceNums(refpos(ii,1),4) ',' niceNums(refpos(ii,2),4) ';'];
end
fra=[fra(1:end-1) '];% ' num2str(hh) ' harmonic'];
disp(fra);
end

% disp('');
% fra=['refpos=['];
% for ii=1:size(refpos,1)/2
%     fra=[fra niceNums(refpos(ii,1),4) ',' niceNums(refpos(ii,2),4) ';'];
% end
% fra=[fra(1:end-1) '];% first harmonic'];
% disp(fra);
% 
% fra=['refpos=[refpos;'];
% for ii=1+size(refpos,1)/2:size(refpos,1)
%     fra=[fra niceNums(refpos(ii,1),4) ',' niceNums(refpos(ii,2),4) ';'];
% end
% fra=[fra(1:end-1) '];% second harmonic'];
% disp(fra);

