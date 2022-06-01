addpath('C:\Users\austi\Documents\alexMultiuse');
clear all;close all;clc;
fol='Datos_Sep2021';
RF='Sep2021';if(exist(RF,'dir')==0),mkdir(RF);end
phs=[512,512];
fov=[-1 1;-1 1];
Ithresh=25;% min counts
eb=4;
ti=1;
rm=floor(eb/2)*(ti);%
Nc=4;
do3comps=0;
pseudocm = [{'kbw'}; {'kgw'}; {'kyw'}; {'krw'}];
laumap=superjet(255,'mubcgly');%seeColormap(laumap);
oldmodphase=0;
if(oldmodphase==1)% old values supplied by Florencia
    if(Nc==3)
        refpos=[59.4 .664;102.8 .653;217.04 .866];
    else
        refpos=[59.4 .664;102.8 .653;217.04 .866;319.5 .942];
        refpos=[refpos;243.4 .217;21.2 .313;23.5 .651;172.6 .854];
    end
    refpos(:,1)=3.1416*refpos(:,1)/180;% pas a rad
    refpos=[refpos(:,2).*sin(refpos(:,1)) refpos(:,2).*cos(refpos(:,1))];% pas a cartesian
else% new values straight in cartesian
    refpos=[0.5890 0.4862;0.6705 -0.1324;-0.5313 -0.6855;-0.6476 0.7044];
    refpos=[refpos;0.4490 0.0989;-0.0887 -0.2148;0.5936 0.2815;-0.8338 0.0713];
end
names={'Laurdan no relajado','Laurdan relajado','Cilia JF549','TZ JF647'};

S=dir([fol filesep '*.lsm']);
N=length(S);
PPs=cell(N,2);
tams=numel(rm+1:256-rm);
ST=zeros(tams,tams,N*2);GT=ST;IT=ST(:,:,1:N);
ensenya('Reading. 00%');
for ii=1:N
    
    [Data,LSMinfo] = lsmread([fol filesep S(ii).name]);
    
    DX=LSMinfo.voxSizeX;
    Ims=permute(Data,[5,4,2,1,3]);%Y,X,C
    %     figure, bar([1:28],sum(sum(permute(Ims,[3 1 2]),3),2));
    
    [S1,G1,I] = fpt(Ims,1);
    [S2,G2,~] = fpt(Ims,2);
    
    [S1,G1]=phasorSmooth(S1,G1,eb,ti);
    [S2,G2]=phasorSmooth(S2,G2,eb,ti);
    
    [S1,G1]=phasorFilterIntensity(S1,G1,I,Ithresh);
    [S2,G2]=phasorFilterIntensity(S2,G2,I,Ithresh);
    
    
    ST(:,:,ii)=S1(rm+1:end-rm,rm+1:end-rm);ST(:,:,ii+N)=S2(rm+1:end-rm,rm+1:end-rm);
    GT(:,:,ii)=G1(rm+1:end-rm,rm+1:end-rm);GT(:,:,ii+N)=G2(rm+1:end-rm,rm+1:end-rm);
    IT(:,:,ii)=I(rm+1:end-rm,rm+1:end-rm);
    auxS=S1;auxG=G1;auxS(isnan(auxS))=[];auxG(isnan(auxG))=[];
    aux=phasorPlot2(auxS,auxG,phs,0,fov);
    PPs{ii,1}=aux;
    auxS=S2;auxG=G2;auxS(isnan(auxS))=[];auxG(isnan(auxG))=[];
    aux=phasorPlot2(auxS,auxG,phs,0,fov);
    PPs{ii,2}=aux;
    percentatge(ii,N);
end

NF=max(max(max(IT)));% normalising factor
comps3=zeros(N,1);
% identifica todo/ un o laltre
for ii=1:N   % each file
    if(contains(lower(S(ii).name),'todo'))
        comps3(ii)=-1;
    else
        if(contains(lower(S(ii).name),'jf549'))
            comps3(ii)=1;
        end
        if(contains(lower(S(ii).name),'jf647'))
            comps3(ii)=2;
        end
    end
end
ensenya('Solving. 00%');
for ii=1:N   % each file
    ima=zeros(tams,tams,Nc);
    warning('off','all');
    % solve 4 comps
    for jj=1:tams % each col
        percentatge((ii-1)*tams*3+jj,N*3*tams);
        for kk=1:tams % each pix
            if(Nc==3)
                Si=[ST(jj,kk,ii)]';Gi=[GT(jj,kk,ii)]';% grab pix components
            end
            if(Nc==4)
                Si=[ST(jj,kk,ii) ST(jj,kk,ii+N)]';Gi=[GT(jj,kk,ii) GT(jj,kk,ii+N)]';% grab pix components in harmonics
            end
            sol=solveNunknowns(Si,Gi,refpos);
            for ll=1:Nc
                ima(jj,kk,ll)=sol(ll);% weight matrix
            end
        end
    end
    
    % 4 component figures
    figurePanels(ima,IT(:,:,ii),refpos(1:Nc,:),PPs{ii,1},Nc,NF,pseudocm,fov,[S(ii).name(1:end-4) '_A_h1'],RF,names);
    %     figurePanels(ima,IT(:,:,ii),refpos(Nc+1:end,:),PPs{ii,2},Nc,NF,pseudocm,fov,[S(ii).name(1:end-4) '_h2'],RF,names);
    figureHists(ima,Nc,[S(ii).name(1:end-4)  '_B' ],RF);
    % 2 component figure - obtined by ignoring components 3 & 4 and renormalising to sum of 1 & 2
    figure2comps(ima,IT(:,:,ii),laumap,[S(ii).name(1:end-4)  '_C_2COMPSfrom4'],names,RF);
    
    % extract balance between 3rd and 4th
    if(comps3(ii)==-1)
        try
            [L3,L4,Re] = segment3n4(ima);imwrite(Re+1,superjet(3,'kyr'),[RF filesep S(ii).name(1:end-4) '_D_mask.png']);
            figure2comps(ima,IT(:,:,ii),laumap,[S(ii).name(1:end-4) '_E_2COMPS_Cilia_'],names,RF,(L3==1));
            figure2comps(ima,IT(:,:,ii),laumap,[S(ii).name(1:end-4) '_F_2COMPS_Base_'],names,RF,(L4==1));
            norma=sum(ima(:,:,1:2),3);
            L1=ima(:,:,1)./norma;% norelax
            L2=ima(:,:,2)./norma;% relax
            aux1=L2;aux1(L3==0)=[];
            aux2=L2;aux2(L4==0)=[];
            gg=cell(0);
            gg{1,1}='JF549';gg{1,2}='JF633';
            for ff=1:max(numel(aux1),numel(aux2))
                try,gg{ff+1,1}=aux1(ff);end
                try,gg{ff+1,2}=aux2(ff);end
            end
            muntaCSV([RF filesep S(ii).name(1:end-4) '_G_RelaxedRatio.csv'],',',gg);
        end
    end
    if(comps3(ii)>0)
        if(comps3(ii)==1)
            % solve 3 comps 1 with JF549
            for jj=1:tams % each col
                percentatge((ii-1)*tams*3+jj+tams,N*3*tams);
                for kk=1:tams % each pix
                    Si=[ST(jj,kk,ii)]';Gi=[GT(jj,kk,ii)]';% grab pix components
                    sol=solveNunknowns(Si,Gi,refpos(1:3,:));
                    for ll=1:3
                        ima(jj,kk,ll)=sol(ll);% weight matrix
                    end
                end
            end
            % 3 component figures 1
            figurePanels(ima(:,:,1:3),IT(:,:,ii),refpos(1:3,:),PPs{ii,1},3,NF,pseudocm,fov,[S(ii).name(1:end-4) '_H__' names{3}],RF,names);
            figure2comps(ima(:,:,1:3),IT(:,:,ii),laumap,[S(ii).name(1:end-4) '_H_2COMPSfrom_' names{3}],names,RF);
        end
        if(comps3(ii)==2)
            % solve 3 comps 1 with JF633
            for jj=1:tams % each col
                percentatge((ii-1)*tams*3+jj+2*tams,N*3*tams);
                for kk=1:tams % each pix
                    Si=[ST(jj,kk,ii)]';Gi=[GT(jj,kk,ii)]';% grab pix components
                    sol=solveNunknowns(Si,Gi,refpos([1,2,4],:));
                    for ll=1:3
                        ima(jj,kk,ll)=sol(ll);% weight matrix
                    end
                end
            end
            % 3 component figures 1
            figurePanels(ima(:,:,1:3),IT(:,:,ii),refpos([1,2,4],:),PPs{ii,1},3,NF,{pseudocm{1},pseudocm{2},pseudocm{4}},fov,[S(ii).name(1:end-4) '_J__' names{4}],RF,{names{1},names{2},names{4}});
            figure2comps(ima(:,:,1:3),IT(:,:,ii),laumap,[S(ii).name(1:end-4) '_J_2COMPSfrom_' names{4}],names,RF);
        end
    end
    
    close all;
    
    %     fra=['''Means'',result,''StDevs'',resultdev'];
    %     for kk=1:Nc
    %         fra=[fra ',''RawFrac' num2str(kk) ''',dadexp{' num2str(kk) '}'];
    %     end
    %    % fra=[fra ',''intensity'',FLIM{' num2str(ii)  '}(:,:,' 1  ')'] ;% add intensity
    %     fra=[fra ',''intensity'',FLIM{' num2str(ii)  '}(:,:,1)'] ;% add intensity
    %     eval(['muntaCSV([''Results/FRACT_' FileNames{ii}(1:end-4) '.csv''],'','',' fra ');']);
    %
end
%     figure,phasorPlot2(ST,GT,phs,2,fov);
%      PP=phasorPlot2(S1,G1,phs,0,fov);



% 1er Armónico:
% Laurdan no relajado: P=59.4   M=0.664
% Laurdan relajado:      P=102.8   M=0.653
% Cilia JF549:                P=217.4   M=0.866
% TZ JF633:                 P=319.5   M=0.942
%
% 2o. Armónico:
% Laurdan no relajado: P=243.4   M=0.217
% Laurdan relajado:      P=21.2     M=0.313
% Cilia JF549:                P=23.5     M=0.651
% TZ JF633:                 P=172.6   M=0.854