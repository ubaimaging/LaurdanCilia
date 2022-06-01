function [mk] = readManualSegmentationMem(ruta,DX,I1,I2,rm,RF)

thFac=.75;% threshold is a factor between the lowest intensity on the cilia ridge and the lowest background around the cilia
I0=imread([ruta]);

if(size(I0,1)~=256)% lack of consistency in experiments, if images have different sizes they need to be equalised
I0=imresize(I0,[256 256],'nearest');
end
nom=ruta(find(ruta==filesep,1,'last')+1:end-4);
radi=ceil(100/(DX*1e9));% swap to nm, numerator is the radius in nm
temp1=getnhood(strel('disk',radi,0));

% 2sigma 68%
% 3sigma 95%
% 4sigma 99.7%
% sigmaFac=3; % full width at X max (insert fraction over unity at which to measure width)
mk=[];
aux=double(I0(:,:,1))-double(I0(:,:,2));
aux=aux';
aux=aux(rm+1:end-rm,rm+1:end-rm);
L=abs(aux);
L(L>0)=1;

% figure, imagesc(L);% figure, imagesc(I0)
if(max(L(:))==0)
ensenya('Could not segment membrane from png','r');
else


    
    mk0=zeros(size(L,1),size(L,2));
    tam=3;% radi de la finestra
    dst=(2/(DX*1e6));% distancei per traçar perfil dintensitat (introduir valor en microns, dst es en pix)
    pinta=4;
    
        
        mk=mk0;
        
        mk(L==1)=1;
        %     figure,imagesc(mk);
       mk=imdilate(mk,temp1);
    
    if(pinta==4)
        ref=[1 0 1;0 1 0;0 0 0];
        II=repmat(NormArray(I1),[1,1,3]);
        Ibg=II(:,:,1);
        II(:,:,2)=II(:,:,2).*(1-mk);
%         figure,imagesc(II);
        [II,MAP] = rgb2ind(II,256);
        

        
        imwrite(uint8(256*NormArray(Ibg)),gray(256),[RF filesep nom '.gif'],'Loopcount',inf);
        imwrite(II,MAP,[RF filesep nom '.gif'],'WriteMode','append');
        imwrite(uint8(mk),[ref(3,:);ref(1,:)],[RF filesep nom '.gif'],'WriteMode','append');
        % mk3=mk2/56;mk3=ceil(mk3);
        % imwrite(uint8(mk3),superjet(3,'rgb'),['seg2/' nom '.gif'],'WriteMode','append');
    end
end
end