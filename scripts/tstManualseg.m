
% test manual segmentation
folder='C:\Users\austi\Documents\1_CiliaLaurdan\Datos_Sep2021';
file='2_todo.lsm';
manual='2_todo c-seleccion.png';


 [Data,LSMinfo] = lsmread([folder filesep file]);
    
    DX=LSMinfo.voxSizeX;
    Ims=permute(Data,[5,4,2,1,3]);%Y,X,C
    %     figure, bar([1:28],sum(sum(permute(Ims,[3 1 2]),3),2));
    

I0=imread([folder filesep manual]);
%I0=permute(I,[2,1,3]);
Is=readManualSegmentation(I0,DX);
mg=[.01];
figure(444);clf;set(gcf,'color','w')
rg=[min(min(sum(Ims,3))),max(max(sum(Ims,3)))];
subNM(1,4,1,mg);imagesc(sum(Ims,3),rg);axis image;title('Intensity');set(gca,'xtick',[],'ytick',[]);
subNM(1,4,2,mg);imagesc(I0);axis image;title('Manual Mask');set(gca,'xtick',[],'ytick',[]);
subNM(1,4,3,mg);imagesc(gray2rgb(Is,superjet(3,'kyr')));axis image;title('Interpreted Mask');set(gca,'xtick',[],'ytick',[]);
subNM(1,4,4,mg);imagesc(sum(Ims,3).*(Is>0),rg);axis image;title('Intensity Masked');set(gca,'xtick',[],'ytick',[]);
% imwrite(Re+1,superjet(3,'kyr'),[RF filesep S(ii).name(1:end-4) '_D_mask.png']);


