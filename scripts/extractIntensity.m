% extract intensity for manual masks

folder='C:\Users\austi\Documents\1_CiliaLaurdan\Datos_Sep2021';


RF=[folder '_Results'];
if(exist(RF,'dir')==0),mkdir(RF);end
S=dir([folder filesep '*.lsm']);
for ii=1:length(S)
    file=S(ii).name;
    [Data,LSMinfo] = lsmread([folder filesep file]);
    
    DX=LSMinfo.voxSizeX;
    Ims=permute(Data,[5,4,2,1,3]);%Y,X,C
    Im=uint8(255*NormArray(sum(Ims,3)));
    imwrite(repmat(Im,[1 1 3]),[RF filesep S(ii).name(1:end-4) '.jpg']);
end
    %     figure, bar([1:28],sum(sum(permute(Ims,[3 1 2]),3),2));
    

    