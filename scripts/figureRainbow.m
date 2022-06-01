function [panel] = figureRainbow(S,G,I,comps,type)
% make rainbow phasor image
% if comps are specified theyre drawn on 

sz=[512,512];
if(type==1)
fov=[-1 1;-1 1];centre=[0 0];centre=[mean(mean(G)) mean(mean(S))];
fla=2;
else
    fov=[0,1;0 .5];centre=[0.5 .25];centre=[mean(mean(G)) mean(mean(S))];
    fla=1;
end



% sz=[128,256];
if((size(S,1)~=size(I,1))||(size(S,2)~=size(I,2))),error('S,G and I should have same sizes!');end
% centre=[.5 .25];%g,s
Gm=G;Sm=S;
Gm(isnan(G))=[];Sm(isnan(S))=[];

%centre=[mean(Gm(:)) mean(Sm(:))];
sz=round(sz);
[x,y]=meshgrid(linspace(fov(1,1),fov(1,2),sz(2)),linspace(fov(2,2),fov(2,1),sz(1)));

if(nargin>4)
    Id=zeros(size(x,1),size(x,2),size(comps,1));
   for ii=1:size(comps,1) 
    ref=sqrt(((x-comps(ii,2)).^2)+((y-comps(ii,1)).^2));
   end
   centre=[mean(comps(:,2)),mean(comps(:,1))];
end



x=centre(1)-x;
y=centre(2)-y;
circ=superjet(1000,'vbclyorpvbclyorpv');
circ=circ(380:879,:);
% circ=[circ(381:end,:);circ(1:380,:)];%start yellow
circ=[circ(121:end,:);circ(1:120,:)];%start blue
%             figure,seeColormap([circ;circ]);
IC=zeros(sz(1),sz(2),3);
angs=atan2(y,x);
angs=angs+3.141692;
angs=round((angs/6.2834)*(size(circ,1)-1))+1;
mod=sqrt(x.*x+y.*y);
% mod=1-mod./sqrt(2);
FRAC=max(max(mod));
mod=1-mod/FRAC;
mod=mod.*mod.*mod.*mod;
mod(mod>1)=1;
%                  mod=mod*.5
%                  angs=angs(end:-1:1,:);
for c1=1:sz(1)
    for c2=1:sz(2)
        for cc=1:3
            IC(c1,c2,cc)=max(circ(angs(c1,c2),cc),mod(c1,c2));
        end
    end
end
IC(IC>1)=1;
figure(444);clf;imagesc(IC)

[PP,~,~,~]=phasorPlot2(S(:),G(:),sz,0,fov);
phasorImage(IC,fov,[],[],[],fla);%(PT,FOV,cm,sty,levs,flag)

if(nargin>3)
       aux=comps;%[~,ord]=sort(aux(:,2));aux=refpos(ord,:);
    aux(end+1,:)=aux(1,:);
    for kk=1:size(comps,1)
    h=line([aux(kk,2) aux(kk+1,2)],[aux(kk,1) aux(kk+1,1)]);set(h,'color','k','linewidth',1);  
    h=line([aux(kk,2) aux(kk+1,2)],[aux(kk,1) aux(kk+1,1)]);set(h,'color','r','linewidth',.5);  
    end
    for kk=1:size(comps,1)
    h=plot(aux(kk,2),aux(kk,1),'ko','markersize',6,'MarkerFaceColor','r');   
    end
    
end
[xc,yc]=meshgrid(linspace(fov(1,1),fov(1,2),sz(2)),linspace(fov(2,2),fov(2,1),sz(1)));
PP=NormArray(imgaussfilt(PP,4));
aux=NormArray(PP);
[~,h]=contour(xc,yc,aux,[.01 .05 .1 .3 .5 .9]);
colormap(zeros(10,3));
[gg,~]=getframe(gcf);

PP=gg;
%         gg=0;
SS=centre(2)-S;
GG=centre(1)-G;
angs=atan2(SS,GG);
angs=angs+3.141692;
angs=round((angs/6.2834)*(size(circ,1)-1))+1;
mod=sqrt(GG.*GG+SS.*SS);
% mod=1-mod./sqrt(2);
mod=1-mod/FRAC;
mod(mod<0)=0;mod(mod>1)=1;
mod=mod.*mod.*mod.*mod;

imcol=zeros(size(I,1),size(I,2),3);
for c1=1:size(imcol,1)
    for c2=1:size(imcol,2)
        if(~isnan(angs(c1,c2)))
            if(~isnan(mod(c1,c2)))
                for cc=1:3
                                                 imcol(c1,c2,cc)=max(circ(angs(c1,c2),cc),mod(c1,c2));
%                     imcol(c1,c2,cc)=circ(angs(c1,c2),cc);
                end
            end
        end
    end
end

I1=imcol;
aux=I(:,:,1);I=I/quantile(aux(:),.99);I(I>1)=1;
gg2=uint8(I.*imcol*256);I2=gg2;
gg2=imresize(gg2,[size(gg,1) round(size(imcol,2)*size(gg,1)/size(imcol,1))]);
%gg2=uint8(gg2*256);

gg1=imresize(imcol,[size(gg,1) round(size(imcol,2)*size(gg,1)/size(imcol,1))]);
gg1=uint8(gg1*256);

% panel=cat(2,gg(:,23:676,:),gg1,gg2);

%%%%%%%%%%%%%%%%%%
gru=20;
I2=cat(1,256*ones(gru,size(I2,2),3),I2,256*ones(gru,size(I2,2),3));
I2=cat(2,256*ones(size(I2,1),gru,3),I2,256*ones(size(I2,1),gru,3));
I2=imresize(I2,[size(I2,1)*size(PP,2)/size(I2,2),size(PP,2)]);
panel=cat(1,PP,I2);
% figure(1),imagesc(panel);axis image;
end