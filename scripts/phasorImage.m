function [] = phasorImage(PT,FOV,cm,sty,levs,flag)
%    phasorImage(PT,FOV,cm,sty,levs,flag)
% stick phasor histogram in field of view with colormap
% sty is 0 for hist or >0 for contour using filtering in pix specified by sty. if sty<0, both are plotted
% lev is the levels where to show the contour plots (over unity)
% if flag=1 the universalcircle is plotted if flag==2 for spectral phasor lines
if(nargin<6),flag=0;end
if(nargin<4)||isempty(sty),sty=0;end
if(nargin<5)||isempty(levs),levs=[.001,.005,.01,.05,.1,.3,.5,.9];if((~isempty(sty))&&(sty<0)),levs=levs(1:4);end;end

if(nargin<3)||isempty(cm)
% cm=superjet(255,'wcbA');
cm=superjet(255,'wSbZctglyGorrmp');
% cm=superjet(255,'weeggllyyGGaak');
% cm=superjet(255,'wyarp');
% cm=superjet(255,'wyarvk');
% cm=superjet(255,'krayw');
%  cm=superjet(256,'wyarANk');% inferno
% cm=superjet(256,'wyazAu');% plasma
end
if(numel(unique(PT))==1),PT(1,1)=unique(PT)+1;end
clins=inverseColor(cm(1,:));
PTl=PT;
if(size(PT,3)>1),sty=0;end
if(sty<0),PT=gray2rgb256(uint8(256*NormArray(PT)),cm);cm(:,1)=clins(1);cm(:,2)=clins(2);cm(:,3)=clins(3);end
if(sty<=0)
 xI=[FOV(1,1:2);FOV(1,1:2)];
    yI=[FOV(2,2) FOV(2,2);FOV(2,1) FOV(2,1)];
    zI=[0 0; 0 0];
    h=surf(xI,yI,zI,'CData',PT,'FaceColor','texturemap','EdgeColor','none');
end
if(sty~=0)
    hold on;
    [x,y]=meshgrid(linspace(FOV(1,1),FOV(1,2),size(PT,2)),linspace(FOV(2,1),FOV(2,2),size(PT,1)));
    I=imgaussfilt(PTl,abs(sty));I=I(end:-1:1,:);
    I=NormArray(I);
    contour(x,y,I,levs);
end
hold on;
    if(flag==1)
     universalCircle(FOV,clins);
    end
    if(flag==2)
     universalCircle(FOV,clins,-12);
    end    
    niceplot;view(0,90);
    axis image;grid off;set(gca,'color','w');
%     if(size(PT,3)==1)
    colormap(cm);
%     end

    xlim(FOV(1,1:2));ylim(FOV(2,1:2));
if(0)
    h=line(FOV(1,:),[FOV(2,1) FOV(2,1)]);set(h,'color',clins);
    h=line(FOV(1,:),[FOV(2,2) FOV(2,2)]);set(h,'color',clins);
    h=line([FOV(1,1) FOV(1,1)],FOV(2,:));set(h,'color',clins);
    h=line([FOV(1,2) FOV(1,2)],FOV(2,:));set(h,'color',clins);
else
    set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[]);
end
end