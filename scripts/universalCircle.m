function [] = universalCircle(fov,coloret,freq)
%  universalCircle(fov,coloret,freq);
% plot universal circle, if field of view is specified also plot borders
% if freq is specified also plot lifetimes values in fibonacci fashion
% if freq is negative, draw spectral circle instead (if non integer, first decimal is number of divisions in radius)
% if freq is a vector, assume it is the values of wavelength, first and last being the range of the transform
if((nargin<2)||(isempty(coloret))),coloret=[0 0 0];end
if((nargin<1)||(isempty(fov))),fov=[0 1;0 .5];end
fs=8;hp={'left','center','right'};vp={'top','mid','bot'};
hold on;
circles=1;% number of circles around the spectral phasor
if(nargin>2)
    if(numel(freq)==1)&&(freq<0)
        sectors=-ceil(freq);
        if(round(freq)~=freq)
circles=round(10*(-freq-floor(-freq)));
        end
    else
        sectors=6;
    end
else
    sectors=6;
end

divra=linspace(0,1,circles+1);divra(1)=[];
rang=2*3.1416*[0:sectors]/sectors;rang=rang(1:sectors);ranglabs=zeros(sectors,1)-1;
if(nargin>=3)&&(numel(freq)>1),ranglabs=freq;rang=2*3.1416*(ranglabs-ranglabs(1))/(ranglabs(end)-ranglabs(1));freq=-1;end
if(nargin<3)||(freq>0)
for ii=1:180
      h=plot(.5*[1+cos((ii)*3.1416/180) 1+cos((ii-1)*3.1416/180)],.5*[sin((ii)*3.1416/180) sin((ii-1)*3.1416/180)]);
      set(h,'color',coloret);   
end
else
    
    sp=3;
 for ii=1:sp:360
     for jj=1:length(divra)
h=plot(divra(jj)*[cos((ii)*3.1416/180) cos((ii-sp)*3.1416/180)],divra(jj)*[sin((ii)*3.1416/180) sin((ii-sp)*3.1416/180)]);
      set(h,'color',coloret);       
     end
%           h=plot(.25*[cos((ii)*3.1416/180) cos((ii-1)*3.1416/180)],.25*[sin((ii)*3.1416/180) sin((ii-1)*3.1416/180)]);
%       set(h,'color',coloret);
%     h=plot(.75*[cos((ii)*3.1416/180) cos((ii-1)*3.1416/180)],.75*[sin((ii)*3.1416/180) sin((ii-1)*3.1416/180)]);
%       set(h,'color',coloret);
%           h=plot(.5*[cos((ii)*3.1416/180) cos((ii-1)*3.1416/180)],.5*[sin((ii)*3.1416/180) sin((ii-1)*3.1416/180)]);
%       set(h,'color',coloret);
 end   
 for ii=1:numel(rang)
h=line([0 cos(rang(ii))],[0,sin(rang(ii))]);set(h,'color',coloret); 
if(ranglabs(ii)>=0)
h=text(.99*cos(rang(ii)),.99*sin(rang(ii)),[num2str(ranglabs(ii)) 'nm']);
ang=180*rang(ii)/3.1416;
if((rang(ii)<3.1416)),vpp=3;else,vpp=1;end
if((rang(ii)<3*3.1416/2)&&(rang(ii)>3.1416/2)),hpp=1;ang=ang-180;else,hpp=3;end
set(h,'color',coloret,'fontname','calibri','fontsize',fs,'horizontalalign',hp{hpp},'verticalalign',vp{vpp},'margin',.1,'Rotation',ang); 
end
 end
end
if(nargin>=3)&&(freq>0)
   taus=[.1 .5 1 2 3 5 8 13 21 34]*1e-9;
   omega=2*3.1416*freq;
   gg=1./(1+(omega*taus).^2);
   ss=(omega*taus)./(1+(omega*taus).^2);
   pseudo=atan2(ss,(gg-.5));
   ra=.01;
   ra=[.5-ra .5+ra];
   extraspace=.02;
%    colordefons=inverseColor(coloret);
   for ii=1:length(taus)
      h=plot([.5+ra(1)*cos(pseudo(ii)) .5+ra(2)*cos(pseudo(ii))],[ra(1)*sin(pseudo(ii)) ra(2)*sin(pseudo(ii))]);
      set(h,'color',coloret);
      coords=[.5+.5*cos(pseudo(ii))  .5*sin(pseudo(ii))];
          hpp=1;vpp=2;exv=0;
    % fancy positioning of text %
    if(coords(1,1)>.5),hpp=3;end
    if(coords(1,2)>.3),hpp=4-hpp;end
    if(coords(1,2)>.461),hpp=2;vpp=1;end
%     if(coords(1,2)>.46),vpp=1;end
     h=text(coords(1,1),coords(1,2),{['  ' num2str(taus(ii)*1e9) 'ns  ' ]});
%        h=text(coords(1,1),coords(1,2),{[' \tau=' num2str(taus(ii)*1e9) 'ns ' ]});
        set(h,'color',coloret,'fontname','calibri','fontsize',fs,'horizontalalign',hp{hpp},'verticalalign',vp{vpp},'margin',.1);
%      set(h,'BackgroundColor',[colordefons .5]);
   end
end
    xl=fov(1,:);
    yl=fov(2,:);
if(nargin==0)
h=line(xl,[yl(1) yl(1)]);set(h,'color',coloret);
h=line(xl,[yl(2) yl(2)]);set(h,'color',coloret);
h=line([xl(1) xl(1)],yl);set(h,'color',coloret);
h=line([xl(2) xl(2)],yl);set(h,'color',coloret);
end
xlim(xl);ylim(yl);
%axis image;
end