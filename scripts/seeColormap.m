function [] = seeColormap(A)

if(max(max(A))>1),A=A/255;end

[a,b]=size(A);
alt=ceil(a/4)+1;
I=zeros(a,alt,b);

for ii=1:b
   
    
    I(:,:,ii)=A(:,ii)*ones(1,alt);
    
end

R=I(:,:,1);
G=I(:,:,2);
B=I(:,:,3);
clear I;

I=zeros(alt,a,3);
I(:,:,1)=R';
I(:,:,2)=G';
I(:,:,3)=B';


mg=.05;
figure(8668);clf;set(gcf,'name','Show Colormap','numbertitle','off');
axes('position',[mg .5+mg 1-2*mg .5-2*mg]);
imagesc((I));%axis equal;
set(gca,'box','on');
%title('low');
%xlabel('high');
set(gca,'xtick',[],'ytick',[],'box','off');
 
axes('position',[mg 0+mg 1-2*mg .5-2*mg]);
hold on;
cc='rgbk';
for ii=1:b
   plot(A(:,ii),cc(ii)); axis tight;
end
xlim([0.5,size(A,1)+.5]);ylim([0,1])    ;
end




