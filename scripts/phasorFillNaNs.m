function [S1] = phasorFillNaNs(S1)
% fill NaNs in image with mean of surrounding closer values
% done iteratively in windows of 3x3 until no gaps are left (requires at least 1 pix to be non NaN!)
fill=find(isnan(S1));
[a,b,c]=size(S1);
while ~isempty(fill)
    N=length(fill);
   for ii=N:-1:1
       [y,x,z]=ind2sub([a,b,c],fill(ii));
       rgy=[max(y-1,1):min(y+1,a)];
       rgx=[max(x-1,1):min(x+1,b)];
       ret=S1(rgy,rgx,z);
       ret(isnan(ret))=[];
       S1(fill(ii))=mean(ret);
       fill(ii)=[];
   end
end

end