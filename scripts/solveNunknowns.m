function [weights1,weights2] = solveNunknowns(Ss,Gs,coefs)
% solve system for n harmonics and N components:
% Sn=coefs(4*(n-1)+[1],1)x+coefs(4*(n-1)*[2],1)y+coefs(4*(n-1)*[3],1)z+coefs(4*(n-1)*[4],1)w;
% Gn=coefs(4*(n-1)+[1],2)x+coefs(4*(n-1)*[2],2)y+coefs(4*(n-1)*[3],2)z+coefs(4*(n-1)*[4],2)w;
% x+y+z+w=1
% the coeficients are in variable coefs, 1st column S, 2nd G, first N lines for
% 1st harmonic, and so on
% return two sets of weights:
% both solve by min squares adding the constraint (used to have different methods)


n=[];M=[];comps=size(coefs,1)/size(Ss,1);
for jj=1:size(Ss,1)
    n=[n Ss(jj) Gs(jj)];
    M=[M;coefs([1:comps]+comps*(jj-1),1:2)'];
end

weights1=zeros(comps,1);weights2=weights1;


    M=[M;ones(1,comps)];
    n=[n';1];
 try
 %I=inv(M'*M); 
 I=M'*M;
 
% 
 %M=I*M'; 
 
 M=I\M';
 
 weights2=M*n;

 end
 weights1=weights2;
 %
% weights1=weights1';
% % weights2=weights1-min(weights1);
% % weights2=weights2/sum(weights2);
% % weights2=NormArray(weights1);
% % weights2=weights2/sum(weights2);

end