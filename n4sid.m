function [A,C,K,vaf] = n4sid(sk,Nid,Nval,s,n,H,G,C,sigma_e)


%n = size(H,1);      % dimension lifted wavefront
%ns = size(G,1);     % dimension lifted sensor slopes

%%
phiPrep=zeros(size(H,1),Nid);

%% For clean data!!!!!!
phiPrep=sk(:,1:Nid);

%% Preprocessing
% 
% M=C*G'*inv(G*C*G'+sigma_e^2*eye(size(G,1)));
% 
% for k = 1:Nid
%     phiPrep(:,k)=M*sk(:,k);
% end

Yh0=[];
Yhs = [];

for i = 1:(Nid - 2*s + 1)
    column = phiPrep(:,i:i+s-1);
    Yh0 = [Yh0 column(:)];
end


for i = s + 1:(Nid - s + 1)
    column = phiPrep(:,i:i+s-1);
    Yhs = [Yhs column(:)];
end

[~,R]=qr([Yh0;Yhs]');
R=R';
R22 = R(1:size(R,1)/2,1:size(R,1)/2);
R32 = R(size(R,1)/2 + 1:end,1:size(R,1)/2);

%Singular value decomposition
    [~, S, V] = svd(R32*inv(R22)*Yh0);
    fig1=figure('units','normalized','outerposition',[0 0 1 1])
    semilogy(S,'+r')
    V = V(:,1:n);
    S = S(1:n,1:n);
    Xs = sqrtm(S)*V';
    
    Xs0 = Xs(1:end-1);
    Xs1 = Xs(2:end);
    
    A=Xs1*pinv(Xs0);
    C=Yhs(1:size(H,1),1:size(Yhs,2)-1)*pinv(Xs0);
    
    W=Xs1-A*Xs0;
    V=Yhs(1:size(H,1),1:size(Yhs,2)-1)-C*Xs0;
    
    Q=W*W'/Nid;
    R=V*V'/Nid;
    S=W*V'/Nid;

end

% function [At, Bt, Ct, Dt, x0t, S] = mysubid(y, u, s, n)
%
% Y=hankel(y(1:s),y(s:end));
% U=hankel(u(1:s),u(s:end));
% Pi=eye(length(U))-U'*inv(U*U')*U;
%
% [~,R]=qr([U;Y]');
% R=R';
% R22=R(s+1:end,s+1:2*s);
%
% [U,S,V] = svd(R22);
% S=diag(S)
% [U2,S2,V2] = svd(Y*Pi*Y');
% S2=sqrt(diag(S2));
% RQSVvsSVD=[S,S2]
%
% Ct=U(1,1:n);
% At=U(1:s-1,1:n)\U(2:s,1:n);
%
% Z=[];
% ZZ=[];
% ZZZ=[];
% for i=1:length(y)
%     Z0=zeros(1,n);
%     for j=1:i-1
%   Z0=Z0+u(j)*Ct*At^(i-j-1);
%     end
%     ZZZ=[ZZZ;u(i)];
%     ZZ=[ZZ;Z0];
%     Z=[Z;Ct*At^(i-1)];
% end
% F=[Z,ZZ,ZZZ];
% Theta=pinv(F)*y;
%
% x0t=Theta(1:n);
% Bt=Theta(n+1:2*n);
% Dt=Theta(end);
%
% end
%
%
%
% end
%
%
%
%
%
%
%
%
