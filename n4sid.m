function [A,C,K,VAF] = n4sid(phiPrep,phiValid,Nid,Nval,s,n)

Yh0=[];
Yhs=[];

%% Hakel matrices
for i = 1:(Nid - 2*s + 1)
    column = phiPrep(:,i:i+s-1);
    Yh0 = [Yh0 column(:)];
end

for i = s + 1:(Nid - s + 1)
    column = phiPrep(:,i:i+s-1);
    Yhs = [Yhs column(:)];
end

%% RQ factorization
[~,R]=qr([Yh0;Yhs]');
R=R';
R22 = R(1:size(R,1)/2,1:size(R,1)/2);
R32 = R(size(R,1)/2 + 1:end,1:size(R,1)/2);

%% Singular value decomposition
[~, S, V] = svd(R32*pinv(R22)*Yh0);

V = V(:,1:n);
S = S(1:n,1:n);
Xs = sqrtm(S)*V';

Xs0 = Xs(:,1:end-1);
Xs1 = Xs(:,2:end);

A=Xs1*pinv(Xs0);
C=Yhs(1:size(phiPrep,1),1:size(Yhs,2)-1)*pinv(Xs0);

W=Xs1-A*Xs0;
V=Yhs(1:size(phiPrep,1),1:size(Yhs,2)-1)-C*Xs0;

Q=W*W'/Nid;
R=V*V'/Nid;
S=W*V'/Nid;

[~,~,K]=dare(A',C',Q,R,S);
K=K';

%% VAF
Xval=zeros(n,Nval);
phiSim=zeros(size(phiPrep,1),Nval);

num=0;
den=0;


%% Simulation using the derived system matrices
for j=Nid+1:Nid+Nval
    Xval(:,j+1-Nid)=(A-K*C)*Xval(:,j-Nid)+K*phiPrep(:,j);
    phiSim(:,j+1-Nid)=C*Xval(:,j+1-Nid);
    num=num+norm(phiValid(:,j-Nid)-phiSim(:,j-Nid))^2;
    den=den+norm(phiValid(:,j-Nid))^2;
end

VAF=100*max(0,1-num/den);



end


