function [var_e] = phiSid(G,H,A,K,C,SNR,lambda,phik)

n = size(H,1);      % dimension lifted wavefront
ns = size(G,1);     % dimension lifted sensor slopes
T = length(phik);   % number of temporal phase points

%%
epsilon=zeros(size(phik));
sk=zeros(ns,T);
u=zeros(n,T);
var_e=zeros(1,T);
Xval=zeros(size(A,1),T);
%% Initial values
epsilon(:,1)=phik(:,1);
var_e(1)=var(epsilon(:,1));

    
%M=inv(H)*C*G'*inv(G*C*G'+sigma_e^2*eye(size(G,1)));

for k = 1:T-1
    
    Xval(:,k+1)=(A-K*C)*Xval(:,k)+K*phik(:,k);
    phiSim(:,k+1)=C*Xval(:,k+1);    
    u(:,k)=inv(H'*H+lambda*eye(size(H,1)))*H'*C*Xval(:,k+1);
    epsilon(:,k+1)=phik(:,k+1)-H*u(:,k);
    var_e(k+1)=var(epsilon(:,k+1));
end

var_e = mean(var_e);

end
