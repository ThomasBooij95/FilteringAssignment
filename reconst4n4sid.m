function [phik]=reconst4n4sid(sk,H,G,C,sigma_e)

%% Preprocessing

M=C*G'*inv(G*C*G'+sigma_e^2*eye(size(G,1)));

for k = 1:Nid
    phiPrep(:,k)=M*sk(:,k);
end