function [ A,C_w,K] = computeKalmanAR(C_phi0,C_phi1,G,sigma_e)
%B = [A*H -H]

C = G;
A = C_phi1*pinv(C_phi0);
C_w = C_phi0 - A*C_phi0*A';
[~,~,K] = dare(A',C',C_w,sigma_e^2*eye(size(G,1)));
K=K';
end