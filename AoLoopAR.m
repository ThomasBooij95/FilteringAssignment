function [ var_e ] = AOloopAR(G,H,sigma_e,A,C_w,phik,K)
% Example of online AO simulation for open_loop measurements
% IN
% phik : incoming turbulence wavefront
% SNR  : Signal to Noise Ratio for the sensor
% H    : influence matrix mapping the wavefront on the mirror
% G    : measurement matrix 
% OUT
% sigma : mean variance of the residual wavefront

n = size(H,1);      % dimension lifted wavefront
ns = size(G,1);     % dimension lifted sensor slopes
T = length(phik);   % number of temporal phase points


%%
epsilon=zeros(size(phik));
sk=zeros(ns,T);
u=zeros(n,T);
var_e=zeros(1,T);

%% Initial values

epsilon(:,1)=phik(:,1);
sk(:,1) = awgn(G*phik(:,1),1/sigma_e^2);
var_e(1)=var(epsilon(:,1));
u(:,1) = inv(H'*H)*H'*([A-K*G , K ]*[epsilon(:,1); sk(:,1)]);  

for k = 2:T
    epsilon(:,k)=phik(:,k)-H*u(:,k-1);
    sk(:,k) = awgn(G*epsilon(:,k),1/sigma_e^2);
    u(:,k) = inv(H'*H)*(H'*[A-K*G , A*H , K ]*[epsilon(:,k); u(:,k-1); sk(:,k)]);
    var_e(k)=var(epsilon(:,k));
end

var_e = mean(var_e);

end
