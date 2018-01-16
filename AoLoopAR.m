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
B=[-H A*H];

epsilon(:,1)=phik(:,1);
sk(:,1) = awgn(G*phik(:,1),1/sigma_e^2);
var_e(1)=var(epsilon(:,1));
u(:,1) = inv(H'*H)*H'*[A-K*G , K ]*[epsilon(:,1); sk(:,1)];  
%M=inv(H)*C*G'*inv(G*C*G'+sigma_e^2*eye(size(G,1)));

for k = 2:T
    epsilon(:,k)=phik(:,k)-H*u(:,k-1);
    sk(:,k) = awgn(G*epsilon(:,k),1/sigma_e^2);
    u(:,k) = inv(H'*H)*(H'*[A-K*G , A*H , K ]*[epsilon(:,k); u(:,k-1); sk(:,k)]);
    %epsilon(:,k+1)=phihat(:,k+1)-H*u(:,k);
    var_e(k)=var((epsilon(:,k)-mean(epsilon(:,k))));
end
% SNR = 15;
% 
% s = zeros(size(G,1), size(phi,2));
% s(:,1) = awgn(G*phi(:,1),SNR);
% 
% u = zeros(size(H,1),size(phi,2));
% u(:,1) = (H'H)\H'((A-K*G)*phi(:,1)+K*s(:,1));
% 
% eps = zeros(size(phi,1),size(phi,2));
% eps(:,1) = phi(:,1);
% for k = 2:size(phi,2)
%     eps(:,k) = phi(:,k) - H*u(:,k-1); % epsilon = phi-phi_DM
%     s(:,k) = awgn(G*eps(:,k),SNR);
%     u(:,k) = (H'H)\H'((A-K*G)*eps(:,k)+A*H*u(:,k-1)+K*s(:,k));
% end
% var_e = var(detrend(eps,'constant'));
% var_eps = mean(var_e);



var_e = mean(var_e);

end
