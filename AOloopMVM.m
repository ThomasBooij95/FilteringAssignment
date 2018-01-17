function [ var_e ] = AOloopMVM(G,H,C,sigma_e,phik)
% Example of online AO simulation for open_loop measurements
% IN
% phik      : incoming turbulence wavefront
% sigma_e   : measurement noise covariance
% C         : variance of the wavefront
% H         : influence matrix mapping the wavefront on the mirror
% G         : measurement matrix
% OUT
% var_e     : variance of the residual wavefront

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

M=inv(H)*C*G'*inv(G*C*G'+sigma_e^2*eye(size(G,1)));

for k = 1:T-1
    u(:,k)=M*sk(:,k);
    epsilon(:,k+1)=phik(:,k+1)-H*u(:,k);
    sk(:,k+1) = awgn(G*phik(:,k+1),1/sigma_e^2);
    var_e(k+1)=var(epsilon(:,k+1));
end

var_e = mean(var_e);

end
