function [ var_e ] = AOloopAR(G,H,C_phi,sigma_e,A,C_w,phik)
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

epsk = zeros(n,T);  % residual wavefront
eps_piston_removed = zeros(n,T); % residual wavefront with mean removed
sk = zeros(ns,T);   % slopes measurements
strehl = zeros(T,1);% strehl ratio
sigma = zeros(T,1);

u      = zeros(n,T) ;% Set initialvalue of the control to zero
phi_DM = zeros(n,T);



for k = 1:T-1
    epsk(:,k+1) = phik(:,k+1);
    eps_piston_removed(:,k+1) = epsk(:,k+1)-mean(epsk(:,k+1)); 
    sigma(k+1) = var(eps_piston_removed(:,k+1));
    SNR = sigma(k+1)/sigma_e^2;
    sk(:,k+1) = awgn(G*epsk(:,k+1),SNR);
    
    u(:,k) = inv(H'*H)*H'*[A-K*C , A*H , K ]*[epsk(:,k);u(:,k-1); sk(:,k)];  
    phi_DM(:,k+1) = H*u(:,k);% function that 
    residual(:,k+1) = epsk(:,k+1) - phi_DM(:,k+1);
    residual_removed(:,k+1) = residual(:,k+1) - mean(residual(:,k+1));
    var_e(:,k+1) = var(residual_removed(:,k+1));
    strehl(k+1) = exp(-sigma(k+1)^2);
end
    
   %strehl = mean(strehl);
var_e = mean(var_e);

end
