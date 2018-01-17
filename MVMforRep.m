for k=1:5000
    C_phi0 = C_phi0+...
            detrend_phi(:,k)*detrend_phi(:,k)';      
end
     C_phi0=C_phi0/length(detrend_phi-1);
   
M=inv(H)*C*G'*inv(G*C*G'+sigma_e^2*eye(size(G,1)));

for k = 1:T-1
    u(:,k)=M*sk(:,k);
    epsilon(:,k+1)=phik(:,k+1)-H*u(:,k);
    sk(:,k+1) = awgn(G*phik(:,k+1),1/sigma_e^2);
end