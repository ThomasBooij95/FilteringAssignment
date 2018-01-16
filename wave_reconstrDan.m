function phi_reconstr = wave_reconstr(phi, G, C_phi, sigma_e, SNR)
    % Inputs:
    %    phi: Atmospheric turbulence measurements.
    %    G: The system matrix G, from the expression
    %       s(k) = G*phi(k) + e(k).
    %    H: The system matrix H, from the expression 
    %       phi_DM(k) = H*u(k-1).
    %    C_phi: the variance of phi(k).
    %    sigma_e: the variance of the noise e(k).
    % Outputs:
    %    phi_reconstr: the reconstructed wavefront.
    
    n = size(G,1);
    s = zeros(n,1);    
    M = C_phi*G'/(G*C_phi*G'+sigma_e^2*eye(size(G,1)));
    
    for k = 1:size(phi,2)
        s(:,k) = awgn(G*phi(:,k), SNR); 
    end
    
    phi_reconstr = M*s;
end