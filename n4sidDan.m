function [A, C, K, vaf] = n4sid(phik, N_id, N_val, s, n)
    % Inputs:
    %    phik: output data of system to be identified.
    %    N_id: number of data points userd for identification.
    %    N_val: number of data points userd for validation.
    %    s: number of rows used in the hankel matrices.
    %    n: the order of the estimated system.
    % Outputs:
    %    A: the estimated systems A matrix.
    %    C: the estimated systems C matrix.
    %    K: the estimated systems Kalman gain.
    %    vaf: Variance accounted for of estimated system.
    
    n_phi = size(phik,1);
    
    %Construct Hankel matrices
    Phi0_ident = [];
    for i = 1:(N_id - 2*s + 1)
        column = phik(:,i:s+i-1);
        Phi0_ident = [Phi0_ident column(:)];
    end

    Phis_ident = [];
    for i = s + 1:(N_id - s + 1)
        column = phik(:,i:s+i-1);
        Phis_ident = [Phis_ident column(:)];
    end

    %QR factorization
    [~, R] = qr([Phi0_ident; Phis_ident]');
    R = R';
    
    R22 = R(1:size(R,1)/2,1:size(R,1)/2);
    R32 = R(size(R,1)/2 + 1:end,1:size(R,1)/2);

    %Singular value decomposition
    [~, S, V] = svd((R32/R22)*Phi0_ident);
    V = V(:,1:n);
    S = S(1:n,1:n);
    Xs = sqrtm(S)*V';
    
    X0 = Xs(:,1:end - 1);
    X1 = Xs(:,2:end);
    Ys_1_N1 = Phis_ident(1:n_phi,1:size(Phis_ident,2) - 1);
    
    %Calculate system matrices
    A = ((X0*X0')\X0*X1')';
    C = ((X0*X0')\X0*Ys_1_N1')';

    %Calculate Kalman gain
    W = X1 - A*X0;
    V = Ys_1_N1 - C*X0;

    Q_cov = W*W'/N_id;
    S_cov = W*V'/N_id;
    R_cov = V*V'/N_id;
    
    [~,~,K] = dare(A',C',Q_cov,R_cov,S_cov);
    K = K';

    %Simulate system and calculate VAF
    X_hat = zeros(size(A,1),size(phik,2));
    Y_hat = zeros(size(phik,1),size(phik,2));
    num = 0; den = 0;

    for i = N_id + 1 : N_id + N_val
        X_hat(:,i+1) = (A-K*C) * X_hat(:,i) + K*phik(:,i); 
        Y_hat(:,i) = C*X_hat(:,i);
        num = num + sqrt((Y_hat(:,i) - phik(:,i))'*(Y_hat(:,i) - phik(:,i)));
        den = den + sqrt(phik(:,i)'*phik(:,i));
    end

    vaf = 100 * max(0,1 - num / den);
end