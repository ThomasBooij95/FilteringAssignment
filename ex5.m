%clear; close all; clc

%Load provided data
load systemMatrices.mat; load turbulenceData.mat;

n = 1; %Number of datasets

%Preprocess and detrend data (remove mean)
phiIdent_detrend = cell(size(phiIdent));
phiSim_detrend = cell(size(phiSim));

C_phi_0_Ident = cell(1,n);
C_phi_0_Sim = cell(1,n);

C_phi_1_Ident = cell(1,n);
C_phi_1_Sim = cell(1,n);

for i = 1:n
    %Detrended data for covariance calculations
    phiIdent_detrend{1,i} = detrend(phiIdent{1,i}, 'constant');
    phiSim_detrend{1,i} = detrend(phiSim{1,i}, 'constant');
    
    C_phi_0_Ident{1,i} = cov(phiIdent_detrend{1,i}',1);
    C_phi_0_Sim{1,i} = cov(phiSim_detrend{1,i}',1);

    n_Ident = size(phiIdent{1,i},1);
    n_Sim = size(phiSim{1,i},1);
    
    C_phi_1_Ident{1,i} = zeros(n_Ident,n_Ident);
    C_phi_1_Sim{1,i} = zeros(n_Sim,n_Sim);
    
    for k=1:4999
        C_phi_1_Ident{1,i} = C_phi_1_Ident{1,i} + ...
            phiIdent_detrend{1,i}(:,k+1)*phiIdent_detrend{1,i}(:,k)';
        
        C_phi_1_Sim{1,i} = C_phi_1_Sim{1,i} + ...
            phiSim_detrend{1,i}(:,k+1)*phiSim_detrend{1,i}(:,k)';
    end
    
    C_phi_1_Ident{1,i} = C_phi_1_Ident{1,i}/4999;
    C_phi_1_Sim{1,i} = C_phi_1_Sim{1,i}/4999;
end

%Wave reconstruction
phik_Ident = cell(size(phiIdent));
phik_Sim = cell(size(phiSim));
sigma_e = 1/sqrt(SNR);

for i = 1:n
    phi_Ident_reconstr = wave_reconstrDan(phiIdent{1,i},...
                        G,C_phi_0_Ident{1,i},sigma_e,SNR);
    phik_Ident{1,i} = phi_Ident_reconstr;
    
    phi_Sim_reconstr = wave_reconstrDan(phiSim{1,i},...
                      G,C_phi_0_Sim{1,i},sigma_e,SNR);
    phik_Sim{1,i} = phi_Sim_reconstr;
end

%System identification
[A, C, K, vaf] = n4sidDan([phik_Ident{1,i} phik_Sim{1,i}], 5000, 5000, 10, 100);
























