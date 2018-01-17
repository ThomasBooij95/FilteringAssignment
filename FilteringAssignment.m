clear all
close all
%clc

load turbulenceData.mat
load systemMatrices.mat

%%
%3.1
%check rank of G
rk=rank(G)
maxrk=min(size(G))
if rk == maxrk
    disp('Matrix G is full rank!')
else
    disp('Matrix G is NOT full rank!')
end

%% Single value decomp
[U,S,V] = svd(G);
SumSV=S(1:rk,1:rk);


%% Plotting of the SVs
% fig1=figure('units','normalized','outerposition',[0 0 1 1])
% semilogy(SumSV,'+r')

%% Loop for all the data
varnc=[];
varMVM=[];
varAR=[];
varSID=[];
var63=[];
var63nocontrol=[];
VAF=[];

% order of the modeled system
n=343;
% number of rows in the Hankel matrix
s=10;

sk={};

for i=1:length(phiIdent)
    %i=1;
    
    detrend_phi=detrend(phiIdent{i},'constant');
    C_phi0=zeros(size(H));
    C_phi1=zeros(size(H));
    
    for k=1:4999
        C_phi0 = C_phi0+...
            detrend_phi(:,k)*detrend_phi(:,k)';
        C_phi1 = C_phi1+...
            detrend_phi(:,k+1)*detrend_phi(:,k)';
    end
    C_phi0 = C_phi0+...
        detrend_phi(:,5000)*detrend_phi(:,5000)';
    C_phi0=C_phi0/length(detrend_phi-1);
    
    C_phi1=C_phi1/length(detrend_phi-2);
    
    %% Definition of noise variance
    %sigma_e=1/sqrt(SNR);
    s1=awgn(G*phiIdent{i},SNR);
    s2=G*phiIdent{i};
    e=s2-s1;
    sigma_e=mean(var(e));
    
    sk{i} = awgn(G*[phiIdent{i},phiSim{i}],SNR);
    %% Clean data!!!!!
    sk{i} =[phiIdent{i},phiSim{i}];
    
    
    
    %% no control
    [ var_nocont ] = AOloop_nocontrol(phiIdent{i},SNR,H,G);
    %% MVM
    [ var_MVM ] = AOloopMVM(G,H,C_phi0,sigma_e,phiIdent{i});
    %% Vector auto-regressive
    [ A,C_w,K] = computeKalmanAR(C_phi0,C_phi1,G,sigma_e);
    [ var_AR ] = AOloopAR(G,H,sigma_e,A,C_w,phiIdent{i},K);
    %% N4SID

    phiasd=awgn(phiIdent{i},SNR);
    [A,C,K,vaf] = n4sid([phiasd,phiSim{i}],phiSim{i},...
        length(phiIdent{i}),length(phiSim{i}),s,n);
    vaf
    [var_SID] = phiSid(G,H,A,K,C,SNR,0,phiIdent{i})
    
    %% 6.3
    [ var_63 ] = AOloop63(G,H,C_phi0,sigma_e,phiIdent{i})
    var63_nocontrol=mean(var(awgn(G*phiIdent{i},SNR)));
    
    %% Collect variances
    
    varnc(i)=var_nocont;
    varMVM(i)=var_MVM;
    varAR(i)=var_AR;
    varSID(i)=var_SID;
    var63(i)=var_63;
    VAF(i)=vaf;
    var63nocontrol(i)=var63_nocontrol;
end

meannc=ones(size(varnc))*mean(varnc);
meanMVM=(ones(size(varMVM))*mean(varMVM));
meanAR=(ones(size(varAR))*mean(varAR));
meanSID=(ones(size(varSID))*mean(varSID));

%% Plot results

fig2=figure('units','normalized','outerposition',[0 0 1 1])
plot(varnc,'ro')
hold on
plot(meannc,'--r')
plot(varMVM,'ko')
plot(meanMVM,'--k')
plot(varAR,'mo')
plot(meanAR,'--m')
plot(varSID,'bo')
plot(meanSID,'--b')
xlabel('Simulation #')
ylabel('Variance')
ylim([0,11])
title('Performance comparison of the different control methods in terms of the mean variance of the wavefront residuals')
legend('No control','No control - mean','MVM','MVM - mean'...
    ,'AR','AR - mean','SID','SID - mean')

fig1=figure('units','normalized','outerposition',[0 0 1 1])
plot(var63,'ro')
hold on
plot(var63nocontrol,'ko')
xlabel('Simulation #')
ylabel('Variance')
ylim([0,3])
title('Performance comparison of the mean variance control of the residual slopes')
legend('MVM','No control')