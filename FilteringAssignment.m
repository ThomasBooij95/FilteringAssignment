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

for i=1:length(phiSim)

detrend_phi=detrend(phiIdent{i},'constant');
C=cov(detrend_phi');
sigma_e=1/sqrt(SNR);

[ var_nocont ] = AOloop_nocontrol(phiIdent{i},SNR,H,G);
[ var_MVM ] = AOloopMVM(G,H,C,sigma_e,phiIdent{i});

varnc(i)=var_nocont;
varMVM(i)=var_MVM;
end

meannc=ones(size(varnc))*mean(varnc);
meanMVM=(ones(size(varMVM))*mean(varMVM));

%% Plot results
fig2=figure('units','normalized','outerposition',[0 0 1 1])
plot(varnc,'ro')
hold on
plot(meannc,'--r')
plot(varMVM,'ko')
plot(meanMVM,'--k')
legend('No control','No control - mean','MVM','MVM - mean')