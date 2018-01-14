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

%Single value decomp
[U,S,V] = svd(G);
SumSV=S(1:rk,1:rk);


%Plotting of the SVs
% fig1=figure('units','normalized','outerposition',[0 0 1 1])
% semilogy(SumSV,'+r')

detrend_phi=detrend(phiIdent{1},'constant');
C=cov(detrend_phi');
sigma_e=1/sqrt(SNR);

[ sigma ] = AOloop_nocontrol(phiIdent{1},SNR,H,G)
[ sigma ] = AOloopMVM(G,H,C,sigma_e,phiIdent{1})