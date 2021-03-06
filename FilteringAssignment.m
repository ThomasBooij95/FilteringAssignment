clear all
close all
clc

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
SumSV=S(1:49,1:49);


%Plotting of the SVs
fig1=figure('units','normalized','outerposition',[0 0 1 1])
semilogy(SumSV,'+r')




%[ sigma ] = AOloop_nocontrol(phiIdent{1},SNR,H,G)