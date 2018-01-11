clear all
close all
clc

load turbulenceData.mat
load systemMatrices.mat

[ sigma ] = AOloop_nocontrol(phiIdent{1},SNR,H,G)