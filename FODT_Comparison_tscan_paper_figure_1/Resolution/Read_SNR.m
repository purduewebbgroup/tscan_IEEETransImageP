%% Brian Bentz 4/2015
%
% Read in SNR

path(path,'../MAT');

N = 18*18;
SNR = loadsnr('snr',N);
Avg_Snr = sum(SNR)/length(SNR)