clc
clear all
close all   

SNR        = 15;
CluLen     = 0;
Band       = 1e4;
CarrFre    = 8e6;
fs         = 10*CarrFre;
M          = 2;      %BPSK
SymNum     = 5000;
SymData    = randi([0, M-1], 1, SymNum); 

%   Generate BPSK modulated signal
[ bpskModulatedSignal, basebandSignal ] = IQMpsk( SymData, SymNum, M, CarrFre, Band, fs );
%   Add noise to the modulated signal
bpskReceive = Channel(real(bpskModulatedSignal), SNR, CluLen);
%   Demodulate the receive signal
[ ModSignal ] = bpskDemodulation( real(bpskReceive), SymNum, M, CarrFre, Band, fs );
subplot(2,1,1);
plot(ModSignal(4790:4820));
subplot(2,1,2);
plot(basebandSignal(4790:4820));




