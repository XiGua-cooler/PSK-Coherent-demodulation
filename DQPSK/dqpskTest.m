clc
clear all
close all   

SNR        = 15;
CluLen     = 0;
Band       = 1e4;
CarrFre    = 8e6;
fs         = 10*CarrFre;
M          = 4;      %DQPSK
SymNum     = 8000;
SymData    = randi([0, M-1], 1, SymNum); 


%   Generate DQPSK modulated signal
[ dqpskModulatedSignal, basebandSignal ] = dqpsk( SymData, SymNum, M, CarrFre, Band, fs );
IQMpsk( SymData, SymNum, M, CarrFre, Band, fs );
%   Add noise to the modulated signal
dqpskReceive = Channel(real(dqpskModulatedSignal), SNR, CluLen);
%   Demodulate the modulated signal
[ ModSignal ] = dqpskDemodulation( real(dqpskReceive), SymNum, M, CarrFre, Band, fs );

scatterplot(ModSignal(6:8000));

% subplot(2,1,1);
% plot(outputData(100:200));
% subplot(2,1,2);
% plot(basebandSignal(100:200));
