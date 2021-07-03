clc
clear all
close all   

SNR        = 15;
CluLen     = 0;
Band       = 1e4;
CarrFre    = 4e6;
fs         = 10*CarrFre;
M          = 4;      %QPSK
SymNum     = 5000;
SymData    = randi([0, M-1], 1, SymNum); 

%   Generate QPPSK modulated signal
[ qpskModulatedSignal, basebandSignal ] = IQMpsk( SymData, SymNum, M, CarrFre, Band, fs );
%   Add noise to the modulated signal
qpskReceive=Channel(real(qpskModulatedSignal), SNR, CluLen);
%   Demodulate the receive signal
[ ModSignal ] = qppskDemodulation( real(qpskReceive), SymNum, M, CarrFre, Band, fs );

subplot(2,1,2);
plot(ModSignal(4950:5000));
subplot(2,1,1);
plot(basebandSignal(4950:5000));





