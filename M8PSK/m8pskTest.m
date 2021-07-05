clc
clear all
close all   

SNR        = 15;
CluLen     = 0;
Band       = 1e4;
CarrFre    = 4e6;
fs         = 10*CarrFre;
M          = 8;      %8PSK
SymNum     = 5000;
SymData    = randi([0, M-1], 1, SymNum);

%   Generate 8PSK modulated signal
[ m8pskModulatedSignal, basebandSignal ] = IQMpsk( SymData, SymNum, M, CarrFre, Band, fs );
%   Add noise to the modulated signal
m8pskReceive=Channel(real(m8pskModulatedSignal), SNR, CluLen);
%   Demodulate the receive signal
[ ModSignal ] = m8pskDemodulation( real(m8pskReceive), SymNum, M, CarrFre, Band, fs );    




