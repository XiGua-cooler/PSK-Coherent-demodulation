function NoiseSig = Channel(Signal, SNR, CluLen)

% h1=[0.005 , 0.009, -0.0024, 0.854, -0.218, 0.049 ,-0.016]; %典型电话信道
% h2=[1 0.5 0.25 0.125];                                                     %普通信道
% h3=[0.0410+0.01091i, 0.0495+0.01231i, 0.0672+0.01701i, 0.0919+0.0235, 0.7920+0.12811i, 0.3960+0.08711i, 0.2715+0.04981i, 0.2291+0.04141i, 0.1287+0.01541i, 0.1032+0.01191i];
% h4=[0.3132, -0.104, 0.8908, 0.3134];   %[0.3132  - o.  104  0.8908  0.3134JC
% hn=h4;
% Signal=conv(hn,Signal);

clu = zeros(1,CluLen);
ModSignal = [clu Signal clu];
power = sum(abs(ModSignal))/size(ModSignal,2);
noise_power = power/(10^(SNR/10));
noise = sqrt(noise_power/2)*(randn(1,size(ModSignal,2))+1i*randn(1,size(ModSignal,2)));
NoiseSig = ModSignal+noise;
end