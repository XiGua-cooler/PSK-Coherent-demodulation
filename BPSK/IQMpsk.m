function [ ModSignal, BebandSignal ] = IQMpsk( SymData, SymNum, M, CarrFre, Band, fs )
%MPSK调制
%SymData 数据
%SymNum 数据长度
%M order
% 载波频率
%Band 带宽
%fs 采样率
    Rb = log2(M)*Band;
    Tb = 1/Rb;
    
    numoflength = log(M)/log(2);
    graymat = 0:M-1;
    for i = 1 : M
         mat1 = (dec2bin(  graymat(i),numoflength  )); 
          mat2 = (dec2bin(  floor(graymat(i)/2),numoflength  )); 
        graymat(i) =     bitxor(   (bin2dec(mat1)),(bin2dec(mat2) )   )  ;  %格雷码对照表 
    end
    
    for i = 1 :length(SymData)
        numofp= find(SymData(i) ==graymat);
        rsdata(i) = cos(1/M*pi + (numofp-1)*2*pi/M) + 1j*sin(1/M*pi + (numofp-1)*2*pi/M);  %psk对应的复值
    end
    
    
    t = 0:1/fs:SymNum*Tb-1/fs; 
    ModSignal = zeros(1,length(t));
    for i = 1:length(rsdata)
        t2 =   floor (i/Band * fs) + 1;
        if(t2 > length(t))
            t2 = length(t);
        end
        t1 =   ceil ((i-1)/Band * fs) +1 ;
        ModSignal(t1:t2) = rsdata(i);
    end
   ModSignal = ModSignal .* exp(1i * 2 * pi * CarrFre * t);
   BebandSignal = SymData;
    

end

