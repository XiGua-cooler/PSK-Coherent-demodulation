function [ ModSignal, BasebandSignal ] = dqpsk(  SymData, SymNum, M, CarrFre, Band, fs )  

%   pi/4 QPSK调制
%   此处显示详细说明

    x = SymData;
    BR = Band;
    Rb = log2(M)*Band/2;
    Tb = 1/Rb;
    STL = SymNum*Tb;
    SR = fs;
    
    if(M ==4)
        numoflength = log(M)/log(2);
        graymat = 0:M-1;
        for i = 1 : M
            mat1 = (dec2bin(  graymat(i),numoflength  )); 
            mat2 = (dec2bin(  floor(graymat(i)/2),numoflength  )); 
            graymat(i) =     bitxor(   (bin2dec(mat1)),(bin2dec(mat2) )   )  ;  %格雷码对照表 
        end

        for i = 1 :length(x)
            numofp= find(x(i) ==graymat);
            rsdata1(i) = cos(1/M*pi + (numofp-1)*2*pi/M) + 1j*sin(1/M*pi + (numofp-1)*2*pi/M);  %1路Qpsk对应的复值
            rsdata2(i) = cos(2/M*pi + (numofp-1)*2*pi/M) + 1j*sin(2/M*pi + (numofp-1)*2*pi/M);  %2路Qpsk对应的复值
            fb1 =  rem((1:length(rsdata1)) ,2  );  %产生2路10交替的方波
            fb2 =  rem(  ( 2:length(rsdata2)+1) ,2  ) ;
           
        end
        rsdata =  rsdata1 .*fb1 +  fb2.* rsdata2;

        t = 0:1/SR:STL-1/SR; 
        ModSignal = zeros(1,length(t));
        for i = 1:length(rsdata)
            t2 =   floor (i/BR * SR) + 1;
            if(t2 > length(t))
                t2 = length(t);
            end
            t1 =   ceil ((i-1)/BR * SR) +1 ;
            ModSignal(t1:t2) = rsdata(i);
        end
   
    else
        ModSignal = x;
    end
    ModSignal = ModSignal .* exp(1i * 2 * pi * CarrFre * t);
    BasebandSignal = SymData;
end

