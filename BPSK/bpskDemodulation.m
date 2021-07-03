%------------------------------------------------------------------------------------------------------------
%                  Introduction to the input and output parameters of the function                          |
%   <SymData>           --- (Datatype = Real number)                                                        |
%                           (intrduction = Coherent demodulation of the input modulated signal.)            |
%   <SymNum>            --- (Datatype = Real number)                                                        |
%                           (intrduction = Baseband data length.)                                           |
%                                                                                                           |
%   <fs>                --- (Datatype = Real number)                                                        |
%                           (intrduction = Specifies the sampling frequency of the ADC, the unit is Hz.)    |
%   <Band>              --- (Datatype = Real number)                                                        |
%                           (intrduction = Specifies the symbol rate.)                                      |
%   <CarrFre>           --- (Datatype = Real number)                                                        |
%                           (intrduction = Specifies the carrier frequency , the unit is Hz)                |
%   <ModSignal>         --- (Datatype = Real number)                                                        |
%                           (intrduction = Demodulated baseband signal.)                                    |
%                                                                                                           |
%                                       Function introduction                                               |
%   This function demodulates bpsk.                                                                         |
%------------------------------------------------------------------------------------------------------------
function [ ModSignal ] = bpskDemodulation( SymData, SymNum, M, CarrFre, Band, fs )

    inputDataLength       = length(SymData);
    carrierFrequency      = CarrFre;
    samplingFrequency     = fs;
    symbolRate            = Band;
    samplingDuration      = SymNum / symbolRate; % Defines the total time of sampling.
    %------------------Down-conversion decimation factor------------------%

    extractionFactor      = floor( (1/symbolRate/4) / (1/samplingFrequency) ); 

    %---------------------------------------------------------------------%
    outputDataLength      = floor( inputDataLength/extractionFactor ) ; 
    %------------------------------------------FIR Filter Design------------------------------------------%

        rp = 2;                           % Passband ripple in dB 
        rs = 60;                          % Stopband ripple in dB
        f = [0.5*CarrFre 1.4*CarrFre];    % Cutoff frequencies
        dev = [(10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)];
        a = [1 0];                        % Desired amplitudes
        % Calculate the order from the parameters using FIRPMORD.
        [ firFilterOrder, firNormalizedFreqPoints, firAmplitudeResponse, firWeight ] = firpmord(f,a,dev,fs);
        % Calculate the coefficients using the FIRPM function.
        firFilterCoff = firpm(firFilterOrder, firNormalizedFreqPoints, firAmplitudeResponse, firWeight);
        % Because the filter order may not be equal to the number of coefficients, you need to re-specify the filter length
        firFilterCoffLength  = length(firFilterCoff);

    %------------------------------------------------END--------------------------------------------------%
    cosOutputData                   = zeros(1, inputDataLength);
    sinOutputData                   = zeros(1, inputDataLength);
    extractionSequenceOutputData    = zeros( 1, floor( inputDataLength/extractionFactor ) );
    bufferFirst                     = zeros( 1, 2 ); % Carrier synchronize loop filter buffer.
    bufferSecond                    = zeros( 1, 2 ); % Carrier synchronize loop filter buffer.
    tempInternalOutput              = zeros( 1, outputDataLength );
    phaseErroDetecerOutputData      = zeros( 1, outputDataLength );
    ncoCosOutputData                = zeros( 1, outputDataLength );
    ncoSinOutputData                = zeros( 1, outputDataLength );
    ncoOutputData                   = zeros( 1, outputDataLength );
    mixerOutputData                 = zeros( 1, outputDataLength );
    ncoCosOutputData(1)             = cos(0);
    ncoSinOutputData(1)             = sin(0);

    interpolatedFilterOutputDataQ   = zeros( 1, outputDataLength);
    interpolatedFilterOutputDataI   = zeros( 1, outputDataLength);
    interpolatedFilterBufferQ       = zeros( 1, 4 );
    interpolatedFilterBufferI       = zeros( 1, 4 );
    errDetecerOutputDataQ           = zeros( 1, outputDataLength);
    errDetecerOutputDataI           = zeros( 1, outputDataLength);
    errDetecerOutputData            = zeros( 1, outputDataLength);
    errDetecerBufferQ               = zeros( 1, 5 );
    errDetecerBufferI               = zeros( 1, 5 );
    symbolLoopFilterTempOutputData  = zeros( 1, outputDataLength);
    loopFilterOutputData            = zeros( 1, outputDataLength);
    timingControlModData            = zeros( 1, outputDataLength);
    timingControlBuffer             = zeros( 1, 2 );
    symbolLoopFilterBuffer          = zeros( 1, 2 );
    mk                              = zeros( 1, outputDataLength); % Defined as symbol synchronized resample clock
    uk                              = zeros( 1, outputDataLength); % Defined as the symbol synchronization interpolation filter coefficient.
    resampleOutputData              = zeros( 1, outputDataLength);
    bpskAngle                       = zeros(1, outputDataLength);  % Defined as the demodulation constellation angle.
    basebandSignal                  = zeros(1, outputDataLength);  % Defined as the demodulated baseband signal.


    %----------------------Generate the Local carrier----------------------

        LocalOscPoint = 0 : 1/samplingFrequency : samplingDuration;
        localSinData = sin((2*pi*carrierFrequency)*LocalOscPoint + 0);
        localCosData = sin((2*pi*carrierFrequency)*LocalOscPoint + pi/2);

    %-------------------------------END------------------------------------
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Down Convertion<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    multiplierCosOutputData = localCosData .* [SymData 0];
    multiplierSinOutputData = localSinData .* [SymData 0];
    %-------------------------------Achieve FIR Filter------------------------------%

        for n = firFilterCoffLength:inputDataLength-firFilterCoffLength
            for k = 1:firFilterCoffLength
                cosOutputData(n) = cosOutputData(n) + firFilterCoff(k) * multiplierCosOutputData(firFilterCoffLength-k+n);
                sinOutputData(n) = sinOutputData(n) + firFilterCoff(k) * multiplierSinOutputData(firFilterCoffLength-k+n);
            end
        end

    %----------------------------------------END------------------------------------%
    cnt = 0;
    k   = 1;
    for n = 1:inputDataLength
        if cnt == 0
            extractionSequenceOutputData(k) = complex( sinOutputData(n), cosOutputData(n) );
            k = k + 1;
        end
        cnt = cnt + 1;
        if cnt == extractionFactor
            cnt = 0;
        end
    end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Carrier Synchronize<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
for n = 1:outputDataLength

    %-----------------------Mixer------------------------%

        % !!! Attention: This is a Complex number multiplication !!!
        mixerOutputData(n) = ncoOutputData(n) * extractionSequenceOutputData(n); 

    %-----------------------END------------------------%

    %-----------------------PhaseErroDetecer------------------------%

        Qin = real(mixerOutputData);
        Iin = imag(mixerOutputData);
            
        phaseErroDetecerOutputData(n) = Qin(n)*sign(Iin(n)) - Iin(n)*sign(Qin(n));

    %-----------------------END------------------------%
    %-----------------------LoopFilter------------------------%

        tempInternalOutput(n) = bufferFirst(2) + (0.015) * phaseErroDetecerOutputData(n);
        bufferFirst(1) = (0.015/32)*phaseErroDetecerOutputData(n) + bufferFirst(2);
        bufferSecond(1) = bufferSecond(2) - tempInternalOutput(n);
        loopFilterOutputData(n) = mod(bufferSecond(2), 6.28);

        bufferFirst(2) = bufferFirst(1);
        bufferSecond(2) = bufferSecond(1);

    %-----------------------END------------------------%
    %-----------------------NCO------------------------%
    
        ncoCosOutputData(n+1) = cos(loopFilterOutputData(n));
        ncoSinOutputData(n+1) = sin(loopFilterOutputData(n));

        ncoOutputData(n+1) = complex( (ncoCosOutputData(n+1)), ncoSinOutputData(n+1) );

    %-----------------------END------------------------%
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Symbol Synchronize<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
for n = 1: outputDataLength
    %-----------------------------------------InterpolatedFilter------------------------------------------%

    %-----------------Updata the coefficient of interpolated filter-----------------%
        Coff(1) = (1/6) *uk(n)^3 - (1/6)*uk(n);
        Coff(2) = (-1/2)*uk(n)^3 + (1/2)*uk(n)^2 +       uk(n);
        Coff(3) = (1/2) *uk(n)^3 -       uk(n)^2 - (1/2)*uk(n) + 1;
        Coff(4) = (-1/6)*uk(n)^3 + (1/2)*uk(n)^2 - (1/3)*uk(n);
    %---------------------------------------------------------------------------------%

        interpolatedFilterBufferQ(1) = real(mixerOutputData(n));
        interpolatedFilterBufferI(1) = imag(mixerOutputData(n));
        interpolatedFilterOutputDataQ(n) = Coff(1)*interpolatedFilterBufferQ(1) + Coff(2)*interpolatedFilterBufferQ(2) ...
                                         + Coff(3)*interpolatedFilterBufferQ(3) + Coff(4)*interpolatedFilterBufferQ(4);
        interpolatedFilterOutputDataI(n) = Coff(1)*interpolatedFilterBufferI(1) + Coff(2)*interpolatedFilterBufferI(2) ...
                                         + Coff(3)*interpolatedFilterBufferI(3) + Coff(4)*interpolatedFilterBufferI(4);
        for k = 1:3
            interpolatedFilterBufferQ(5-k) = interpolatedFilterBufferQ(4-k);
            interpolatedFilterBufferI(5-k) = interpolatedFilterBufferI(4-k);
        end

    %------------------------------------------------END--------------------------------------------------%
    %---------------------------------------------ErrDetecer-----------------------------------------------%

        errDetecerBufferQ(1) = interpolatedFilterOutputDataQ(n);
        errDetecerBufferI(1) = interpolatedFilterOutputDataI(n);
        errDetecerOutputDataQ(n) = (errDetecerBufferQ(1)-errDetecerBufferQ(5)) * errDetecerBufferQ(3);
        errDetecerOutputDataI(n) = (errDetecerBufferI(1)-errDetecerBufferI(5)) * errDetecerBufferI(3);
        errDetecerOutputData(n) = errDetecerOutputDataQ(n) + errDetecerOutputDataI(n);
        for k = 1:4
            errDetecerBufferQ(6-k) = errDetecerBufferQ(5-k);
            errDetecerBufferI(6-k) = errDetecerBufferI(5-k);
        end

    %------------------------------------------------END--------------------------------------------------%
    %--------------------------------------------LoopFilter----------------------------------------------%

        if mk(n) == 1
            symbolLoopFilterBuffer(1) = (1/128)*errDetecerOutputData(n);
            symbolLoopFilterTempOutputData(n) = (symbolLoopFilterBuffer(1) + symbolLoopFilterBuffer(2) + errDetecerOutputData(n)) * 1.35e-4; 
            symbolLoopFilterBuffer(2) = symbolLoopFilterBuffer(1);
        end

    %------------------------------------------------END-------------------------------------------------%
    %-------------------------------------------TimingControl---------------------------------------------%

        timingControlModData(n) = mod(timingControlBuffer(2), 1);
        timingControlBuffer(1) = timingControlModData(n) - symbolLoopFilterTempOutputData(n) - (1/4);
        if timingControlBuffer(2) < timingControlModData(n)
            mk(n+1) = 1;
        else
            mk(n+1) = 0;
        end
        if mk(n) == 1
            uk(n+1) = timingControlModData(n);
        end
        timingControlBuffer(2) = timingControlBuffer(1);

    %------------------------------------------------END-------------------------------------------------%
end

%-------------------------------------------Resample--------------------------------------------%

    cnt = 1;
    for n = 1: outputDataLength
        if mk(n) == 1
            resampleOutputData(cnt) = complex( interpolatedFilterOutputDataQ(n), interpolatedFilterOutputDataI(n) );
            cnt = cnt + 1;
        end 
    end

%----------------------------------------------END-----------------------------------------------%
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Symbol Synchronize<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>BPSK Decoding<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    for n = 1: cnt-1
        bpskAngle(n) = atan( imag(resampleOutputData(n))/real(resampleOutputData(n)) );
        if real(resampleOutputData(n)) > 0
            bpskAngle(n) = bpskAngle(n) + pi/2;
        elseif real(resampleOutputData(n)) < 0
            bpskAngle(n) = bpskAngle(n) + 3*pi/2;
        else
            if imag(resampleOutputData(n)) < 0
                bpskAngle(n) = bpskAngle(n) + pi/2;
            end
        end
        if (bpskAngle(n) > pi/4) && (bpskAngle(n) < 3*pi/4)
            basebandSignal(n) = 1;
        else
            basebandSignal(n) = 0;
        end
    end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ModSignal = basebandSignal(1:cnt-1); % This is the output of the final demodulated signal 
end