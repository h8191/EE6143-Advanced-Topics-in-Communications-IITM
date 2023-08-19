%{
This script is used to simulate an OFDM system with a AWGN channel
according to 5G specifications. Further, LDPC soft decoding and 
Rate matching are implemented to improve BER.

When run, the script produces plots of BER vs Eb/N0 for a given QAMorder
and its corresponding code rates

Author: Harshavardhan P
Expected Runtime: < 60s (per QAMorder)
%} 

clc; 
close all;
clear all;

%% params
QAMorder = 16;

fullScreen      = true; %maximize plots to full screen
saveFigures     = false; %save output plots to current directory
printRunTime    = false; %print the time taken for simulation to complete
tic()

%%LDPC and rate matching params
numBits             = 42000; %number of message bits
maxIterLDPCDecoding = 25;    %max iterations before the ldpc decoder terminates
demod_output_type   = 'approxllr'; %approxllr or llr 
% DL-SCH coding parameters
rv                  = 0;     % Redundancy version, 0-3
nlayers             = 1;     % Number of layers, 1-4 for a transport block

%% Parameters to generate Eb/N0 values for good ldpc plots
init_EbN0               = -20 ;%initial value of Eb/N0 in (dB)
init_step_size_EbN0     =   2 ;%initial value by which Eb/N0 (dB) is increased
BER_threshold           = 1e-3;%minimum BER necessary for the simulation to stop
BER_ylim                = 0.4* BER_threshold;% min ylim of plots

%%rates used for different QAMorder 
rates     = {};
rates{16} = [434/1024, 616/1024];
rates{64} = [466/1024, 873/1024];

%% Defining the 5G NR parameters
CPLen1 = 352;                   % length of the cyclic prefix for the first OFDM symbol
CPLen2 = 288;                   % length of the cyclic prefix for the OFDM symbols 2-14
FFTsize = 4096;                 % size of the FFT being used
SCS = 30e3;                     % Subcarrier spacing
SamplingRate = SCS*FFTsize;     % symbol/sample rate for time domain data
numREperRB = 12;                % number of resource elements per resource block
QAMorderList = [QAMorder];      % the value of M in M-QAM, modulation used for data
QAMencoding = 'gray';           % QAM encoding type
EbN0_list_theory = -20:0.5:25;  % Range of Eb/No in dB for theoretical awgn ber curves

%% Defining Resource Grid (RG) parameters
numRB = 50;                             % number of resource blocks allocated 
numRE = numRB*numREperRB;               % number of REs available for transmission 
%numOFDMSymb will be computed after looking at length of rate matched output

%% Output variables
BER_vs_QAM_theory_awgn  = zeros(length(EbN0_list_theory), length(QAMorderList));
EbN0_vs_QAM             = {};%EbN0 list for each QAM and its code rates
BER_vs_QAM              = {};%BER  list for each QAM and its code rates

% Iterating over different Modulation orders
for M_index = 1:length(QAMorderList)
    QAMorder = QAMorderList(M_index);
    modulation = string(QAMorder) + 'QAM';% Modulation scheme 16QAM, 64QAM
    %% setting seed (to get identical results for each run)
    rng(123);

    for rate_index = 1:length(rates{QAMorder})
        rate = rates{QAMorder}(rate_index);  % Target code rate, 0<R<1

        cbsInfo = nrDLSCHInfo(numBits,rate);
        %disp('DL-SCH coding parameters')
        %disp(cbsInfo)
        
        %% Generating transmit data

        % Random transport block data generation
        TXbits = randi([0 1],numBits,1);

        % Transport block CRC attachment
        tbIn = nrCRCEncode(TXbits,cbsInfo.CRC);

        % Code block segmentation and CRC attachment
        cbsIn = nrCodeBlockSegmentLDPC(tbIn,cbsInfo.BGN);

        % LDPC encoding
        enc = nrLDPCEncode(cbsIn,cbsInfo.BGN);

        % Rate matching and code block concatenation
        outlen          = ceil(numBits/rate);
        dataBits        = nrRateMatchLDPC(enc,outlen,rv,modulation,nlayers);

        numDataBits     = length(dataBits);% length of rate matched output 
        numQAMSymb      = length(dataBits) / log2(QAMorder);
        numOFDMSymb     = ceil(numQAMSymb / numRE);

        %fill the empty REs of resource grid with zeros(bits)
        numFillerBits   = numOFDMSymb * numRE * log2(QAMorder) - numDataBits;
        OFDM_input_bits = [dataBits; zeros(numFillerBits,1)];
        
        %% Generating transmit data

        %generating QAM modulated symbols
        modulatedSymbols = qammod(OFDM_input_bits, QAMorder, QAMencoding, ...
                    "InputType", "bit", "UnitAveragePower", true);

        %% Generating the resource grid from modulatedSymbols
        RG = reshape(modulatedSymbols, numRE, numOFDMSymb);

        %% Generating the FFT grid

        FFTgrid = zeros(FFTsize, numOFDMSymb);

        % Generating the suitable index mapping. This is necessary because of the 
        % IFFT implementation in MATLAB. One needs to rearrange [0,a,b,c,d,e,f,0] 
        % to [d,e,f,0,0,a,b,c]. This occurs only along the first (RE) dimension.
        orgIndexSet = 0:numRE-1;
        newIndexSet = mod(orgIndexSet - numRE/2, FFTsize) + 1;
        FFTgrid(newIndexSet,:) = RG;

        % Plotting the generated FFT grid
        % plotResourceGrid(abs(FFTgrid), "FFT grid", "OFDM symbols", "REs");

        %% Generating the time domain data

        % IFFT
        normalizingFactor = sqrt(FFTsize);
        timeDataParallel = normalizingFactor*ifft(FFTgrid, FFTsize, 1);

        % Add CP and time domain sequence
        timeDataSerial = [];

        % Add CP according to the OFDM symbol number
        for i = 1:numOFDMSymb
            if mod(i-1,14) == 0
                CPlength = CPLen1;
            else
                CPlength = CPLen2;
            end
            timeDataSymbol = timeDataParallel(:,i);
            timeDataSerial = [timeDataSerial; timeDataSymbol(end-CPlength + 1:end,:);...
                timeDataSymbol];
        end

        % Change input shape 
        timeDataSerial = squeeze(timeDataSerial);
        
        %% initialzing parameters for Eb/N0 values generation
        EbN0_vs_QAM{QAMorder, rate_index} = []; %list to store Eb Values
        BER_vs_QAM{ QAMorder, rate_index} = []; %list to store BER Values
        
        iteration           = 1;
        BER                 = 1;    
        BER_old             = 0;    
        EbN0                = init_EbN0;
        EbN0_step_size      = init_step_size_EbN0;
        EbN0_upper_bound    = 100;                %min Eb/N0 for which BER == 0

        % Iterating over Eb/N0 values
        threshold_achieved = false;
        while ~threshold_achieved
            SNRdB_uncoded = EbN0 + 10*log10(log2(QAMorder));

            %%SNR at QAM modulator and demodulator
            %scaling the SNR by rate. (SNR decreases since rate < 1)
            SNRdB_coded   = SNRdB_uncoded + 10*log10(rate);

            %%SNR at AWGN channel
            %OFDM decreases SNR since noise bandwidth increases
            SNRdB_OFDM    = SNRdB_coded - 10*log10(FFTsize/numRE);

            %% Transmitting the time domain data through the AWGN channel
            channelOutput = awgn(timeDataSerial, SNRdB_OFDM, "measured", 1234);

            %% Receiver
            RXinput = channelOutput;

            %% CP strip and serial to parallel
            CPcumulative = 0;
            RXtimeDataParallel = zeros(FFTsize, numOFDMSymb);
            for i = 1:numOFDMSymb
                if mod(i-1, 14) == 0
                    CPlength = CPLen1;
                else
                    CPlength = CPLen2;
                end
                CPcumulative = CPcumulative + CPlength;
                RXtimeDataSymbol = RXinput(CPcumulative + (i-1)*FFTsize + 1: CPcumulative +...
                    i*FFTsize,:);
                RXtimeDataParallel(:,i) = RXtimeDataSymbol;
            end

            %% FFT
            % Scaling performed again because IFFT was also scaled
            RXFFTgrid = 1/normalizingFactor*fft(RXtimeDataParallel, FFTsize, 1);

            %% Mapping FFT grid to RB grid
            % This is necessary to undo the mapping done before hand and to rearrange 
            % [d,e,f,0,0,a,b,c] to [0,a,b,c,d,e,f,0]. This occurs only along the 1st
            % (RE) dimension.

            RXfreqData = RXFFTgrid(newIndexSet, :);

            %% Obtaining the transmitted bits

            RXfreqSerial = reshape(RXfreqData, [],1);

            %% dropping filler qam symbols in RG
            RXfreqSerial = RXfreqSerial(1:numQAMSymb);

            %% soft bits (LogLikelihoodRatio) generation
            RXllr  = qamdemod(RXfreqSerial, QAMorder, QAMencoding, ...
                'OutputType', demod_output_type, ...
                'UnitAveragePower',true,'NoiseVariance', 10 ^ (-SNRdB_coded / 10)) ;

            %dropping fillerbits
            RXllr  = RXllr(1:numDataBits);

            % Rate recovery
            raterec = nrRateRecoverLDPC(RXllr,numBits,rate,rv,modulation,nlayers);

            % LDPC decoding
            decBits = nrLDPCDecode(raterec,cbsInfo.BGN,maxIterLDPCDecoding);

            % Code block desegmentation and CRC decoding
            [blk,blkErr] = nrCodeBlockDesegmentLDPC(decBits,cbsInfo.BGN,numBits+cbsInfo.L);

            % Transport block CRC decoding
            [RXbits, tbErr] = nrCRCDecode(blk,cbsInfo.CRC);

            %% BER
            assert (length(RXbits) == numBits & length(TXbits) == numBits, ...
                'TXbits and/or RXbits have different sizes than expected')
            BER = sum(RXbits ~= TXbits)/numBits;

            %fprintf("EbN0 iter %d, M %d, rate %f EbN0 %f BER %f BERold %f stepSize %d\n" ...
            %    ,iteration, QAMorder, rate, EbN0, BER, BER_old, EbN0_step_size);
            iteration = iteration + 1;

            threshold_achieved  = BER > 0 & BER < BER_threshold;%BER_threshold will be visible in plot
            step_size_too_large = (BER == 0 || BER_old/BER > 5);%BER too low or large dip in BER
            if EbN0 >= 0 && EbN0_step_size >= init_step_size_EbN0 && init_step_size_EbN0 > 1
                EbN0_step_size = 1;
            end

            %% controlling drop/dip in BER per step(for smooth plots)
            if (~ step_size_too_large) || (threshold_achieved)
                EbN0_vs_QAM{QAMorder, rate_index} = [EbN0_vs_QAM{QAMorder, rate_index}, EbN0];
                BER_vs_QAM{ QAMorder, rate_index} = [BER_vs_QAM{ QAMorder, rate_index},  BER];
                EbN0 = EbN0 + EbN0_step_size;
                BER_old = BER;
            else    
                if BER == 0 % BER too low to be measured by simulation
                    EbN0_upper_bound = min(EbN0_upper_bound, EbN0);
                end
                
                EbN0 = EbN0 - EbN0_step_size;   %return to previous EbN0
                while true
                    EbN0_step_size = EbN0_step_size / 2;        %decrease step size
                    if EbN0 + EbN0_step_size < EbN0_upper_bound %more likely to get non zero BER
                        EbN0 = EbN0 + EbN0_step_size;
                        break
                    end
                end
            end  
        end%end of while loop(EbNo generating loop)
    %fprintf("iterations %d %d %d\n", iteration, QAMorder, rate_index);

    %% BER theoretical AWGN curves
    BER_vs_QAM_theory_awgn(:,M_index) =  berawgn(EbN0_list_theory,'qam', QAMorder);
    end %end of rate loop
end

if printRunTime
    toc()
end

%% Plotting the obtained results

legendList = {};
fig = figure;
plot_count = 1;

%OFDM LDPC BER curves
for M_index = 1:length(QAMorderList)
    QAMorder = QAMorderList(M_index);
    for rate_index = 1:length(rates{QAMorder})
        semilogy(EbN0_vs_QAM{QAMorder, rate_index}, ...  
                  BER_vs_QAM{QAMorder, rate_index}, ...
            "o--", 'LineWidth', 2);
        legendValue = string(QAMorderList(M_index)) + ...
            sprintf(" QAM rate = %d/1024", rates{QAMorder}(rate_index) * 1024);
        legendList{1,plot_count} = legendValue;
        plot_count = plot_count + 1;
        hold on; 
    end
end

%theoretical AWGN BER curves 
for M_index = 1:length(QAMorderList)
    semilogy(EbN0_list_theory, BER_vs_QAM_theory_awgn(:,M_index),'LineWidth', 2);
    legendValue = string(QAMorderList(M_index)) + " QAM theoretical AWGN";
    legendList{1,plot_count} = legendValue;
    plot_count = plot_count + 1;
end

grid on
ylim([BER_ylim,1]);
xlabel("$10\log_{10}\left(\frac{E_b}{N_0}\right)$", "Interpreter", "latex", 'FontSize', 22);
ylabel("BER", 'FontSize', 22);
titleValue = sprintf("%d QAM (%s) OFDM LDPC",QAMorder,demod_output_type);
title(titleValue,'FontSize', 22);
legend(legendList);
xticks(-20:1:20)

if fullScreen
    set(gcf, 'Position', get(0, 'Screensize')); %maximize plot size 
end

if saveFigures
    saveas(gcf,sprintf('./OFDM_LDPC_%d_QAM.jpg',QAMorder) );
    %saveas(gcf,sprintf('./OFDM_LDPC_%d_QAM.png',QAMorder) );
end
