%expected runtime for code 20secs
clc;clear all; close all;

%----------------------------------------------------------%
%        symb2bit map     |   xy coordinates   |    symb2ind
%        gray coded       |                    |
%    00 01     ==  0 1    |  (-d, d)  (d, d)   |    00  01
%    10 11         2 3    |  (-d,-d)  (d,-d)   |    10  11
%     
%     x coordinate = ( 2*xindex - (rootM+1) )*d
%----------------------------------------------------------%

Eb      = 1;%Energy per bit
mArr    = [4, 16, 64, 256];%symbols per constellation
N0Arr   = logspace(-2.1,-0.2,12);%noise variances 
SNRBarr = 10*log10(Eb./N0Arr);% SNR per bit
results_empirical   = zeros(length(mArr),length(N0Arr),3);
results_theoritical = zeros(length(mArr),length(N0Arr),3);
nSymbPerSimArr = [100000,100000,100000,100000];%number of symbols per simulation for each M
%%%%%nSymbPerSimArr = round(nSymbPerSimArr/10);

for iM = 1:length(mArr)
    M             = mArr(iM);
    rootM         = sqrt(M);
    bitsPerSymb   = log2(M);
    sumOfBits     = sumBitArrFunc(bitsPerSymb);%array containing sum of bits for each number 
                                               %in its binary representation
    constellation = symb2BitGrayMap(M);% gray coding constellations for minimum BER for same BER
    d             = halfMinDistMqam(Eb,M);%half of min distance b/w two symbols in constellation
    nSymbPerSim   = nSymbPerSimArr(iM); % number of symbols for each simulation
    
    %d = 1;
    for iN0 = 1:length(N0Arr)
        N0   = N0Arr(iN0); 
        noiseStd = sqrt(N0/2); %noise standard deviation
        %generate random symbols by randomly choosing x index and y index
        xind  = randi(rootM,1,nSymbPerSim);
        yind  = randi(rootM,1,nSymbPerSim);
        xcoor = ( 2*xind - (rootM+1) )*d;% convert indices to coordinates
        ycoor = ( 2*yind - (rootM+1) )*d;
        xcoor = xcoor + noiseStd*randn(1,nSymbPerSim);%add noise to xcoordinates of symbols
        ycoor = ycoor + noiseStd*randn(1,nSymbPerSim);%add noise to ycoordinates of symbols
        %%%[min(xcoor), max(xcoor), min(ycoor), max(ycoor)]
        
        % find the nearest symbol for each noisy received symbol
        [xind_est, yind_est] = nearestSymbInd(xcoor,ycoor,d,rootM);
        Eb_by_N0_dB = 10*log10(Eb/N0);
        ser = sum(or(xind ~= xind_est,yind ~= yind_est))/nSymbPerSim; % symbol error rate
        
        totalBitErrors = 0;
        totalBits      = bitsPerSymb*nSymbPerSim;
        for iSymb=1:nSymbPerSim
            actualBits     = constellation(xind(iSymb),yind(iSymb));%get the bits for the correct symbol
            predictedBits  = constellation(xind_est(iSymb),yind_est(iSymb));% get bits for predicted symbol
            bitErrors      = bitxor(actualBits, predictedBits);%find the bits that flipped
            totalBitErrors = totalBitErrors + sumOfBits(bitErrors+1);%add the number of bits flipped to total bit errors
        end
        ber = totalBitErrors/totalBits;%bit error rate
        results_empirical(iM,iN0,1)   = Eb_by_N0_dB;%SNR per bit
        results_empirical(iM,iN0,2)   = ber;%bit error rate
        results_empirical(iM,iN0,3)   = ser;%symbol error rate
        
        tmp = ser_theoritical(d,noiseStd,M);
        results_theoritical(iM,iN0,1) = Eb_by_N0_dB;%SNR per bit
        results_theoritical(iM,iN0,2) = tmp/bitsPerSymb;% assuming only one bit error per symbol error
        results_theoritical(iM,iN0,3) = tmp; % theoritical symbol error
        %%%%fprintf("%d %f %f %d %f\n",M,ser, Eb_by_N0_dB, totalBitErrors,ber);
    end
end

%%results_empirical;
for iM = 1:length(mArr)
    semilogy(results_empirical(iM,:,1),results_empirical(iM,:,2));% empirical BER
    hold on;
    semilogy(results_theoritical(iM,:,1),results_theoritical(iM,:,2));%theoretical BER
end

grid on;
ylim([1e-4,1.1]);
title('BER vs E_b/N_o (in dB) for M QAM');
xlabel('10log_{10}(E_b/N_o)');
ylabel('Bit Error Rate(BER)');
legend('M=4e','M=4t','M=16e','M=16t','M=64e','M=64t','M=256e','M=256t');

function y = qfunc_custom(x)
    %custom q function since matlabs inbuilt qfunc needs communication toolbox
    y = 0.5*erfc(x/sqrt(2));
end

function pError = ser_theoritical(d,noiseStd,M)
    q = qfunc_custom(d/noiseStd); % gaussian tail distribution integration
    rootM = sqrt(M);
    pCorrect = (4/M)*((1-q)^2); % symbols on corners
    pCorrect = pCorrect + (4*(rootM-2)/M)*(1-q)*(1-2*q);%symbols on sides
    pCorrect = pCorrect + (((rootM-2)^2)/M)*((1-2*q)^2);%symbols in the center
    pError   = 1 - pCorrect;
end

function [xind,yind] = nearestSymbInd(x,y,d,rootM)
    %return the index of the nearest symbol
    k  = rootM + 1; % constant
    x1 = x./d;
    y1 = y./d;
    xind = round((x1 + k)./2);% rounding off received symbol index to some integer
    yind = round((y1 + k)./2);
    xind = max(min(xind,rootM),1);% bounding the index between 1 and rootM
    yind = max(min(yind,rootM),1);
end

function x = sumBitArrFunc(M)
    M = min(M,10);
    x = zeros(1,2^M);
    x(2) = 1;
    for i=1:M-1
        j = 2^i;
        x(j+1:2*j) = x(1:j) + 1;
    end
end

function symb2bit = symb2BitGrayMap(M)
    % function to gray code symbols so that any two adjacent symbols differ
    % only by one bit
    
    % M should be even power of two
    rootM = round(sqrt(M));
    graycode = grayCodeGenerator(rootM);
    symb2bit = zeros(rootM,rootM);
    count    = 0;
    for i1 = 1:rootM
        for i2 = 1:rootM
            symb2bit(i1,i2) = rootM*graycode(i1) + graycode(i2);
        end
    end
        
end

function graycode = grayCodeGenerator(x)
    % generates x bit gray code
	if x<=1
		graycode = [0,1];
    else
        l   = grayCodeGenerator(x-1);
        tmp = zeros(1,2*length(l));
        tmp(1:length(l))           = l;
        tmp(length(l)+1:length(tmp)) = 2^(x-1) + l(length(l):-1:1);
        graycode = tmp;
    end
end
    
function d = halfMinDistMqam(Eb,M)
	% function to calculate the half of min distance(0.5*2d) for a M-QAM
	% M is an even power of 2  Ex: M = 4,16,64,256,..
    rootM = sqrt(M);
    % calculating average symbol energy if mindistance = 2d\
    % Es = E[R^2] = E[X^2] + E[Y^2] = 2*E[X^2]
    % E[X^2]*rootM = (rootM-1)^2 + .. +  3^2 + 1^2 + 1^2 + 3^2 + .. + (rootM-1)^2 
    % Eb = Es/(log2(M)) = 2*E[X^2]/(log2(M))
    % d  = log2(M)/(2*E[X^2])
    
    E_xsqr = 2*sum( (1:2:rootM).^2 )/rootM;
    E_rsqr = 2* E_xsqr;
    d = sqrt(Eb*log2(M)/E_rsqr);
end

