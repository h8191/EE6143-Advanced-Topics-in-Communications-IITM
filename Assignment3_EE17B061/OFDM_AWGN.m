% % %Expected run time < 30sec
timing_offset = 5;
CFO           = 0;
M_CFO         = 4; % constellation for which CFO effects are plotted
n_iter_CFO    = 5; % calculating BER(CFO) multiple times for accurate plots

%struct containing parameters of OFDM 
OFDM_params                      = struct('n_slots', 1);
OFDM_params.n_ofdm_symb_per_slot = 14;
OFDM_params.n_RB_per_ofdm_symb   = 50;
OFDM_params.n_RE_per_ofdm_symb   = OFDM_params.n_RB_per_ofdm_symb * 12;
OFDM_params.fft_size             = 4096;

Eb       = 1;                            %Energy per bit
M_arr    = [4, 16, 64, 256];             %symbols per constellation
N0_arr   =  logspace(-2.2,0.0,25);      % specifying N0 ranges

%arrays to store BER and SER, theoretical and empirical
results_empirical    = zeros(length(M_arr),length(N0_arr),5);
results_theoretical  = zeros(length(M_arr),length(N0_arr),3);

%% code for generating ber curves with and without timing offset corrected
for ind_M = 1:length(M_arr)
    modulation_params = gen_modulation_params_M_QAM(M_arr(ind_M), Eb); %M-QAM parameters
    
    %generate random message words in range 0 to M - 1
    message_length = OFDM_params.n_ofdm_symb_per_slot * OFDM_params.n_RE_per_ofdm_symb;
    message_bits_tx = randi([0,modulation_params.M-1], 1 ,message_length);
    % each of these M numbers(0 to M-1) is equally likely and these
    % M numbers cover all the log2(M) bit numbers (since M is a power
    % of 2). So each of the log2(M) bits is equally
    
    % start of transmitter logic
    [tx_sequence, ifft_output_arr] = transmitter_logic(message_bits_tx, OFDM_params, modulation_params);
    %transmitted sequence is same for a given M for any SNR
    
    for ind_N0 = 1:length(N0_arr)
        N0          = N0_arr(ind_N0);  %N0 value
        noiseStd    = sqrt(N0/2);      %noise standard deviation per each dimension
        Eb_by_N0_dB = 10*log10(Eb/N0); %SNR per bit
        
        % data transmission through AWGN channel 
        rx_sequence = AWGN_channel(tx_sequence, noiseStd);
        
        % introducing timing offset
        rx_sequence = [zeros(1,timing_offset), rx_sequence];
        
% % %         % uncomment this to include the effect of CFO in timing offset curves
% % %         n = 1:length(rx_sequence);
% % %         cfo_effect_factor = exp(1j*2*pi*CFO*n/OFDM_params.fft_size);
% % %         rx_sequence = rx_sequence .* cfo_effect_factor;
    
        % demodulation with timing offset correction
        [message_bits_rx1] = receiver_logic(rx_sequence, timing_offset, ...
            OFDM_params, modulation_params);
        
        % demodulation with uncorrected timing offset 
        [message_bits_rx2] = receiver_logic(rx_sequence, 0, ...
            OFDM_params, modulation_params);
        
        %calculating ser and ber by comparing received and transmitted bits
        % ber and ser calculation(for timing offset corrected)
        [ser_e, ber_e, ser_t, ber_t] = calculate_error_rates(...
            message_bits_tx, message_bits_rx1, noiseStd, modulation_params);
        results_theoretical(ind_M,ind_N0,:)   = [Eb_by_N0_dB, ser_t, ber_t];
        results_empirical(ind_M,ind_N0,1:3)   = [Eb_by_N0_dB, ser_e, ber_e];
        
        % ber and ser calculation(for uncorrected timing offset )
        [ser_e, ber_e, ser_t, ber_t] = calculate_error_rates(...
            message_bits_tx, message_bits_rx2, noiseStd, modulation_params);
        results_empirical(ind_M,ind_N0,4:5)   = [ ser_e, ber_e];
    end
end

%% 16 QAM highest SNR (timing offset received symbol plots)
low_noise_variance_16QAM         = 10^(-2.2);
modulation_params                = gen_modulation_params_M_QAM(16,Eb);

%generate random message words in range 0 to M - 1 (each bit equally likely) 
message_length  = OFDM_params.n_ofdm_symb_per_slot * OFDM_params.n_RE_per_ofdm_symb;
message_bits_tx = randi([0,modulation_params.M-1], 1 ,message_length);

[tx_sequence] = transmitter_logic(message_bits_tx, OFDM_params, modulation_params);
rx_sequence   = AWGN_channel(tx_sequence, sqrt(low_noise_variance_16QAM));

%creating timing offset
rx_sequence = [zeros(1,timing_offset), rx_sequence];

%demodulation with timing offset corrected
[~, ~ , resource_grid_symb_rx1] = receiver_logic(rx_sequence, timing_offset, ...
    OFDM_params, modulation_params);

%demodulation without correcting timing_offset
[~, ~, resource_grid_symb_rx2] = receiver_logic(rx_sequence, 0, ...
    OFDM_params, modulation_params);

%% section 3.4.2 (Effect of CFO on 4-QAM)
CFO_values         = [0.0005, 0.001, 0.0015, 0.002, 0.003];
modulation_params  = gen_modulation_params_M_QAM(M_CFO,Eb);
CFO_results = zeros(length(CFO_values)+1, length(N0_arr), 3);

%generate random message words in range 0 to M - 1
message_length = OFDM_params.n_ofdm_symb_per_slot * OFDM_params.n_RE_per_ofdm_symb;
message_bits_tx = randi([0,modulation_params.M-1], 1 ,message_length);

[tx_sequence] = transmitter_logic(message_bits_tx, OFDM_params, ...
    modulation_params);

for iter=1:n_iter_CFO
    for ind_CFO = 1:length(CFO_values)
        CFO = CFO_values(ind_CFO);
        for ind_N0 = 1:length(N0_arr)
            N0          = N0_arr(ind_N0);
            noiseStd    = sqrt(N0/2);     %noise standard deviation
            Eb_by_N0_dB = 10*log10(Eb/N0);%SNR per bit
            
            rx_sequence = AWGN_channel(tx_sequence, noiseStd);
            
            % CFO_effect
            n = 1:length(rx_sequence);
            cfo_effect_factor = exp(1j*2*pi*CFO*n/OFDM_params.fft_size);
            rx_sequence = rx_sequence .* cfo_effect_factor;
            
            % receiver
            [message_bits_rx] = receiver_logic(rx_sequence, 0, ...
                OFDM_params, modulation_params);
            
            % ber and ser calculation
            [ser_e, ber_e, ser_t, ber_t] = calculate_error_rates(...
                message_bits_tx, message_bits_rx, noiseStd, modulation_params);
            
            CFO_results(ind_CFO,ind_N0,1) = CFO_results(ind_CFO,ind_N0,1) + Eb_by_N0_dB;
            CFO_results(ind_CFO,ind_N0,2) = CFO_results(ind_CFO,ind_N0,2) + ser_e;
            CFO_results(ind_CFO,ind_N0,3) = CFO_results(ind_CFO,ind_N0,3) + ber_e;
            CFO_results(end,ind_N0,:)     = [Eb_by_N0_dB, ser_t, ber_t];
        end
    end
end
CFO_results(1:end-1,:,:) = CFO_results(1:end-1,:,:)/n_iter_CFO;

plot1 = plot_ber_curves(results_empirical, results_theoretical, M_arr);
% saveas(gcf,'pics/BER_timing_offset.jpg');
% saveas(gcf,'pics/BER_timing_offset.png');
% saveas(gcf,'pics/BER_timing_offset.pdf');

plot2 = plot_received_symbols(resource_grid_symb_rx1, resource_grid_symb_rx2);
% saveas(gcf,'pics/received_symbols.jpg');
% saveas(gcf,'pics/received_symbols.png');
% saveas(gcf,'pics/received_symbols.pdf');

plot3 = plot_CFO_results(CFO_results, CFO_values, M_CFO, n_iter_CFO);
% saveas(gcf,'pics/CFO_plots.jpg');
% saveas(gcf,'pics/CFO_plots.png');
% saveas(gcf,'pics/CFO_plots.pdf');

% figure()
% [p,f,t] = pspectrum( ifft_output_arr(1,:) );
% plot(f, p)

%% functions used 
function handle = plot_ber_curves(results_empirical, results_theoretical, M_arr)
handle = figure();
set(gcf,'position',[0,0,900,450])
subplot(1,2,1)
for ind_M = 1:length(M_arr)
    semilogy(results_empirical(ind_M,:,1),results_empirical(ind_M,:,3) ...
        , 'DisplayName', sprintf('M=%de', M_arr(ind_M)));
    hold on;
    semilogy(results_theoretical(ind_M,:,1),results_theoretical(ind_M,:,3) ...
        , 'DisplayName', sprintf('M=%dt', M_arr(ind_M)));
end
title(sprintf('BER vs E_b/N_o (in dB) for M QAM \n with timing offset corrected (no CFO)'));
xlabel('10log_{10}(E_b/N_o)');
ylabel('Bit Error Rate(BER)');
grid on;
legend();
ylim([1e-4,1.1]);

subplot(1,2,2)
for ind_M = 1:length(M_arr)
    semilogy(results_empirical(ind_M,:,1),results_empirical(ind_M,:,5) ...
        , 'DisplayName', sprintf('M=%de', M_arr(ind_M)));
    hold on;
    semilogy(results_theoretical(ind_M,:,1),results_theoretical(ind_M,:,3) ...
        , 'DisplayName', sprintf('M=%dt', M_arr(ind_M)));
end
title(sprintf('BER vs E_b/N_o (in dB) for M QAM \n uncorrected timing offset (no CFO)'));
xlabel('10log_{10}(E_b/N_o)');
ylabel('Bit Error Rate(BER)');
legend();
grid on;
ylim([1e-4,1.1]);

end

function handle = plot_received_symbols(resource_grid_symb_rx1,resource_grid_symb_rx2)
handle = figure();
real_part1  = real(reshape(resource_grid_symb_rx1.',1,[]));
imag_part1  = imag(reshape(resource_grid_symb_rx1.',1,[]));
real_part2  = real(reshape(resource_grid_symb_rx2.',1,[]));
imag_part2  = imag(reshape(resource_grid_symb_rx2.',1,[]));

subplot(1,2,1);
scatter(real_part1, imag_part1);
title({'received symbols with timing offset corrected' '(no CFO)'});
xlabel('real part'); ylabel('imag part');grid on;
subplot(1,2,2);

scatter(real_part2, imag_part2);
title(sprintf('received symbols with uncorrected timing offset \n (no CFO)'));
xlabel('real part'); ylabel('imag part');grid on;

set(gcf,'position',[0,0,900,400])
end

function handle = plot_CFO_results(CFO_results, CFO_values, M_CFO, n_iter_CFO)
handle = figure();
semilogy(CFO_results(end,:,1), CFO_results(end,:,3), ...
                    'DisplayName','Theoretical');%theoretical BER
hold on;
for ind_CFO = 1:length(CFO_values)
    CFO_value = CFO_values(ind_CFO);
    semilogy(CFO_results(ind_CFO,:,1), CFO_results(ind_CFO,:,3), ...
            'DisplayName',sprintf('CFO = %f', CFO_value) );%empirical BER
end

title(sprintf('BER vs E_b/N_o (in dB) for %d QAM for different CFO \n (no timing offset) %d iterations', ...
            M_CFO, n_iter_CFO));
xlabel('10log_{10}(E_b/N_o)');
ylabel('Bit Error Rate(BER)');
legend
grid on;
ylim([1e-4,1.1]);

end

function modulation_params = gen_modulation_params_M_QAM(M,Eb)
    modulation_params                = struct('M',M);
    modulation_params.root_M         = sqrt(M);
    modulation_params.bits_per_symb  = log2(M);
    modulation_params.d              = half_of_min_dist_of_M_QAM(Eb,M);

    [modulation_params.symb2bits_map, ...
        modulation_params.bits2symb_map] = gen_Symb2BitGrayMaps(M);
    
    modulation_params.sum_of_bits_arr = generate_sum_of_bits_arr(...
                                          modulation_params.bits_per_symb);
end

function [tx_sequence, ifft_output_arr, resource_grid_symb_tx, ...
    message_symb_tx] = transmitter_logic(message_bits_tx, OFDM_params, ...
    modulation_params)
%getting necessary variables
root_M  = modulation_params.root_M;
d       = modulation_params.d;
bits2symb_map = modulation_params.bits2symb_map;

n_RE_per_ofdm_symb   =  OFDM_params.n_RE_per_ofdm_symb;
n_ofdm_symb_per_slot =  OFDM_params.n_ofdm_symb_per_slot;
fft_size             =  OFDM_params.fft_size;

%convert the message bits into complex symbols of M-QAM
symb_indices = bits2symb_map(message_bits_tx+1, :);
symb_coors = ( 2*symb_indices - (root_M+1) )*d;% convert indices to coordinates
message_symb_tx = symb_coors(:,1) + symb_coors(:,2)*1i;% complex msg symbols

%populate the resource grid with generated symbols
resource_grid_symb_tx = reshape(message_symb_tx, ...
    [n_RE_per_ofdm_symb, n_ofdm_symb_per_slot]).';
%reshaping message_symb array row_wise( by doing transpose of columnwise reshape)

%pad with zeros and do ifft
ifft_output_arr = ifft_custom(resource_grid_symb_tx, fft_size);

%cp addition
tx_sequence = zeros(1,61440);
l_used = 0;
for ind=1:n_ofdm_symb_per_slot
    cp_size = 288 + 64*(ind == 1);
    len_tmp = fft_size + cp_size;
    tx_sequence(l_used+1:l_used+len_tmp) = cp_add(ifft_output_arr(ind,:), cp_size);
    l_used  = l_used + len_tmp;
end
end

function [message_bits_rx, message_symb_rx, resource_grid_symb_rx, ...
    fft_input_arr] = receiver_logic(rx_sequence, timing_offset, ...
    OFDM_params, modulation_params)
%parameters used
root_M  = modulation_params.root_M;
d       = modulation_params.d;
symb2bits_map = modulation_params.symb2bits_map;

n_ofdm_symb_per_slot  = OFDM_params.n_ofdm_symb_per_slot;
fft_size              = OFDM_params.fft_size;
n_RE_per_ofdm_symb    = OFDM_params.n_RE_per_ofdm_symb;

%cp removal
fft_input_arr = zeros(n_ofdm_symb_per_slot, fft_size);
l_used = 0;
for ind = 1:n_ofdm_symb_per_slot
    cp_size = 288 + 64*(ind == 1);
    len_tmp = fft_size + cp_size;
    fft_input_arr(ind,:) = cp_remove(rx_sequence(l_used+1:l_used+len_tmp),cp_size);
    l_used = l_used + len_tmp;
end

%do fft and discard sub carriers where we put zeros and timing offset correction
resource_grid_symb_rx = fft_custom(fft_input_arr, n_RE_per_ofdm_symb, timing_offset);

%received message symbols from resource grid
message_symb_rx  = reshape(resource_grid_symb_rx.', 1, []);

%nearest neighbour decoding of symbols
real_part  = real(message_symb_rx);% row coordinate of symbol in constellation
imag_part  = imag(message_symb_rx);% col coordinate of symbol in constellation
%find indices of the nearest symbol in the constellation
[row_indices, col_indices] = get_indices_of_nearest_symbol(real_part,imag_part,d,root_M);

%converting the nearest symbol to bits
linear_indices  = root_M * (col_indices - 1) + row_indices;
message_bits_rx = symb2bits_map(linear_indices);% linear indexing of matlab
end

function rx_sequence = AWGN_channel(tx_sequence, noiseStd)
complex_noise    = noiseStd * complex( randn(1,length(tx_sequence)), ...
                                       randn(1,length(tx_sequence)));
rx_sequence      = tx_sequence + complex_noise;
end

function [ser_e, ber_e, ser_t, ber_t] = calculate_error_rates(...
    message_bits_tx, message_bits_rx, ...
    noiseStd, modulation_params)

error_magnitudes  = bitxor(message_bits_tx, message_bits_rx);
weight_of_errors  = modulation_params.sum_of_bits_arr(error_magnitudes + 1);
ser_e             = mean(weight_of_errors > 0, 'all');
ber_e             = mean(weight_of_errors, 'all')/modulation_params.bits_per_symb;

ser_t = SER_theoretical(modulation_params.d,noiseStd,modulation_params.M);
ber_t = ser_t/modulation_params.bits_per_symb;
end


function y = cp_add(x,cp_length)
%concatenate last cp_length symbols of x to the front of x
y = [x(end -cp_length+1 : end), x];
end

function y = cp_remove(x,cp_length)
y = x(cp_length + 1:end);
end

function y_zeros_removed = fft_custom(arr,final_size,timing_offset)
fft_size = size(arr, 2);%fft size
y = fft(arr, fft_size, 2)/sqrt(fft_size); % row wise fft
if timing_offset ~= 0 %timing offset correction
    correction_factor = exp( (2j*pi*timing_offset/fft_size) * (0 : fft_size - 1) );
    y  = y .* correction_factor;
end
y_zeros_removed = [y(:, fft_size - final_size/2 + 1 : fft_size), ...
                            y(:, 1:final_size/2)]; %resource grid received 
end

function x = ifft_custom(data,fft_size)
len_data = size(data, 2); % data is 2D resource grid 
n_rows_data = size(data, 1);
%pad zeros to increase size of input to ifft. message bits are mapped
%using the mapping given in spec
data_padded = [data(:, len_data/2 + 1 :len_data), ...
               zeros(n_rows_data, fft_size - len_data), data(:,1:len_data/2) ];
x = ifft(data_padded, fft_size, 2)*sqrt(fft_size); % row wise ifft
end

function x = generate_sum_of_bits_arr(nBits)
%returns an array A such that A[i + 1] = number of ones in binary representation of i %

% limiting the number of bits to avoid accidentally exhausing RAM
assert(nBits<=20, "nBits is too many for graycode generation")

x      = zeros(1,2 ^ nBits);
x(2) = 1;
for ind = 1 : nBits - 1
    i = 2^ind;
    x(i+1:2*i) = x(1:i) + 1;
end
end

function graycode = generate_graycode(nBits)
% generates x bit gray code

% limiting the number of bits to avoid accidentally exhausing RAM
assert(nBits<=20, "nBits is too many for graycode generation")

len_graycode    = 2 ^ nBits;
graycode        = zeros(1, len_graycode);
graycode(1 : 2) = [0, 1];

for exponent = 1 : nBits - 1
    i = 2 ^ exponent;
    graycode(i + 1 : 2*i) = bitxor( graycode(i:-1:1), i) ;
end
end

function [xind,yind] = get_indices_of_nearest_symbol(x,y,d,root_M)
%return the index of the nearest symbol
k  = root_M + 1; % constant
x1 = x./d;
y1 = y./d;
xind = round((x1 + k)./2);% rounding off received symbol index to some integer
yind = round((y1 + k)./2);
xind = max(min(xind,root_M),1);% bounding the index between 1 and rootM
yind = max(min(yind,root_M),1);
end

function [symb2bits,bits2symb] = gen_Symb2BitGrayMaps(M)
% function to gray code symbols so that any two adjacent symbols differ
% only by one bit

% M should be 2 power (even number)
rootM    = sqrt(M);
nbits    = log2(rootM);
graycode = generate_graycode(nbits);
symb2bits = zeros(rootM,rootM); %symbol indices to graycoded bits
bits2symb = zeros(M,2);         %bits to symbols indices

for i1 = 1:rootM
    for i2 = 1:rootM
        g1 = graycode(i1); g2 = graycode(i2);
        graycodeword = rootM*g1 + g2;
        symb2bits(i1,i2) = graycodeword;
        bits2symb(graycodeword + 1,:) = [i1,i2];
    end
end
end

function d = half_of_min_dist_of_M_QAM(Eb,M)
% function to calculate the half of min distance(0.5*2d) for a M-QAM
% M is an even power of 2  Ex: M = 4,16,64,256,..
rootM = sqrt(M);
% calculating average symbol energy if mindistance = 2d\
% Es/(d*d) = E[R^2] = E[X^2] + E[Y^2] = 2*E[X^2]
% E[X^2]*rootM = (rootM-1)^2 + .. +  3^2 + 1^2 + 1^2 + 3^2 + .. + (rootM-1)^2
% Eb  = Es/(log2(M)) = 2*E[X^2]*(d*d)/(log2(M))
% d*d = log2(M)/(2*E[X^2])

E_xsqr = 2*sum( (1:2:rootM).^2 )/rootM;
E_rsqr = 2* E_xsqr;
d      = sqrt(Eb*log2(M)/E_rsqr);
end

function y = qfunc_custom(x)
%custom q function since matlabs inbuilt qfunc needs communication toolbox
y = 0.5*erfc(x/sqrt(2));
end

function pError = SER_theoretical(d,noiseStd,M)
q = qfunc_custom(d/noiseStd); % gaussian tail distribution integration
rootM = sqrt(M); %length of side of square
pCorrect = (4/M)*((1-q)^2);                         % symbols on corners
pCorrect = pCorrect + (4*(rootM-2)/M)*(1-q)*(1-2*q);%symbols on sides
pCorrect = pCorrect + (((rootM-2)^2)/M)*((1-2*q)^2);%other symbols
pError   = 1 - pCorrect;
end

