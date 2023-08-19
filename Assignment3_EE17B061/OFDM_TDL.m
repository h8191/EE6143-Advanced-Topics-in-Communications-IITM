% % %expected run time < 60s

fd            = 0;
n_iter        = 20;

%struct containing parameters of OFDM 
OFDM_params                      = struct('n_slots', 1);
OFDM_params.n_ofdm_symb_per_slot = 14;
OFDM_params.n_RB_per_ofdm_symb   = 50;
OFDM_params.n_RE_per_ofdm_symb   = OFDM_params.n_RB_per_ofdm_symb * 12;
OFDM_params.fft_size             = 4096;
OFDM_params.pilot_value          = (1 + 1j)/2;

Eb       = 1;                            %Energy per bit
M_arr    = [4, 16, 64, 256];             %symbols per constellation
N0_arr   =  logspace(-5.0,-0.0,15);      % specifying N0 ranges

%arrays to store BER and SER, theoretical and empirical
average_results_emp  = zeros(length(M_arr),length(N0_arr),3);
results_empirical    = zeros(length(M_arr),length(N0_arr),3);
results_theoretical  = zeros(length(M_arr),length(N0_arr),3);


%tdl initialisation
tdl = nrTDLChannel('MaximumDopplerShift', fd, ...
'DelayProfile','TDL-C','DelaySpread',30e-9, ...
'NumTransmitAntennas',1, 'NumReceiveAntennas',1, ...
'MIMOCorrelation', 'Low','SampleRate',4096*30e3);


%% multiple iterations 
for iter = 1:n_iter
    for ind_M = 1:length(M_arr)
        modulation_params = gen_modulation_params_M_QAM(M_arr(ind_M), Eb); %M-QAM parameters

        %generate random message words in range 0 to M - 1
        message_length = 0.75 * OFDM_params.n_ofdm_symb_per_slot * ...
            OFDM_params.n_RE_per_ofdm_symb; % 0.75 so we can insert pilots in the other 0.25
        message_bits_tx = randi([0,modulation_params.M-1], 1 ,message_length);
        % each of these M numbers(0 to M-1) is equally likely and these
        % M numbers cover all the log2(M) bit numbers (since M is a power
        % of 2). So each of the log2(M) bits is equally

        % start of transmitter logic
        [tx_sequence] = transmitter_logic(message_bits_tx, OFDM_params, modulation_params);
        %transmitted sequence is same for a given M for any SNR
        
        tdl_output = tdl(tx_sequence.').';
        
        for ind_N0 = 1:length(N0_arr)
            N0          = N0_arr(ind_N0);  %N0 value
            noiseStd    = sqrt(N0/2);      %noise standard deviation per each dimension
            Eb_by_N0_dB = 10*log10(Eb/N0); %SNR per bit
            
            % white noise addition
            %tdl output and transmitted sequence are same for every SNR
            %but changes for different M-QAM
            rx_sequence = AWGN_channel(tdl_output, noiseStd);
%             rx_sequence = AWGN_channel(tx_sequence, noiseStd);
            
            % demodulation with uncorrected timing offset (since no timing offset introduced)
            [message_bits_rx] = receiver_logic(rx_sequence, 0, OFDM_params, modulation_params);
            
            %calculating ser and ber by comparing received and transmitted bits
            [ser_e, ber_e, ser_t, ber_t] = calculate_error_rates(...
                message_bits_tx, message_bits_rx, noiseStd, modulation_params);
            results_theoretical(ind_M,ind_N0,:)   = [Eb_by_N0_dB, ser_t, ber_t];
            results_empirical(ind_M,ind_N0,:)     = [Eb_by_N0_dB, ser_e, ber_e];
        end
    end
    average_results_emp(:,:,1) = results_empirical(:,:,1);
    average_results_emp(:,:,2) = average_results_emp(:,:,2) + results_empirical(:,:,2);
    average_results_emp(:,:,3) = average_results_emp(:,:,3) + results_empirical(:,:,3);
end
average_results_emp(:,:,2:3) = average_results_emp(:,:,2:3)/n_iter;

plot_avg_ber_curves(average_results_emp, results_theoretical, M_arr, n_iter, fd);
% saveas(gcf,sprintf('pics2/TDL_%dHz.jpg',fd));
% saveas(gcf,sprintf('pics2/TDL_%dHz.png',fd));
% saveas(gcf,sprintf('pics2/TDL_%dHz.pdf',fd));

%% functions used 
function [] = plot_avg_ber_curves(average_results_emp, results_theoretical, M_arr, n_iter, fd)
figure()
set(gcf,'position',[0,0,800,500])

for ind_M = 1:length(M_arr)
    semilogy(average_results_emp(ind_M,:,1),average_results_emp(ind_M,:,3) ...
        , 'DisplayName', sprintf('M=%dtdl', M_arr(ind_M)));
    hold on;
    semilogy(results_theoretical(ind_M,:,1),results_theoretical(ind_M,:,3) ...
        , 'DisplayName', sprintf('M=%dth', M_arr(ind_M)));
end

title(sprintf('BER vs E_b/N_o (in dB) for M QAM \n tdl channel %d iterations fd = %d Hz',n_iter,fd));
xlabel('10log_{10}(E_b/N_o)');
ylabel('Bit Error Rate(BER)');
legend()%,'Location', 'southwest');
grid on;
ylim([1e-4,1.1]);

end

function symb_with_pilots = add_pilots(message_symb, pilot_value, n_RE)
symb_with_pilots = zeros(1,4*length(message_symb)/3);
% disp(numel(symb_with_pilots));
l1 = 0; l2 = 0;
ind    = 1;
while l2 < length(message_symb)
    if mod(ind,2) ==  1
        symb_with_pilots(l1 + 1: 2: l1 + n_RE) = pilot_value;
        symb_with_pilots(l1 + 2: 2: l1 + n_RE) = message_symb(l2 + 1 : l2 + n_RE/2);
        l2 = l2 + n_RE/2;
    else
        symb_with_pilots(l1 + 1 : l1 + n_RE) = message_symb(l2 + 1 : l2 + n_RE);
        l2 = l2 + n_RE;
    end
    l1  = l1  + n_RE;
    ind = ind + 1;
end
end

function symb_without_pilots = remove_pilots(message_symb, n_RE)
symb_without_pilots = zeros(1,0.75*length(message_symb));
l1 = 0; l2 = 0;
ind    = 1;
while l2 < length(message_symb)
    if mod(ind,2) ==  1
        symb_without_pilots(l1 + 1 : l1 + n_RE/2) = message_symb(l2 + 2: 2: l2 + n_RE);
        l1 = l1 + n_RE/2;
    else
        symb_without_pilots(l1 + 1 : l1 + n_RE) = message_symb(l2 + 1 : l2 + n_RE);
        l1 = l1 + n_RE;
    end
    l2  = l2  + n_RE;
    ind = ind + 1;
end
end

function channel_estimates = get_channel_estimates(resource_grid, ...
                                    n_ofdm_symb_per_slot, n_RE, pilot_value)
    assert( mod(n_ofdm_symb_per_slot,2) == 0)
    assert( mod(n_RE,2) == 0)
    tmp                     = zeros(n_ofdm_symb_per_slot/2,   n_RE);
    channel_estimates       = zeros(n_ofdm_symb_per_slot,   n_RE);
    tmp(:, 1:2:end)   = resource_grid(1:2:end, 1:2:end) ./ pilot_value;
    %for subcarriers without pilots take average of estimates from
    %neighbouring subcarriers with pilots
    tmp(:, 2:2:end-2) = 0.5* ( tmp(:, 1:2:end-3) + tmp(:, 3:2:end-1));
    tmp(:, end)       = tmp(:, end-1);
    
    % taking average of estimates for the same channel over time
%     n_estimates = (1:n_ofdm_symb_per_slot/2).';
%     tmp = cumsum(tmp,1)./n_estimates;    
    channel_estimates(1:2:end)      = tmp;
    channel_estimates(2:2:end)      = tmp;
%     channel_estimates(2:2:end-2) = 0.5* (  channel_estimates(1:2:end-3) ...
%                                          + channel_estimates(3:2:end-1));
%     channel_estimates(end)       = channel_estimates(end-1);

%     channel_estimates(:) = 1;
%     channel_estimates;
end


function output = zero_forcing(channel_estimates, resource_grid)
    output = resource_grid./channel_estimates;
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
pilot_value          =  OFDM_params.pilot_value;

%convert the message bits into complex symbols of M-QAM
symb_indices = bits2symb_map(message_bits_tx+1, :);
symb_coors = ( 2*symb_indices - (root_M+1) )*d;% convert indices to coordinates
message_symb_tx = symb_coors(:,1) + symb_coors(:,2)*1i;% complex msg symbols
message_symb_tx = add_pilots(message_symb_tx, pilot_value, n_RE_per_ofdm_symb);

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
pilot_value           = OFDM_params.pilot_value;

%cp removal
fft_input_arr = zeros(n_ofdm_symb_per_slot, fft_size);
l_used = 0;
for ind = 1:n_ofdm_symb_per_slot
    cp_size = 288 + 64*(ind == 1);
    len_tmp = fft_size + cp_size;
    fft_input_arr(ind,:) = cp_remove(rx_sequence(l_used+1:l_used+len_tmp),cp_size);
    l_used = l_used + len_tmp;
end

%do fft and discard sub carriers where we put zeros 
resource_grid_symb_rx = fft_custom(fft_input_arr, n_RE_per_ofdm_symb, timing_offset);

channel_estimates = get_channel_estimates(resource_grid_symb_rx, ...
                         n_ofdm_symb_per_slot, n_RE_per_ofdm_symb, pilot_value);
resource_grid_symb_rx = zero_forcing(channel_estimates, resource_grid_symb_rx);

%received message symbols from resource grid
message_symb_rx  = reshape(resource_grid_symb_rx.', 1, []);
message_symb_rx  = remove_pilots(message_symb_rx,n_RE_per_ofdm_symb);

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

