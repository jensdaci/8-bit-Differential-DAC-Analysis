clear; close all;
clc;

%Change to the directory for your spectrum file
%directory = '/Users/jens/Documents/MATLAB/ELEN E6316/';
%file comes from Cadence calculator output table, and has 1 header line
file = 'tt_J_sweep_data.csv';
%file = 'J_vs_SFDR_original.csv';

fs = 5e9;
filename = strcat(file);
fullData = readtable(filename);
freq = table2array(fullData(:,1));
numOfSweepElements = width(fullData(1,2:end));

SFDRs = zeros(numOfSweepElements,1);
SNDR_dBs = zeros(numOfSweepElements,1);
ENOBs = zeros(numOfSweepElements,1);

for k=1:numOfSweepElements
    
    temp0 = table2array(fullData(:,k+1));
    
    %throw away the DC bin.  May need to be careful here, if the signal
    %power is located in bin 2
    temp0(1) = temp0(2);

    %The data read from the file is the DAC spectrum that
    %includes the rolloff due to the Zero-order hold function
    fft_dB20_ZOH = temp0;

    %Zero-order hold attenuation in the first (dc) frequency bin is 0
    ZOH(1) = 0;

    %initialize
    fft_dB20_2(1) = fft_dB20_ZOH(1) - ZOH(1);

    for i = 2:length(freq)
        %compute the attenuation from the Zero-Order Hold at each frequency
        ZOH(i) = 20*log10(abs(sin(pi*freq(i)/fs) / (pi*freq(i)/fs)));
        %remove this attenuation from our output spectrum
        fft_dB20_2(i) = fft_dB20_ZOH(i) - ZOH(i);        
    end

    %transpose the matrix from the CSV file
    fft_dB20 = fft_dB20_2';

    fft_mag = 10.^(fft_dB20/20);
    fft_mag_sq = fft_mag.^2;

    max_fft = min(fft_dB20);

    %find the signal (assumed to be the largest tone)
    for i = 1:length(fft_dB20)
        if (fft_dB20(i) >= max_fft)
            max_fft = fft_dB20(i);
            sig_bin = i;
            max_fft_freq = freq(i);
        end
    end


    %initialize the distortion tone to be the smallest power in any fft bin
    max_fft_dist_tone = min(fft_dB20);

    total_power = 0;

    %Find the 2nd largest tone to compute the SFDR.  We skip over the bin that has the signal
    for i = 1:length(fft_dB20)
        total_power = total_power + fft_mag_sq(i);
        if ( i == sig_bin)  %may need to be careful, if the signal power is located in more than one frequency bin
            % do nothing, we want to skip over this bin
        elseif (fft_dB20(i) >= max_fft_dist_tone)
            max_fft_dist_tone = fft_dB20(i);
            dist_bin = i;
            max_fft_dist_freq = freq(i);
        end
    end

    SFDR = max_fft - max_fft_dist_tone;

    sig_power = fft_mag_sq(sig_bin);

    %Noise and distortion power is found by subtracting the signal power from
    %the total power
    ND = total_power - sig_power;

    SNDR = sig_power / ND;  %signal power divided by all the other power
    SNDR_dB = 10*log10(SNDR);
    ENOB = (SNDR_dB-1.76)/6.02;

    %compute the signal amplitude (zero-to-peak and peak-to-peak)
    sig_amp = 10^(max_fft/20);
    sig_Vpp = 2*sig_amp

    %Report the results
    ENOB
    SNDR_dB
    SFDR
    
    SFDRs(k) = SFDR;
    SNDR_dBs(k) = SNDR_dB;
    ENOBs(k) = ENOB;

end

columnNames = fullData.Properties.VariableNames(2:end);
Js = zeros(numOfSweepElements,1);
JStrings = strings(numOfSweepElements,1);
for j=1:numOfSweepElements
    
    JString = extractBetween(columnNames{j}, 'J','_');
    JStrings(j) = JString{1};
    
    Js(j) = str2double(JString);
end

plot(Js, SFDRs); xlabel('J'); ylabel('SFDR (dB)'); title('J vs. SFDR when M=2048 and fsamp=400M');