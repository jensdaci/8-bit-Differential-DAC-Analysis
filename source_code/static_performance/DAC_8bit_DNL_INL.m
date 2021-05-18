% Course Name:  Analog-Digital Interfaces
% Course Code:  ELEN E6316
% Assignment:   PROJECT MEASUREMENTS

clear;
clc;
format long;
N = 8; % 8-bit DAC
V_FS = 1.5; % Full Scale Voltage

filename = 'ss_8-bit_DAC_DNL_INL.csv';
DAC_level = csvread(filename);

% Plotting the DAC Code vs. Analog Output Value
figure(1);
a = DAC_level;
d = (0:1:255);
plot(d, a)
xlim([0 257]);
ylim([-0.760 0.760]);
title("Digital Input Code vs. Analog Output Value")
xlabel("Digital Input Code");
ylabel("Analog Output Value (V)");

% Finding the effective LSB size
V_LSB_actual = (DAC_level(256) - DAC_level(1)) / (2^N - 1);
fprintf('Effective LSB = %.5f mV \n\n', V_LSB_actual * 10^3);

% Finding the Offset Error
V_out_code0 = DAC_level(1);
V_LSB_ideal = V_FS / 2^N;
Offset_Error = V_out_code0 / V_LSB_ideal;
fprintf('Offset Error = %.5f LSBs \n', Offset_Error);

% Finding the Full Scale Error
V_out_all_1_actual = DAC_level(256);
V_out_all_1_ideal = (2^N-1) * V_LSB_ideal;  % 511*VLSB
FullScale_Error = (V_out_all_1_actual - V_out_all_1_ideal) / V_LSB_ideal;
fprintf('Full Scale Error = %.5f LSBs \n\n', FullScale_Error);

% Calculating the DNL and INL 
DNL = zeros(256, 1);
DNL(1) = 0; 

INL = zeros(257, 1);
INL(1) = 0;

for i=2:256
    DNL(i) = ((DAC_level(i) - DAC_level(i-1)) / V_LSB_actual) - 1;
    INL(i+1) = sum(DNL);
end    

% --- DNL ---
% Plotting the DNL vs. DAC Code
figure(2);
plot((0:1:255), DNL)
xlim([0 258]);
title("DNL vs. DAC Code (+0.00098 / -0.00101, std dev 0.00021)");
xlabel("Digital Input Code");
ylabel("DNL");

% Standard deviation of the DNL
DNL_std_dev = std(DNL);
fprintf('Standard Deviation (DNL) = %.5f \n', DNL_std_dev);

% Minimum and Maximum values of DNL 
DNL_Min = min(DNL);
fprintf('Minimum DNL = %.5f LSBs \n', DNL_Min);
DNL_Max = max(DNL);
fprintf('Maximum DNL = %.5f LSBs \n', DNL_Max);
fprintf('The DAC is monotonic. Minimum DNL > -1 LSB. \n\n')

% --- INL ---
% Plotting the INL vs. DAC Code
figure(3);
plot((0:1:256), INL)
xlim([0 258]);
title("INL vs. DAC Code (+0.00706 / -0.01193, std dev 0.00651)");
xlabel("Digital Input Code");
ylabel("INL");

% Standard deviation of the INL
INL_std_dev = std(INL);
fprintf('Standard Deviation (INL) = %.5f \n', INL_std_dev);

% Minimum and Maximum values of INL
INL_Min = min(INL);
fprintf('Minimum INL = %.5f LSBs \n', INL_Min);
INL_Max = max(INL);
fprintf('Maximum INL = %.5f LSBs \n', INL_Max);


