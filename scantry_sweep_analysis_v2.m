close all
clc

% define the folder, and the filename of the HDF5 file for reading
% dropbox_loc = 'D:\Dropbox';
% dropbox_loc = 'C:\Users\Marc\Dropbox';
folder = 'C:\Users\RKPC\Dropbox\Toronto Team\Calibration Data\571-T1000H550\Scan Data\01-H825_hydrophone_calibration';
filename = '571_T1000H550_sweep_825kHz_01.hdf5';

save_location = 'C:\Users\RKPC\Dropbox\Toronto Team\Calibration Data\323-T1500H750\other (matlab sweep)';


% folder = 'fus_instruments\marc\sites\childrens_copley_rk50'
% filename = '508_T1490H750_sweep_1.hdf5';

% hyd_sens = 5.819E-008;  % V/Pa

% specify the full path to the HDF5 data file
h5_filename = fullfile(folder, filename);

h5_info = h5info(h5_filename);

% load the first burst to understand data lengths
input_mV = h5read(h5_filename, '/Scan/Input voltage amplitude (mV)');
len_input_mV = length(input_mV);

max_mV = h5read(h5_filename, '/Scan/Max output pressure (Pa)'); %positive axis
len_max_mV = length(max_mV);

min_mV = h5read(h5_filename, '/Scan/Min output pressure (Pa)'); %negative axis
len_min_mV = length(min_mV);

fwd_pwr_V = h5read(h5_filename, '/Scan/Forward power meter waveforms (V)');

rev_pwr_V = h5read(h5_filename, '/Scan/Reverse power meter waveforms (V)');

raw_press_Pa = h5read(h5_filename, '/Scan/Raw pressure waveforms (Pa)');
len_raw_press_Pa = length(raw_press_Pa);

metadata = h5read(h5_filename, '/Scan/Scan metadata');

% gets the hydrophone sensitivity multiplier
my_string = metadata{4};
% length(my_string);
outputs = strsplit(my_string);
hyd_sens = outputs{3};
hyd_sens = str2double(hyd_sens);

% converts the raw output voltage back from pressure
outout_mv = raw_press_Pa(:, 1).*hyd_sens;

% get sampling period
my_string = metadata{3};
% length(my_string);
outputs = strsplit(my_string);
sampling_period = outputs{4};
sampling_period = str2double(sampling_period);

sampling_frequency = 1/(sampling_period*1e-9);

Y = fft(outout_mv);
P2 = abs(Y/len_raw_press_Pa);
P1 = P2(1:len_raw_press_Pa/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure;
hold on;
plot(input_mV, abs(min_mV), '^');
plot(input_mV, abs(max_mV), 'o');
grid on;
xlabel('Input, mV');
ylabel('Output, Pa');
hold off;
%% 
ts = ((0:length(raw_press_Pa(:, 1)) - 1) .* sampling_period + (35008e-9))';

figure;
hold on;
plot(ts, raw_press_Pa(:, 1).*hyd_sens);
grid on;
xlabel('Time, ns');
ylabel('Output, mV');
hold off;
%% 

f = sampling_frequency/len_raw_press_Pa*(0:(len_raw_press_Pa/2));

figure;
hold on;
plot(f, P1);
grid on;
ylabel('fft Amplitude');
xlabel('f, Hz')
hold off;


% neg_pressure = abs(min_mV) .* 1e-3 ./ hyd_sens .* 1e-6;  % MPa

% neg_pressure = abs(min_mV) .* 1e-3 .* 1e-6;

neg_pressure = abs(min_mV) .* 1e-6;

figure;
hold on;
plot(input_mV, neg_pressure, 'o');
% plot(input_mV, abs(max_mV), 'o');
hold off;
grid on;
xlabel('Input, mV');
ylabel('Peak Negative Pressure, MPa');

x = input_mV;           % Create Data
y = neg_pressure;           % Create Data
B = x(:) \ y(:);                  % Regression Of Line Through Origin
xfit = [0, 2000];
yfit = xfit(:) * B;                  % Calculate Fitted Line
figure;
plot(x, y, 'bo')                % Plot Data
hold on
plot(xfit, yfit, '-r')             % Plot Fitted Regression Line
hold off
grid
axis([0  max(xfit)    0  max(yfit)])

yCalc1 = B * x;

% calculate R-squared
R2 = 1 - sum((y - yCalc1).^2) / sum((y - mean(y)).^2);

% fprintf('vol2press = %0.4f MPa/Vpp\n', B * 1e3);


% de-rate the voltage for different electronics box by -6dB:

gain_EB50_cal = 10 * log10(1.768 / 0.022);  % S/N 2183, NP-977 
gain_EB50_sys = 10 * log10(1.670 / 0.023);  % S/N 4063

gain_diff = gain_EB50_cal - gain_EB50_sys;

x = input_mV * 10^(gain_diff / 20);           % Create Data
y = neg_pressure;           % Create Data
B = x(:) \ y(:);                  % Regression Of Line Through Origin
xfit = [0, 1500 * 10^(gain_diff / 20)];
yfit = xfit(:) * B;                  % Calculate Fitted Line
figure;
plot(x, y, 'bo')                % Plot Data
hold on
plot(xfit, yfit, '-r')             % Plot Fitted Regression Line
hold off
grid
axis([0  max(xfit)    0  max(yfit)])

yCalc1 = B * x;

% calculate R-squared
R2 = 1 - sum((y - yCalc1).^2) / sum((y - mean(y)).^2);

% fprintf('vol2press = %0.4f MPa/Vpp\n', B * 10);
fprintf('vol2press = %0.4f MPa/Vpp\n', B * 1e3);

% calculate electrical power and voltage across the transducer
% EB-50 S/N 2183

gain_EB50 = gain_EB50_cal;  % at 1.0 MHz

v_in = input_mV * 1e-3;  % in Vpp
v_out = v_in .* 10.^(gain_EB50 ./ 20.0);

data_mtx = zeros(length(neg_pressure), 3);

data_mtx(:, 1) = neg_pressure;
data_mtx(:, 2) = v_out;
data_mtx(:, 3) = v_out .^ 2 ./ 8.0 / 50.0;

% writematrix(data_mtx,'M_tab.txt','Delimiter','comma');
% type 'M_tab.txt'

% header_arr = {'Peak Negative Pressure (MPa)', ...
%     'Voltage Across the Transducer (Vpp)', ...
%     'Electrical Power (W)'
% };
% 
% save_filename = 'sweep_524-T1570H750_1570kHz.txt';
% fid = fopen(fullfile(save_location, save_filename), 'wt');
% 
% [n_rows, n_cols] = size(data_mtx);
% 
% for col_hdr = 1 : n_cols
%     fprintf(fid, '%s', header_arr{col_hdr});
%     if col_hdr ~= n_cols
%         fprintf(fid, ', ');
%     end
% end
% fprintf(fid, '\n');
% 
% save_data_filename = 'sweep_524-T1570H750_1570kHz_DATA.txt';
% fid_data = fopen(fullfile(save_location, save_data_filename), 'wt');
% 
% for row = 1 : n_rows
%     for col = 1 : n_cols        
%         fprintf(fid, '%0.3f', data_mtx(row, col));   
%         fprintf(fid_data, '%0.3f', data_mtx(row, col));   
%         if col ~= n_cols
%             fprintf(fid, ', ');
%             fprintf(fid_data, ', ');
%         end        
%     end
%     fprintf(fid, '\n');
%     fprintf(fid_data, '\n');
% end
% fclose(fid);
% fclose(fid_data);

