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
hyd_sens_multiplier = outputs{3};
hyd_sens_multiplier = str2double(hyd_sens_multiplier);

% converts the raw output pressure to voltage
output_mv = zeros(len_raw_press_Pa, 5);
Y = zeros(len_raw_press_Pa, 5);
P1_mat = zeros(floor(len_raw_press_Pa/2+1), 5);
P2_mat = zeros(len_raw_press_Pa, 5);
for i = 1:5
    output_mv(:, i) = raw_press_Pa(:, i).*hyd_sens_multiplier;
    Y(:, i) = fft(output_mv(:, i));
    % creates the positive-half interval for the fft graph to get the full amplitude
    P2 = abs(Y(:, i)/len_raw_press_Pa);
    P1 = P2(floor(1:len_raw_press_Pa/2+1));
    P1(2:end-1) = 2*P1(2:end-1);
    P1_mat(:, i) = P1;
    P2_mat(:, i) = P2;
end

% get sampling period
my_string = metadata{3};
% length(my_string);
outputs = strsplit(my_string);
sampling_period = outputs{4};
sampling_period = str2double(sampling_period);

sampling_frequency = 1/(sampling_period*1e-9);


% this figure graphs the input voltage going into the transducer (source), and its
% correlated output pressure (seen by the hydrophone). Its abs min and max because (+) pressure is
% compression, while (-) pressure is rarefaction
% 1111 1 1 1 1  1  1 1 1 1 1111 1 1 1 1  1  1  1 1 1 1111
% a non-linear trend is observed following a certain voltage.
% note this is not what we are using to observe hydrophone sensitivity,
% just to observe the effect of ovltage increase on pressure.
figure;
hold on;
plot(input_mV, abs(min_mV), '^');
plot(input_mV, abs(max_mV), 'o');
grid on;
xlabel('Input, mV');
ylabel('Output, Pa');
hold off;
%% 
ts = zeros(len_raw_press_Pa, 5);
for i = 1:5
    ts(:, i) = ((0:length(raw_press_Pa(:, i)) - 1) .* sampling_period + (35008e-9))';
    % this plots the time characetristic of the hydrophone voltage amplitude at the 
    % input voltages 10mV, 20 mV, 30 mV, 40 mV, and 50 mV to the transducer.
    % It is multiplied by the hyrdophone sensitivity to yield the voltage value
    % of the hydrophone instead of the pressure it observes. It then stores
    % it in a matrix.
    figure;
    hold on;
    plot(ts(:, i), raw_press_Pa(:, i).*hyd_sens_multiplier);
    grid on;
    xlabel('Time (ns) ');
    ylabel('Output (mV)');
    title(['Time Profile at ', num2str(i*10), ' mV']);
    hold off;
end

%% 
    % the first half gets the resolution/freq. spacing of the fft. The second
    % part yields the length of the positive field. multiply that and you have
    % your corresponding frequencies in Hz.
    f = sampling_frequency/len_raw_press_Pa*(0:(len_raw_press_Pa/2));
for i = 1:5
    % This plots the FFT. Here the ampltidue is important (voltage seen at the
    % hydrophone).
    figure;
    hold on;
    plot(f, P1_mat(:, i));
    grid on;
    ylabel('fft Amplitude');
    xlabel('f (Hz)');
    title(['FFT at ', num2str(i*10), ' mV']);
    hold off;
end
%%
max_volt = zeros(5, 1);
resonant_freq = zeros(5, 1);
index_resonance = zeros(5, 1);

% collecting the voltages that occur in the hydrophone at the resonant f of the transducer
for i = 1:5
    [max_volt(i, 1), index_resonance(i, 1)] = max(P1_mat(:, i));
    resonant_freq(i, 1) = f(index_resonance(i, 1));
end;

% plots the voltages seen by the hydrophone over the input voltages to the
% transducer.
figure;
plot(input_mV, max_volt, 'o');
xlabel('Input Voltage to the Transducer, mV');
ylabel('Output seen by Hydrophone, V');

%%
% Now we want to extraxt the pressure measuerments from the source
% transducer to convert the voltage input to pressure
% loading the fils into scantry_sweep_analysis_v1.m, the sweep does not
% start until 50 mV. Because the relation is relatively linear, we can use
% 100 mV and divide the corresponding pressure by 10 to get the pressure at
% 10 mV.

% get the source transducer data

folder_source = 'C:\Users\RKPC\Dropbox\Toronto Team\Calibration Data\01-H825\Scan Data';
filename_source = '01_TH825_sweep_825kHz_05.hdf5';

save_location2 = 'C:\Users\Marc\Dropbox\fus_instruments\marc\calibration_data\Hydrophone Scan\524_T1570H750_QUEENS';


% folder = 'fus_instruments\marc\sites\childrens_copley_rk50'
% filename = '508_T1490H750_sweep_1.hdf5';

% hyd_sens = 5.819E-008;  % V/Pa

% specify the full path to the HDF5 data file
h5_filename_source = fullfile(folder_source, filename_source);

h5_info_source = h5info(h5_filename_source);

% load the first burst to understand data lengths
input_mV_source = h5read(h5_filename_source, '/Scan/Input voltage amplitude (mV)');
len_input_mV_source = length(input_mV_source);

min_mV_source = h5read(h5_filename_source, '/Scan/Min output pressure (Pa)');
len_min_mV_source = length(min_mV_source);

metadata_source = h5read(h5_filename_source, '/Scan/Scan metadata');

neg_pressure_source = abs(min_mV_source);

figure;
hold on;
plot(input_mV_source, neg_pressure_source, 'o');
% plot(input_mV, abs(max_mV), 'o');
hold off;
grid on;
xlabel('Input, mV');
ylabel('Peak Negative Pressure, Pa');

input_press_source = zeros(5, 1);
input_press_source(1, 1) = (neg_pressure_source(find(input_mV_source==100)))/10;
input_press_source(2, 1) = (neg_pressure_source(find(input_mV_source==200)))/10;
input_press_source(3, 1) = (neg_pressure_source(find(input_mV_source==300)))/10;
input_press_source(4, 1) = (neg_pressure_source(find(input_mV_source==400)))/10;
input_press_source(5, 1) = (neg_pressure_source(find(input_mV_source==500)))/10;

%%
% now we want to plot the hydrophone voltage over the transducer source
% pressure it exerted.

figure;
hold on
plot(input_press_source, max_volt, 'o');
plot(input_press_source, max_volt, '--');
xlabel('Pressure of Source, Pa');
ylabel('Voltage seen on Hydrophone, V')

hyd_sens = polyfit(input_press_source, max_volt, 1);
hyd_sens = hyd_sens(1)*1e+9;

% hydrophone sensitivity
fprintf('Hydrophone Sensitivity = %0.4f mV/MPa\n', hyd_sens);

%% DONT WORRY ABOUT THIS FOR NOW

% % neg_pressure = abs(min_mV) .* 1e-3 ./ hyd_sens .* 1e-6;  % MPa
% 
% % neg_pressure = abs(min_mV) .* 1e-3 .* 1e-6;
% 
% neg_pressure = abs(min_mV) .* 1e-6;
% 
% figure;
% hold on;
% plot(input_mV, neg_pressure, 'o');
% % plot(input_mV, abs(max_mV), 'o');
% hold off;
% grid on;
% xlabel('Input, mV');
% ylabel('Peak Negative Pressure, MPa');
% 
% x = input_mV;           % Create Data
% y = neg_pressure;           % Create Data
% B = x(:) \ y(:);                  % Regression Of Line Through Origin
% xfit = [0, 2000];
% yfit = xfit(:) * B;                  % Calculate Fitted Line
% figure;
% plot(x, y, 'bo')                % Plot Data
% hold on
% plot(xfit, yfit, '-r')             % Plot Fitted Regression Line
% hold off
% grid
% axis([0  max(xfit)    0  max(yfit)])
% 
% yCalc1 = B * x;
% 
% % calculate R-squared
% R2 = 1 - sum((y - yCalc1).^2) / sum((y - mean(y)).^2);
% 
% % fprintf('vol2press = %0.4f MPa/Vpp\n', B * 1e3);
% 
% 
% % de-rate the voltage for different electronics box by -6dB:
% 
% gain_EB50_cal = 10 * log10(1.768 / 0.022);  % S/N 2183, NP-977 
% gain_EB50_sys = 10 * log10(1.670 / 0.023);  % S/N 4063
% 
% gain_diff = gain_EB50_cal - gain_EB50_sys;
% 
% x = input_mV * 10^(gain_diff / 20);           % Create Data
% y = neg_pressure;           % Create Data
% B = x(:) \ y(:);                  % Regression Of Line Through Origin
% xfit = [0, 1500 * 10^(gain_diff / 20)];
% yfit = xfit(:) * B;                  % Calculate Fitted Line
% figure;
% plot(x, y, 'bo')                % Plot Data
% hold on
% plot(xfit, yfit, '-r')             % Plot Fitted Regression Line
% hold off
% grid
% axis([0  max(xfit)    0  max(yfit)])
% 
% yCalc1 = B * x;
% 
% % calculate R-squared
% R2 = 1 - sum((y - yCalc1).^2) / sum((y - mean(y)).^2);
% 
% % fprintf('vol2press = %0.4f MPa/Vpp\n', B * 10);
% fprintf('vol2press = %0.4f MPa/Vpp\n', B * 1e3);
% 
% % calculate electrical power and voltage across the transducer
% % EB-50 S/N 2183
% 
% gain_EB50 = gain_EB50_cal;  % at 1.0 MHz
% 
% v_in = input_mV * 1e-3;  % in Vpp
% v_out = v_in .* 10.^(gain_EB50 ./ 20.0);
% 
% data_mtx = zeros(length(neg_pressure), 3);
% 
% data_mtx(:, 1) = neg_pressure;
% data_mtx(:, 2) = v_out;
% data_mtx(:, 3) = v_out .^ 2 ./ 8.0 / 50.0;
% 
% % writematrix(data_mtx,'M_tab.txt','Delimiter','comma');
% % type 'M_tab.txt'
% 
% % header_arr = {'Peak Negative Pressure (MPa)', ...
% %     'Voltage Across the Transducer (Vpp)', ...
% %     'Electrical Power (W)'
% % };
% % 
% % save_filename = 'sweep_524-T1570H750_1570kHz.txt';
% % fid = fopen(fullfile(save_location, save_filename), 'wt');
% % 
% % [n_rows, n_cols] = size(data_mtx);
% % 
% % for col_hdr = 1 : n_cols
% %     fprintf(fid, '%s', header_arr{col_hdr});
% %     if col_hdr ~= n_cols
% %         fprintf(fid, ', ');
% %     end
% % end
% % fprintf(fid, '\n');
% % 
% % save_data_filename = 'sweep_524-T1570H750_1570kHz_DATA.txt';
% % fid_data = fopen(fullfile(save_location, save_data_filename), 'wt');
% % 
% % for row = 1 : n_rows
% %     for col = 1 : n_cols        
% %         fprintf(fid, '%0.3f', data_mtx(row, col));   
% %         fprintf(fid_data, '%0.3f', data_mtx(row, col));   
% %         if col ~= n_cols
% %             fprintf(fid, ', ');
% %             fprintf(fid_data, ', ');
% %         end        
% %     end
% %     fprintf(fid, '\n');
% %     fprintf(fid_data, '\n');
% % end
% % fclose(fid);
% % fclose(fid_data);

