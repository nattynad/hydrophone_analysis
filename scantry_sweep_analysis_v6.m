close all
clc

% define the folder, and the filename of the HDF5 file for reading
% dropbox_loc = 'D:\Dropbox';
% dropbox_loc = 'C:\Users\Marc\Dropbox';
folder = 'C:\Users\RKPC\Dropbox\Toronto Team\Calibration Data\600-T1650H825\Scan Data\01-H825_hydrophone_calibration';
fileList = {'600_T1650H825_sweep_825kHz_5_cycles_01.hdf5', '600_T1650H825_sweep_825kHz_5_cycles_02.hdf5', ...
    '600_T1650H825_sweep_825kHz_5_cycles_03.hdf5', '600_T1650H825_sweep_825kHz_5_cycles_04.hdf5', ...
    '600_T1650H825_sweep_825kHz_5_cycles_05.hdf5'};

save_location = 'C:\Users\RKPC\Dropbox\Toronto Team\Calibration Data\323-T1500H750\other (matlab sweep)';

% Extraxt the pressure measuerments from the source
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

hyd_sens_array = [];
%%
for n = 1:length(fileList)    
    % specify the full path to the HDF5 data file
    h5_filename = fullfile(folder, fileList{n});
    
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
    column_raw_press_Pa = size(raw_press_Pa, 2);
    
    % the only thing that changes between scans is the raw_press_Pa variable.
    % The files just append. Therefore just rewrite that var at pos 1-5.
    
    if column_raw_press_Pa > len_input_mV
        raw_press_Pa(:, 1:len_input_mV) = raw_press_Pa(:, end-(len_input_mV-1):end);
    end
    
    metadata = h5read(h5_filename, '/Scan/Scan metadata');
    
    % get sampling period
    my_string = metadata{3};
    % length(my_string);
    outputs = strsplit(my_string);
    sampling_period = str2double(outputs{4});
    
    sampling_frequency = 1/(sampling_period*1e-9);
    
    %gets the trigger delay (how much we delayed the graph)
    my_string = metadata{3};
    outputs = strsplit(my_string);
    trig_delay = str2double(outputs{11});
    
    % gets the hydrophone sensitivity multiplier
    my_string = metadata{4};
    % length(my_string);
    outputs = strsplit(my_string);
    hyd_sens_multiplier = str2double(outputs{3});
    
    % converts the raw output voltage back from pressure
    output_mv = zeros(len_raw_press_Pa, len_input_mV);
    Y = zeros(len_raw_press_Pa, len_input_mV);
    P1_mat = zeros(floor(len_raw_press_Pa/2+1), len_input_mV);
    P2_mat = zeros(len_raw_press_Pa, len_input_mV);
    
    for i = 1:len_input_mV
        output_mv(:, i) = raw_press_Pa(:, i).*hyd_sens_multiplier;
        Y(:, i) = fft(output_mv(:, i));
        % creates the positive-half interval for the fft graph to get the full amplitude
        P2 = abs(Y(:, i)/len_raw_press_Pa);
        P1 = P2(floor(1:len_raw_press_Pa/2+1));
        P1(2:end-1) = 2*P1(2:end-1);
        P1_mat(:, i) = P1;
        P2_mat(:, i) = P2;
    end
    
    % this figure graphs the input voltage going into the transducer (source), and its
    % correlated output pressure (seen by the hydrophone). Its abs min and max because (+) pressure is
    % compression, while (-) pressure is rarefaction
    % 1111 1 1 1 1  1  1 1 1 1 1111 1 1 1 1  1  1  1 1 1 1111
    % a non-linear trend is observed following a certain voltage.
    % note this is not what we are using to observe hydrophone sensitivity,
    % just to observe the effect of ovltage increase on pressure.
    % figure;
    % hold on;
    % plot(input_mV, abs(min_mV), '^');
    % plot(input_mV, abs(max_mV), 'o');
    % grid on;
    % xlabel('Input, mV');
    % ylabel('Output, Pa');
    % legend('Negative Pressure (Min)', 'Positive Pressure (Max)');
    % title('Sweep %d', n);
    % hold off;
    %% 
    ts = zeros(len_raw_press_Pa, len_input_mV);
    for i = 1:len_input_mV
        ts(:, i) = ((0:length(raw_press_Pa(:, i)) - 1) .* sampling_period + (trig_delay*1e-9))';
        % this plots the time characetristic of the transducers voltage amplitude at the 
        % input voltages 10mV, 20 mV, 30 mV, 40 mV, and 50 mV.
        % It is multiplied by the hyrdophone sensitivity to yield the voltage value
        % of the hydrophone instead of the pressure it observes. It then stores
        % it in a matrix
        % figure;
        % hold on;
        % plot(ts(:, i), raw_press_Pa(:, i).*hyd_sens_multiplier);
        % grid on;
        % xlabel('Time (ns) ');
        % ylabel('Output (mV)');
        % titleString = sprintf('Time Profile at %d mV, Sweep %d', i*10, n);
        % title(titleString);
        % hold off;
    end
    
    %% 
        % the first half gets the resolution/freq. spacing of the fft. The second
        % part yields the length of the positive field. multiply that and you have
        % your corresponding frequencies in Hz.
        f = sampling_frequency/len_raw_press_Pa*(0:(len_raw_press_Pa/2));
    for i = 1:len_input_mV
        % plots the FFT. Here the ampltidue is important (voltage seen at the
        % hydrophone).
        % figure;
        % hold on;
        % plot(f, P1_mat(:, i));
        % grid on;
        % ylabel('fft Amplitude');
        % xlabel('f (Hz)');
        % titleString = sprintf('FFT at %d mV, Sweep %d', i*10, n);
        % title(titleString);
        % hold off;
    end
    
    %%
    max_volt = [];
    
    % this subtracts a frequency by 825kHz and finds the minimum difference. the
    % minimum difference means it it the number closest to 825kHz
    [~, idx] = min(abs(f - 825000));
    % closest_y_val = P1_mat(idx,:);
    
    % collecting the voltages that occur in the hydrophone at the resonant f of the transducers
    for i = 1:length(input_mV)
        max_volt(i, 1)= P1_mat(idx, i);
    end;
    
    % figure;
    % plot(input_mV, max_volt, 'o');
    % xlabel('Input Voltage to the Transducer, mV');
    % ylabel('Output seen by Hydrophone, V');

    %%    
    % figure;
    % hold on;
    % plot(input_mV_source, neg_pressure_source, 'o');
    % % plot(input_mV, abs(max_mV), 'o');
    % hold off;
    % grid on;
    % xlabel('Input, mV');
    % ylabel('Peak Negative Pressure, Pa');
    
    input_press_source = zeros(len_input_mV, 1);
    
    k = 100;
    for i = 1:len_input_mV
        input_press_source(i) = (neg_pressure_source(find(input_mV_source==k)))/10;
        k = k + 100;
    end
    
    %%
    % now we want to plot the hydrophone voltage over the transducer source
    % pressure it exerted.
    
    hyd_sens = polyfit(input_press_source, max_volt, 1);
    xfit = linspace(min(input_press_source), max(input_press_source), 10000);
    yfit = polyval(hyd_sens, xfit);
    hyd_sens = hyd_sens(1)*1e+9;
    
    figure;
    hold on
    plot(input_press_source, max_volt, 'o');
    %plot(input_press_source, max_volt, '--');
    Hfit = plot(xfit, yfit, '-', "MarkerFaceColor",[0.8500 0.3250 0.0980]);
    xlabel('Pressure of Source, Pa');
    ylabel('Voltage seen on Hydrophone, V');
    legendString = sprintf('Sensitivity = %.1f mV/MPa', hyd_sens);
    legend(Hfit, legendString);
    titleString = sprintf('Sensitivity, Sweep %d', n);
    title(titleString);
    hold off
    
    fprintf('Sweep %d, Hydrophone Sensitivity = %0.1f mV/MPa\n', n, hyd_sens);

    hyd_sens_array = [hyd_sens_array; hyd_sens];
end

[std_dev_deci, final_hyd_sens] = std(hyd_sens_array);
std_dev = round(std_dev_deci);
if std_dev == 0
    std_dev = 1;
end
fprintf('\nFinal Hydrophone Sensitivity = %0.0f ± %d mV/MPa\n\n', final_hyd_sens, std_dev);

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

