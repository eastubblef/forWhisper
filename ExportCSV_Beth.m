%% This script saves the csv file
%from Junchol 9.15.17
%can also just do this within jrclust and then copy over to RAID2
%To access RAID2 from laptop, just drag the folder into the terminal command window

% addpath(genpath('C:\Users\labadmin\Desktop\jrclust'));                % jrclust folder
% cd('Z:\parkj\NeuralData\WR23_acc_dms_090117')                         % RAID directory
% jrc3 exportcsv WR23_acc_dms_090117_g0_t0.nidq_hh2_2probes_park.prm    % spike timing, cluster number and max. site locations are saved in a tabular format
% csvFileName = 'WR23_acc_dms_090117_g0_t0.nidq_hh2_2probes_park';      % name of the csv file to be read
% csv = csvread(csvFileName); % read the csv file (1st col: spike timing, 2nd col: cluster id, 3rd col: max site)

%% For my files:
% addpath(genpath('C:\Users\labadmin\Desktop\jrclust'));            % jrclust folder
% cd('Z:\Shared\Beth\ImecProcessedData')                            % RAID directory
% jrc3 exportcsv Vgateight_17_04_27_g0_t0.imec.ap_imec3_opt3.prm    % spike timing, cluster number and max. site locations are saved in a tabular format

% csvFileName = 'Vgateight_17_04_27_g0_t0.imec.ap_imec3_opt3.csv';    % name of the csv file to be read
csvFileName = 'Vgateight_170428_g0_t0.imec.ap_imec3_opt3.csv';        % name of the csv file to be read

csv = csvread(csvFileName);                                           % read the csv file (1st col: spike timing, 2nd col: cluster id, 3rd col: max site)

