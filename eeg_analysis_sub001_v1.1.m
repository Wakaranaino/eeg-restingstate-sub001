clear; close all; clc

% -------------------------------------------------------------------------
% EEG Analysis Script – Subject 001
% Version: v1.1 | Last updated: June 2, 2025
% Updates: Added automated artifact rejection via clean_rawdata() and ICLabel
%          (automatically detects and removes artifacts: blinks, muscle, ECG)
% -------------------------------------------------------------------------


% EEG Analysis Script for Subject 001
% Manual-style processing turned into script
% Pipeline: Load → ICA → Component Rejection → PSD → Bandpower → Topomap
% Author: Vermut Gao
% Date: Jun 1st, 2025

eeglab;

% Load original dataset
EEG = pop_loadset('filename', 'sub-001_task-eyesclosed_eeg.set', 'filepath', '/Users/bohago/Desktop/neuro analysis/EEG resting-state_Jun 2025/sub-001/eeg/');
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);

% Automated Pre-Cleaning (via clean_rawdata)
EEG = clean_rawdata(EEG, 5, -1, 0.85, -1, -1, -1);

% Apply 1–40 Hz bandpass filter
EEG = pop_eegfiltnew(EEG, 'locutoff', 1, 'hicutoff', 40);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

% Re-reference to average
EEG = pop_reref(EEG, []);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

% Run ICA (extended Infomax)
EEG = pop_runica(EEG, 'extended', 1);
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

% Label components using ICLabel
EEG = pop_iclabel(EEG, 'default');
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

% Reject selected components (based on ICLabel visual inspection)
% EEG = pop_subcomp(EEG, [1 2 3 5 10 14 16 18 19], 0);
% [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

% Reject components with Brain probability < 0.5 (automated ICLabel threshold)
artifactICs = find(EEG.etc.ic_classification.ICLabel.classifications(:,1) < 0.5); 
EEG = pop_subcomp(EEG, artifactICs, 0);

% Power Spectral Density (PSD) Computation
specdata = pop_spectopo(EEG, 1, [0 599798], 'EEG', 'freqrange',[2 25], 'electrodes','off', 'freqs',[6 10 22]);
set(gcf, 'Position', [100, 100, 800, 600]);  % ensures scalp maps are visible

% Save the specdata to CSV
writematrix(specdata, 'sub001_psd_matrix.csv');

% Compute bandpower for Theta, Alpha, Beta bands across all EEG channels
% Assumes cleaned EEG dataset is already loaded into EEGLAB as 'EEG'

% Band definitions (Hz)
theta_band = [4 7];
alpha_band = [8 13];
beta_band  = [13 30];

fs = EEG.srate;  % Sampling rate

% Initialize bandpower arrays
theta_power = zeros(1, EEG.nbchan);
alpha_power = zeros(1, EEG.nbchan);
beta_power  = zeros(1, EEG.nbchan);

% Loop through channels
for ch = 1:EEG.nbchan
    data = EEG.data(ch, :);
    [pxx, f] = pwelch(data, hamming(1024), [], [], fs);

    theta_power(ch) = bandpower(pxx, f, theta_band, 'psd');
    alpha_power(ch) = bandpower(pxx, f, alpha_band, 'psd');
    beta_power(ch)  = bandpower(pxx, f, beta_band, 'psd');
end

% Create result table
channel_labels = {EEG.chanlocs.labels};
T = table(channel_labels', theta_power', alpha_power', beta_power', ...
    'VariableNames', {'Channel', 'ThetaPower', 'AlphaPower', 'BetaPower'});
disp(T)

writetable(T, 'sub001_bandpower_matrix.csv');

% Plot topographic maps
figure('Position', [100, 100, 1200, 400]);

subplot(1,3,1);
topoplot(theta_power, EEG.chanlocs, 'maplimits', 'maxmin', 'electrodes', 'on');
title('Theta Power (4–7 Hz)');
colorbar;

subplot(1,3,2);
topoplot(alpha_power, EEG.chanlocs, 'maplimits', 'maxmin', 'electrodes', 'on');
title('Alpha Power (8–13 Hz)');
colorbar;

subplot(1,3,3);
topoplot(beta_power, EEG.chanlocs, 'maplimits', 'maxmin', 'electrodes', 'on');
title('Beta Power (13–30 Hz)');
colorbar;

% Save final dataset
EEG = pop_saveset(EEG, 'filename','sub001_cleaned_final.set');

% Reload final dataset into EEGLAB GUI
EEG = pop_loadset('filename','sub001_cleaned_final.set');
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw;

disp('Saved: sub001_psd_matrix.csv');
disp('Saved: sub001_bandpower_matrix.csv');
