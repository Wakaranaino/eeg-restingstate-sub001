% ------------------------------------------------------------------------------------------
% Version: v2.1 | Last updated: June 9, 2025
% Updates:
% • Added interpolation of removed channels (based on original montage) using spherical method.
% • Logged interpolated channel labels per subject to .mat file for transparency.
% • Updated QC log to include both count and list of interpolated channels.
% Purpose:
% To ensure consistent channel count across subjects for group-level topographic analysis.
% ------------------------------------------------------------------------------------------

% EEG Analysis Script for Subject 001
% Manual-style processing turned into script
% Pipeline: Load → ICA → Component Rejection → PSD → Bandpower → Topomap
% Author: Vermut Gao
% Date: Jun 1st, 2025

% =================== MULTI-SUBJECT LOOP SETUP =====================
clear; close all; clc

subject_list = {'sub-004','sub-001','sub-002','sub-003','sub-037','sub-038','sub-039','sub-040','sub-066','sub-067','sub-068','sub-069'};  % Add more subject IDs here
base_input_dir = '/Users/bohago/Desktop/neuro analysis/EEG resting-state_Jun 2025/';
output_dir = 'outputs';  % Common output directory for all subjects

% Create output folder if needed
if ~isfolder(output_dir)
    status = mkdir(output_dir);
    if ~status
        error('Failed to create output folder. Check permissions or path.');
    end
end

% Set band definitions once for all subjects
band_defs.theta = [4 7];
band_defs.alpha = [8 13];
band_defs.beta  = [13 30];
band_names = fieldnames(band_defs);

% Master QC file (initialize only once)
qc_file = fullfile(output_dir, 'qc_summary_all_subjects.csv');
if ~isfile(qc_file)
    fid = fopen(qc_file, 'w');
    if fid == -1
        error('Failed to open QC summary file for writing.');
    end
    fprintf(fid, 'Subject,ICsRejected,ThetaMean,AlphaMean,BetaMean\n');
    fclose(fid);
end

% Launch EEGLAB once
eeglab;

% ======================= MAIN SUBJECT LOOP ========================
for s = 1:length(subject_list)
    subject_id = subject_list{s};
    disp(['--- Processing: ' subject_id ' ---']);
    output_prefix = [subject_id '_v2.1'];
    input_path = fullfile(base_input_dir, subject_id, 'eeg');
    input_filename = [subject_id '_task-eyesclosed_eeg.set'];

    % Create per-subject output directory
    subject_output_dir = fullfile(output_dir, subject_id);
    if ~isfolder(subject_output_dir)
        mkdir(subject_output_dir);
    end

% Skip preprocessing if cleaned set already exists
load_cleaned_set = exist(fullfile(subject_output_dir, [output_prefix '_cleaned_final.set']), 'file');

if load_cleaned_set
    EEG = pop_loadset('filename', [output_prefix '_cleaned_final.set'], 'filepath', subject_output_dir);
    artifactICs = EEG.etc.artifactICs;
    [ALLEEG, EEG, CURRENTSET] = eeg_store([], EEG, 0);
else
    % Load original dataset
    EEG = pop_loadset('filename', input_filename, 'filepath', input_path);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);

    % Save original channel locations for interpolation
    original_chanlocs = EEG.chanlocs;  

    % Automated artifact rejection using `clean_rawdata()` (ASR 4 SD, soft 0.5 Hz high-pass)
    EEG = clean_rawdata(EEG, 5, [0.25 0.75], 0.85, 4, 4, 0.2);

    % Apply low-pass filter at 40 Hz (high-pass already handled by clean_rawdata)
    EEG = pop_eegfiltnew(EEG, 'hicutoff', 40);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

    % Re-reference to average
    EEG = pop_reref(EEG, []);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

    % Detect which channels were removed
    % BEFORE interpolation
    remaining_labels = {EEG.chanlocs.labels};
    interp_flags = ~ismember({original_chanlocs.labels}, remaining_labels);
    interpolated_channels = {original_chanlocs(interp_flags).labels};

    % Save log BEFORE interpolation
    interp_log.subject_id = subject_id;
    interp_log.channel_labels = {original_chanlocs.labels};
    interp_log.interpolated_channels = interpolated_channels;
    interp_log.interp_flags = interp_flags;
    save(fullfile(subject_output_dir, [output_prefix '_interpolation_log.mat']), '-struct', 'interp_log');

    % Interpolate missing channels for group analysis
    EEG = pop_interp(EEG, original_chanlocs, 'spherical');

    % Run ICA (extended Infomax)
    EEG = pop_runica(EEG, 'extended', 1);
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

    % Label components using ICLabel
    EEG = pop_iclabel(EEG, 'default');
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

    % Reject components with Brain probability < 0.5 (automated ICLabel threshold)
    artifactICs = find(EEG.etc.ic_classification.ICLabel.classifications(:,1) < 0.5); 
    EEG.etc.artifactICs = artifactICs;  % <-- Save to EEG struct
    EEG = pop_subcomp(EEG, artifactICs, 0);

    % Save final dataset
    EEG = pop_saveset(EEG, 'filename', [output_prefix '_cleaned_final.set'], 'filepath', subject_output_dir, 'savemode', 'onefile');
end


% Power Spectral Density (PSD) Computation
[psdAll, freqs] = spectopo(EEG.data, 0, EEG.srate, 'plot', 'off', 'freqrange', [2 25]);
writematrix(psdAll, fullfile(subject_output_dir, [output_prefix '_psd_matrix.csv']));

% PSD scalp maps
close all;
specdata = pop_spectopo(EEG, 1, [0 EEG.pnts], 'EEG', 'freqrange',[2 25], 'electrodes','off', 'freqs',[6 10 22]);
saveas(gcf, fullfile(subject_output_dir, [output_prefix '_psd_scalp.png']));

% Compute bandpower for Theta, Alpha, Beta bands across all EEG channels

fs = EEG.srate;
bandpower_struct = struct();

for i = 1:length(band_names)
    band = band_names{i};
    range = band_defs.(band);
    power_vals = zeros(1, EEG.nbchan);

    for ch = 1:EEG.nbchan
        data = EEG.data(ch, :);
        [pxx, f] = pwelch(data, hamming(1024), [], [], fs);
        power_vals(ch) = bandpower(pxx, f, range, 'psd');
    end

    bandpower_struct.(band) = power_vals;
end

% Create Bandpower Table Dynamically
channel_labels = {EEG.chanlocs.labels}';
T = table(channel_labels, 'VariableNames', {'Channel'});

for i = 1:length(band_names)
    band = band_names{i};
    col_data = bandpower_struct.(band)';
    col_name = [upper(band(1)) band(2:end) 'Power'];  % e.g., ThetaPower
    T.(col_name) = col_data;
end

% Write to CSV
writetable(T, fullfile(subject_output_dir, [output_prefix '_bandpower_matrix.csv']));

disp(['Subject: ' subject_id]);
disp(['ICs rejected: ' num2str(length(artifactICs))]);

% =====QC Logging (Console)=======
for i = 1:length(band_names)
    band = band_names{i};
    band_mean = mean(bandpower_struct.(band));
    disp([upper(band(1)) band(2:end) ' Power Mean: ' num2str(band_mean)]);
end

% Normalize power for comparison across subjects
total_power = mean(bandpower_struct.theta + bandpower_struct.alpha + bandpower_struct.beta);
norm_struct.theta = mean(bandpower_struct.theta) / total_power;
norm_struct.alpha = mean(bandpower_struct.alpha) / total_power;
norm_struct.beta  = mean(bandpower_struct.beta)  / total_power;

disp(['Normalized Alpha: ' num2str(norm_struct.alpha)]);

save(fullfile(subject_output_dir, [output_prefix '_summary.mat']), 'bandpower_struct', 'norm_struct', 'band_defs', 'psdAll', 'freqs', 'artifactICs');

log_file = fullfile(subject_output_dir, [output_prefix '_qc_log.txt']);
fid = fopen(log_file, 'w');
if fid == -1
    error('Could not open QC log file for writing.');
end

% ----- Load interpolated channels log if exists -----
interp_log_path = fullfile(subject_output_dir, [output_prefix '_interpolation_log.mat']);
if isfile(interp_log_path)
    interp_data = load(interp_log_path);  % loads as struct
    interpolated_channels = interp_data.interpolated_channels;
else
    interpolated_channels = {};
end

fprintf(fid, 'Subject: %s\n', subject_id);
fprintf(fid, 'ICs rejected: %d\n', length(artifactICs));
fprintf(fid, 'Interpolated Channels: %d\n', length(interpolated_channels));

if ~isempty(interpolated_channels)
    fprintf(fid, 'Channels: %s\n', strjoin(interpolated_channels, ', '));
end

for i = 1:length(band_names)
    band = band_names{i};
    band_mean = mean(bandpower_struct.(band));
    fprintf(fid, '%s Power Mean: %.4f\n', [upper(band(1)) band(2:end)], band_mean);
end

fclose(fid);

% Save QC-relevant data in .mat format for later batch/group analysis
save(fullfile(subject_output_dir, [output_prefix '_summary.mat']), 'bandpower_struct', 'band_defs', 'psdAll', 'freqs', 'artifactICs');

% ===== Append QC to master log file =====
existingQC = readtable(qc_file);
if ~any(strcmp(existingQC.Subject, subject_id))
    fid = fopen(qc_file, 'a');
    fprintf(fid, '%s,%d,%.4f,%.4f,%.4f\n', subject_id, length(artifactICs), norm_struct.theta, norm_struct.alpha, norm_struct.beta);
    fclose(fid);
else
    disp(['✔ QC already logged for ' subject_id ', skipping append.']);
end

% Plot topographic maps for each band
close all;
figure('Position', [100, 100, 1200, 400]);
for i = 1:length(band_names)
    subplot(1, length(band_names), i);
    topoplot(bandpower_struct.(band_names{i}), EEG.chanlocs, 'maplimits', 'maxmin', 'electrodes', 'on');
    title([upper(band_names{i}(1)) band_names{i}(2:end) ' Power (' num2str(band_defs.(band_names{i})(1)) '–' num2str(band_defs.(band_names{i})(2)) ' Hz)']);
    colorbar;
end
saveas(gcf, fullfile(subject_output_dir, [output_prefix '_bandpower.png']));

% Combined figure with spectopo
close all;
figure('Position', [100, 100, 1200, 800]);

% Subplot 1: Full PSD using EEGLAB
subplot(2,2,1);
pop_spectopo(EEG, 1, [0 EEG.pnts], 'EEG', 'freqrange', [2 25], 'electrodes', 'off');
title('PSD (2–25 Hz)');
colorbar;

% Subplots 2–4: Bandpower Topomaps (Theta, Alpha, Beta)
for i = 1:length(band_names)
    subplot(2,2,i+1);  % i+1 to start from subplot 2
    topoplot(bandpower_struct.(band_names{i}), EEG.chanlocs, ...
        'maplimits', 'maxmin', 'electrodes', 'on');
    title([upper(band_names{i}(1)) band_names{i}(2:end) ...
        ' Power (' num2str(band_defs.(band_names{i})(1)) ...
        '–' num2str(band_defs.(band_names{i})(2)) ' Hz)']);
    colorbar;
end

saveas(gcf, fullfile(subject_output_dir, [output_prefix '_topomap_summary.png']));


disp(['Done with ' subject_id]);
end

eeglab redraw;
close all;

pause(1);  % Makes the loop friendlier in logs


% ------------------------------------------------------------------------------------------
% Version: v2.0 | Last updated: June 7, 2025
% Updates:
% • Wrapped single-subject pipeline into a multi-subject loop
% • Modularized per-subject output handling (paths, filenames, folders)
% • Added group-level QC summary logging (CSV + normalized bandpower)
% • Implemented bandpower normalization for inter-subject comparison
% • Included external subject integrity checker (NaN, empty, bad channels)
% • Prepared for group-level stats & visualization 
% ------------------------------------------------------------------------------------------
% ------------------------------------------------------------------------------------------
% Version: v1.3 | Last updated: June 6, 2025
% Updates:
% • Standardized parameter setup for multi-subject conversion
% • Unified output path handling (fullfile, pop_saveset)
% • Added per-subject QC logging (.txt + .mat)
% • Appended master QC summary to .csv for group-level tracking
% • Refactored all bandpower-related code into modular structures:
%   - Defined band_defs and band_names for cleaner frequency management
%   - Used bandpower_struct for per-band storage and plotting
%   - Automated topomap titles and QC summary to adapt to new bands
% • Added load_cleaned_set switch for skipping preprocessing on reruns
% ------------------------------------------------------------------------------------------

% ------------------------------------------------------------------------------------------
% Version: v1.2 | Last updated: June 3, 2025
% Updates:
% • Upgraded artifact cleaning to research-level parameters:
%     clean_rawdata(EEG, 5, [0.25 0.75], 0.85, 4, 4, 0.2)
% • Applied low-pass filtering with pop_eegfiltnew (cutoff: 40 Hz)
% • Automated PSD computation and bandpower topomap generation
% • Used spectopo() for cleaner full-channel PSD matrix export
% ------------------------------------------------------------------------------------------

% ------------------------------------------------------------------------------------------
% Previous Versions:
% v1.1 (June 2, 2025) – Initial automation added: clean_rawdata + ICLabel
% ------------------------------------------------------------------------------------------