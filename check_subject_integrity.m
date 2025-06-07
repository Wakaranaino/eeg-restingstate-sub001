
% ======================= QC INTEGRITY CHECK ============================
clear; clc

qc_file = 'outputs/qc_summary_all_subjects.csv';
summary_dir = 'outputs/';

% Load QC summary
T = readtable(qc_file);

% Check for NaNs or zeros in normalized bandpower
bad_idx = any(isnan(T{:, 3:end}) | T{:, 3:end} == 0, 2);
if any(bad_idx)
    disp('Subjects with NaN or zero bandpower values:');
    disp(T.Subject(bad_idx));
else
    disp('No NaNs or zeros in bandpower values.');
end

% Check if ICs rejected is zero (may imply ICA failed)
zeroICs = T.ICsRejected == 0;
if any(zeroICs)
    disp('Subjects with 0 ICs rejected (check ICA output):');
    disp(T.Subject(zeroICs));
else
    disp('All subjects rejected at least 1 IC.');
end

% Load per-subject .mat to check channel count consistency
bad_chan_subjects = {};
expected_channels = [];

for i = 1:height(T)
    subj_id = T.Subject{i};
    mat_file = fullfile(summary_dir, subj_id, [subj_id '_v2.0_summary.mat']);

    if isfile(mat_file)
        S = load(mat_file);
        chan_count = length(S.bandpower_struct.theta);

        % On first subject, set expected
        if isempty(expected_channels)
            expected_channels = chan_count;
        elseif chan_count ~= expected_channels
            bad_chan_subjects{end+1} = subj_id;
        end
    else
        warning(['Missing summary file for ' subj_id]);
    end
end

if isempty(bad_chan_subjects)
    disp(['All subjects have consistent channel count: ' num2str(expected_channels)]);
else
    disp('Channel count mismatch found in:');
    disp(bad_chan_subjects');
end