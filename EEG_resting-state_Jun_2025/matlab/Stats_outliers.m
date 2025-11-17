% detect_outliers_iqr.m
% IQR-based outlier scan per band x group, with subject IDs and CSV export.

clear; clc;

infile  = 'outputs/group_merged_summary.csv';
stats_dir = 'outputs/stats';
if ~isfolder(stats_dir), mkdir(stats_dir); end

merged = readtable(infile);

% Columns to check (must exist in merged)
bands  = {'ThetaMean','AlphaMean','BetaMean'};
% Groups are coded 'C','A','F' in your files
groups = {'C','A','F'};
pretty  = struct('C','CN','A','AD','F','FTD');

% Collect rows for a CSV
rows_Subject = {};
rows_Group   = {};
rows_Band    = {};
rows_Value   = [];
rows_Lower   = [];
rows_Upper   = [];
rows_IsOut   = [];

fprintf('--- IQR outlier scan per band × group ---\n');

for b = 1:numel(bands)
    band = bands{b};
    fprintf('\n=== %s ===\n', band);
    
    for gi = 1:numel(groups)
        g = groups{gi};
        mask_g = strcmp(merged.Group, g);
        vals = merged.(band)(mask_g);
        subs = merged.Subject(mask_g);

        if isempty(vals)
            fprintf('%s: 0/0 possible outlier(s)\n', pretty.(g));
            continue;
        end
        
        Q1 = quantile(vals, 0.25);
        Q3 = quantile(vals, 0.75);
        IQR = Q3 - Q1;
        lower = Q1 - 1.5*IQR;
        upper = Q3 + 1.5*IQR;
        out_idx = (vals < lower) | (vals > upper);

        n_out = sum(out_idx);
        fprintf('%s: %d/%d possible outlier(s)\n', pretty.(g), n_out, numel(vals));

        if n_out > 0
            % Print subject IDs
            flagged_subs = subs(out_idx);
            disp(flagged_subs);
        end
        
        % Append to CSV rows (for all subjects so you can review thresholds)
        rows_Subject = [rows_Subject; subs]; %#ok<AGROW>
        rows_Group   = [rows_Group; repmat({pretty.(g)}, numel(subs), 1)]; %#ok<AGROW>
        rows_Band    = [rows_Band; repmat({band}, numel(subs), 1)]; %#ok<AGROW>
        rows_Value   = [rows_Value; vals]; %#ok<AGROW>
        rows_Lower   = [rows_Lower; repmat(lower, numel(subs), 1)]; %#ok<AGROW>
        rows_Upper   = [rows_Upper; repmat(upper, numel(subs), 1)]; %#ok<AGROW>
        rows_IsOut   = [rows_IsOut; out_idx]; %#ok<AGROW>
    end
end

% Build and save the tidy table
T = table(rows_Subject, rows_Group, rows_Band, rows_Value, rows_Lower, rows_Upper, rows_IsOut, ...
    'VariableNames', {'Subject','Group','Band','Value','Lower','Upper','IsOutlier'});

outfile = fullfile(stats_dir, 'outliers_iqr_all.csv');
writetable(T, outfile);
fprintf('\n[INFO] Wrote detailed outlier table → %s\n', outfile);