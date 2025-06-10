% ------------------------------------------------------------------------------------------
% Script: eeg_group_analysis_v1.m
% Version: v1  |  Last Updated: June 10, 2025
%
% Description:
% This script performs group-level EEG analysis of normalized bandpower data.
% Includes statistical testing, effect size computation, and topographic plotting.
%
% Key Features:
%   • QC integration from `qc_summary_all_subjects.csv`
%   • Barplots of group-level bandpower (CN, AD, FTD)
%   • One-way ANOVA with post-hoc Tukey HSD comparisons per band
%   • Significance stars on barplots and summary tables
%   • Cohens d effect size calculations across all group pairs
%   • Final summary table: [Band, ANOVA p, Post-hoc Significant, Max Cohen’s d]
%   • Group topographic plots averaged by band and diagnosis group
%
% Output Directory:
%   → All results saved in `outputs/stats/`
%
% Notes:
%   - Requires individual subject bandpower files named: `sub-XXX_v2.1_bandpower_matrix.csv`
%   - Topomap requires `*_v2.1_cleaned_final.set` files with consistent channel order
% ------------------------------------------------------------------------------------------

close all;

% Create stats subfolder inside 'outputs' 
stats_dir = fullfile('outputs', 'stats');
if ~exist(stats_dir, 'dir')
    mkdir(stats_dir);
end

% Load QC + Group Info
qc_file = 'outputs/qc_summary_all_subjects.csv';
group_file = 'participants.tsv';

qc_table = readtable(qc_file);
group_table = readtable(group_file, 'FileType', 'text', 'Delimiter', '\t');

% Merge by subject ID
merged = innerjoin(qc_table, group_table, 'LeftKeys', 'Subject', 'RightKeys', 'participant_id');

% Basic Sanity Check
disp('Subjects per Group:');
disp(groupcounts(merged.Group));

% Print subject counts to console
n_CN  = sum(strcmp(merged.Group, 'C'));
n_AD  = sum(strcmp(merged.Group, 'A'));
n_FTD = sum(strcmp(merged.Group, 'F'));

fprintf('Control (C): %d subjects\n', n_CN);
fprintf('AD (A):      %d subjects\n', n_AD);
fprintf('FTD (F):     %d subjects\n', n_FTD);

% Plot Normalized Bandpower by Group
bands = {'ThetaMean', 'AlphaMean', 'BetaMean'};
group_labels = unique(merged.Group);

for i = 1:length(bands)
    figure;
    boxchart(categorical(merged.Group), merged.(bands{i}));
    title(['Normalized ' bands{i} ' by Group']);
    ylabel('Normalized Power (0–1)');
    xlabel('Group');
end


% Create logical group indices
isCN = strcmp(merged.Group, 'C');
isAD = strcmp(merged.Group, 'A');
isFTD = strcmp(merged.Group, 'F');

% Compute Group Averages + Error Bars
group_labels = {'Control', 'AD', 'FTD'};
group_masks = {isCN, isAD, isFTD};
band_labels = {'Theta', 'Alpha', 'Beta'};

group_means = zeros(3,3);
group_stes  = zeros(3,3);  % Standard Error of the Mean

for g = 1:3
    mask = group_masks{g};
    group_means(g,1) = mean(merged.ThetaMean(mask));
    group_means(g,2) = mean(merged.AlphaMean(mask));
    group_means(g,3) = mean(merged.BetaMean(mask));
    
    group_stes(g,1) = std(merged.ThetaMean(mask)) / sqrt(sum(mask));
    group_stes(g,2) = std(merged.AlphaMean(mask)) / sqrt(sum(mask));
    group_stes(g,3) = std(merged.BetaMean(mask)) / sqrt(sum(mask));
end


% Plot Bar Graph with Error Bars
figure('Position', [100 100 600 400]);
bar_handle = bar(group_means);
hold on;
ngroups = size(group_means, 1);
nbars = size(group_means, 2);

% Add error bars
for b = 1:nbars
    x = bar_handle(b).XEndPoints;
    errorbar(x, group_means(:,b), group_stes(:,b), 'k', 'linestyle', 'none', 'LineWidth', 1);
end

x_labels = {sprintf('Control (n=%d)', n_CN), sprintf('AD (n=%d)', n_AD), sprintf('FTD (n=%d)', n_FTD)};
set(gca, 'XTickLabel', x_labels);
legend(band_labels, 'Location', 'northoutside', 'Orientation', 'horizontal');
ylabel('Normalized Bandpower');
title('Group-Level Bandpower Comparison');
ylim([0 1]);


T_stats = array2table(group_means, 'VariableNames', band_labels, 'RowNames', group_labels);
writetable(T_stats, fullfile(stats_dir, 'group_bandpower_means.csv'), 'WriteRowNames', true);

% ==========Step: ANOVA and Post-Hoc Tests for Group Differences=============
anova_results = table('Size', [length(bands), 2], 'VariableTypes', {'double', 'cell'}, 'VariableNames', {'ANOVA_p', 'PostHoc'}, 'RowNames', bands);

for i = 1:length(bands)
    band = bands{i};
    
    % Extract data and group labels
    y = merged.(band);
    g = merged.Group;

    % Run one-way ANOVA
    [p, tbl, stats] = anova1(y, g, 'off');
    anova_results.ANOVA_p(i) = p;

    % Run post-hoc comparisons 
    posthoc = multcompare(stats, 'Display', 'off');

    % Convert numeric group indices to group names
    group_names = stats.gnames;  % e.g., {'C', 'A', 'F'}
    G1 = group_names(posthoc(:,1));
    G2 = group_names(posthoc(:,2));

    anova_results.PostHoc{i} = table(G1, G2, posthoc(:,3), posthoc(:,4), posthoc(:,5), posthoc(:,6), 'VariableNames', {'Group1', 'Group2', 'LowerCI', 'MeanDiff', 'UpperCI', 'pValue'});
end

% Save ANOVA p-values
writetable(anova_results(:, 'ANOVA_p'), fullfile(stats_dir, 'anova_bandpower_pvalues.csv'), 'WriteRowNames', true);

% Save post-hoc table for all bands
all_posthoc = table();  % empty table

for i = 1:length(bands)
    band = bands{i};
    T_posthoc = anova_results.PostHoc{i};
    T_posthoc.Band = repmat({band}, height(T_posthoc), 1);
    all_posthoc = [all_posthoc; T_posthoc];
end

writetable(all_posthoc, fullfile(stats_dir, 'posthoc_all_bands.csv'));

disp('--- ANOVA Results ---');
disp(anova_results(:, 'ANOVA_p'));

% Add Stars for Post-hoc Significance
% Optional: Define Y-axis buffer
y_offset = 0.05;  % vertical space above bars for the stars

% Coordinates for significance lines
group_idx = [1 2 3];  % CN = 1, AD = 2, FTD = 3

for b = 1:length(bands)
    T_posthoc = anova_results.PostHoc{b};  % one table per band
    band_name = bands{b};

    % Get Y max for this band to place stars above that group’s bar
    max_y = max(group_means(:, b) + group_stes(:, b));

    for r = 1:height(T_posthoc)
        g1 = find(strcmp(T_posthoc.Group1(r), {'C','A','F'}));  % index from CN/AD/FTD
        g2 = find(strcmp(T_posthoc.Group2(r), {'C','A','F'}));

        % Only mark if significant
        pval = T_posthoc.pValue(r);
        if pval < 0.05
            % Determine number of stars
            if pval < 0.001
                stars = '***';
            elseif pval < 0.01
                stars = '**';
            else
                stars = '*';
            end

            % Horizontal line position
            x1 = group_idx(g1); x2 = group_idx(g2);
            y = max_y + y_offset;

            % Draw line
            plot([x1 x2], [y y], 'k-', 'LineWidth', 1);
            % Draw stars
            text(mean([x1 x2]), y + 0.01, stars, 'HorizontalAlignment', 'center', 'FontSize', 14);
            % Bump up for next annotation
            max_y = y + 0.05;
        end
    end
end

% Export to file
saveas(gcf, fullfile(stats_dir, 'group_bandpower_comparison.png'));
exportgraphics(gcf, 'outputs/stats/group_bandpower_comparison.pdf', 'ContentType', 'vector');

% Cohen’s d (All Group Comparisons)
comparisons = {'AD_vs_CN', 'FTD_vs_CN', 'FTD_vs_AD'};
all_d = zeros(length(comparisons), length(bands));

for i = 1:length(bands)
    band = bands{i};

    % AD vs CN
    x1 = merged{isAD, band}; x2 = merged{isCN, band};
    s1 = std(x1); s2 = std(x2); n1 = length(x1); n2 = length(x2);
    pooled_sd = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1 + n2 - 2));
    all_d(1,i) = (mean(x1) - mean(x2)) / pooled_sd;

    % FTD vs CN
    x1 = merged{isFTD, band}; x2 = merged{isCN, band};
    s1 = std(x1); s2 = std(x2); n1 = length(x1); n2 = length(x2);
    pooled_sd = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1 + n2 - 2));
    all_d(2,i) = (mean(x1) - mean(x2)) / pooled_sd;

    % FTD vs AD
    x1 = merged{isFTD, band}; x2 = merged{isAD, band};
    s1 = std(x1); s2 = std(x2); n1 = length(x1); n2 = length(x2);
    pooled_sd = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1 + n2 - 2));
    all_d(3,i) = (mean(x1) - mean(x2)) / pooled_sd;
end

T_d_all = array2table(all_d, 'VariableNames', bands, 'RowNames', comparisons);
writetable(T_d_all, fullfile(stats_dir, 'effect_size_all_pairs.csv'), 'WriteRowNames', true);

% Save Merged Data
writetable(merged, 'outputs/group_merged_summary.csv');

% group_bandpower_summary
% Define readable group labels
group_map = {'C', 'A', 'F'};
group_names = {'CN', 'AD', 'FTD'};

summary_rows = {};

for b = 1:length(bands)
    band = bands{b};
    raw_p = anova_results.ANOVA_p(b);
    if raw_p < 0.001
        anova_p = sprintf('%.3g***', raw_p);
    elseif raw_p < 0.01
        anova_p = sprintf('%.3g**', raw_p);
    elseif raw_p < 0.05
        anova_p = sprintf('%.3g*', raw_p);
    else
        anova_p = sprintf('%.3g', raw_p);
    end

    % Post-hoc (get readable significant comparisons)
    T_posthoc = anova_results.PostHoc{b};
    sig_rows = T_posthoc.pValue < 0.05;

    if any(sig_rows)
        sig_comparisons = arrayfun(@(r) [group_names{T_posthoc.Group1(r)} ' vs ' group_names{T_posthoc.Group2(r)}], find(sig_rows), 'UniformOutput', false);
        posthoc_sig = strjoin(sig_comparisons, ', ');
    else
        posthoc_sig = '–';
    end

    % Max Cohens d (from all_d table)
    max_d = max(abs(all_d(:, b)));

    summary_rows = [summary_rows; {band, anova_p, posthoc_sig, max_d}];
end

T_summary = cell2table(summary_rows, 'VariableNames', {'Band', 'ANOVA_p', 'PostHoc_Significant', 'Max_Cohens_d'});

% Save to file
writetable(T_summary, fullfile(stats_dir, 'bandpower_stats_summary.csv'));

% ------------------- Group Topographic Summary (Combined by Group) -------------------
subjects = merged.Subject;
groups   = merged.Group;  % 'C', 'A', 'F'

band_cols = {'ThetaPower', 'AlphaPower', 'BetaPower'};
band_titles = {'Theta (4–7 Hz)', 'Alpha (8–13 Hz)', 'Beta (13–30 Hz)'};
group_codes = {'C', 'A', 'F'};
group_names = {'Control', 'AD', 'FTD'};

% Load reference chanlocs
ref_sub = subjects{1};
setfile = fullfile('outputs', ref_sub, [ref_sub '_v2.1_cleaned_final.set']);
EEG = pop_loadset(setfile);
chanlocs = EEG.chanlocs;
expected_n = length(chanlocs);

for g = 1:length(group_codes)
    group_code = group_codes{g};
    group_subjects = subjects(strcmp(groups, group_code));
    
    figure('Position', [100 100 1200 400]);
    for b = 1:length(band_cols)
        band = band_cols{b};
        all_vals = [];

        for s = 1:length(group_subjects)
            subj_id = group_subjects{s};
            band_file = fullfile('outputs', subj_id, [subj_id '_v2.1_bandpower_matrix.csv']);
            if isfile(band_file)
                T = readtable(band_file);
                vals = T.(band);
                if length(vals) == expected_n
                    all_vals = [all_vals; vals'];
                end
            end
        end

        subplot(1, 3, b);
        if isempty(all_vals)
            text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center');
        else
            mean_vals = mean(all_vals, 1);
            topoplot(mean_vals, chanlocs, 'maplimits', 'maxmin', 'electrodes', 'on');
            title([group_names{g} ' - ' band_titles{b}]);
            colorbar;
        end
    end

    saveas(gcf, fullfile(stats_dir, ['topomap_group_' lower(group_code) '.png']));
end
