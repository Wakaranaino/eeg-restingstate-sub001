close all; clc;

%% ============================================================
%  0. Setup / load tables
%% ============================================================
stats_dir = fullfile('outputs', 'stats');
if ~exist(stats_dir, 'dir')
    mkdir(stats_dir);
end

qc_file    = 'outputs/qc_summary_all_subjects.csv';
group_file = 'participants.tsv';

qc_table    = readtable(qc_file);
group_table = readtable(group_file, 'FileType','text','Delimiter','\t');

% Merge QC with grouping info
merged = innerjoin(qc_table, group_table, ...
    'LeftKeys','Subject','RightKeys','participant_id');

disp('Subjects per Group:');
disp(groupcounts(merged.Group));

isCN  = strcmp(merged.Group, 'C');
isAD  = strcmp(merged.Group, 'A');
isFTD = strcmp(merged.Group, 'F');

n_CN  = sum(isCN);
n_AD  = sum(isAD);
n_FTD = sum(isFTD);

fprintf('Control (C): %d subjects\n', n_CN);
fprintf('AD (A):      %d subjects\n', n_AD);
fprintf('FTD (F):     %d subjects\n', n_FTD);

% Band columns in merged (already normalized 0–1):
bands        = {'ThetaMean','AlphaMean','BetaMean'};
band_labels  = {'Theta','Alpha','Beta'};     
group_codes  = {'C','A','F'};                


%% ============================================================
%  1. CLEAN per-band distribution plots, 1 figure per band
%% ============================================================
for bi = 1:numel(bands)
    this_band = bands{bi};
    this_label = band_labels{bi};

    fig_box = figure('Position',[100 100 600 400]);

    % boxchart by diagnostic group
    boxchart( categorical(merged.Group), merged.(this_band), ...
        'BoxFaceColor',[0.8 0.85 1.0], ...
        'LineWidth',1.2, ...
        'WhiskerLineColor',[0 0 0], ...
        'MarkerStyle','none');  % turn off default scatter markers

    hold on;

    % add light jittered subject dots so it's not just a box
    % we jitter manually to avoid MATLAB's default legend spam
    cats = categories(categorical(merged.Group)); 
    for ci = 1:numel(cats)
        mask_here = strcmp(merged.Group, cats{ci});
        x_pos = ci + (rand(sum(mask_here),1)-0.5)*0.15; % jitter
        y_pos = merged.(this_band)(mask_here);
        scatter(x_pos, y_pos, 20, 'k', 'filled', ...
            'MarkerFaceAlpha',0.4, 'MarkerEdgeColor','none', ...
            'HandleVisibility','off');
    end

    xlabel('Group');
    ylabel('Normalized Power (0–1)');
    title(['Normalized ' this_label ' by Group']);
    ylim([0 1]); % keep consistent visual scale across bands

    % save, then close so nothing leaks into next plot
    out_png = fullfile(stats_dir, ['boxplot_' this_band '.png']);
    saveas(fig_box, out_png);
    close(fig_box);
end


%% ============================================================
%  2. Compute group means / SEM for summary bar plot + stats
%% ============================================================
maskC = isCN;
maskA = isAD;
maskF = isFTD;

group_means = zeros(3,3); % rows: C,A,F  cols: Theta,Alpha,Beta
group_stes  = zeros(3,3);

for bi = 1:numel(bands)
    Y = merged.(bands{bi});
    group_means(1,bi) = mean(Y(maskC));
    group_means(2,bi) = mean(Y(maskA));
    group_means(3,bi) = mean(Y(maskF));

    group_stes(1,bi)  = std(Y(maskC)) / sqrt(sum(maskC));
    group_stes(2,bi)  = std(Y(maskA)) / sqrt(sum(maskA));
    group_stes(3,bi)  = std(Y(maskF)) / sqrt(sum(maskF));
end


%% ============================================================
%  3. Run ANOVA + posthoc for each band (C,A,F)
%     save ANOVA p, posthoc tables, etc.
%% ============================================================
anova_results = table('Size',[numel(bands),2], ...
    'VariableTypes',{'double','cell'}, ...
    'VariableNames',{'ANOVA_p','PostHoc'}, ...
    'RowNames',bands);

all_posthoc = table(); % collect all bands’ posthocs

for bi = 1:numel(bands)
    Y = merged.(bands{bi});   % data values
    G = merged.Group;         % 'C','A','F'

    [p,~,stats] = anova1(Y, G, 'off');
    anova_results.ANOVA_p(bi) = p;

    ph = multcompare(stats, 'Display','off');
    % ph columns: [i j lower mean upper pval]
    % stats.gnames gives the label for group index i/j

    gnames = stats.gnames;  % e.g. {'A','C','F'} depending on MATLAB's internal order
    G1 = gnames(ph(:,1));
    G2 = gnames(ph(:,2));

    T_posthoc = table(G1, G2, ph(:,3), ph(:,4), ph(:,5), ph(:,6), ...
        'VariableNames',{'Group1','Group2','LowerCI','MeanDiff','UpperCI','pValue'});

    anova_results.PostHoc{bi} = T_posthoc;

    BandCol = repmat(bands(bi), height(T_posthoc), 1);
    this_block = [T_posthoc table(BandCol,'VariableNames',{'Band'})];
    all_posthoc = [all_posthoc; this_block];
end

% Save ANOVA table + all posthoc to disk
writetable(anova_results(:, 'ANOVA_p'), ...
    fullfile(stats_dir,'anova_bandpower_pvalues.csv'), ...
    'WriteRowNames',true);

writetable(all_posthoc, ...
    fullfile(stats_dir,'posthoc_all_bands.csv'));

disp('--- ANOVA p-values by band ---');
disp(anova_results(:, 'ANOVA_p'));


%% ============================================================
%  4. MAIN SUMMARY FIGURE
%     - One grouped bar plot:
%         x-axis: CN, AD, FTD
%         bar colors: Theta / Alpha / Beta
%     - SEM error bars
%     - Significance markers ONLY vs CN:
%           AD vs CN
%           FTD vs CN
%       (NO AD vs FTD)
%     - Legend should only show Theta/Alpha/Beta. No "data1,data2..."
%% ============================================================

fig_sum = figure('Position',[100 100 750 480]);

% First draw the bars
bh = bar(group_means); hold on;

% Give each band a stable color:
band_rgb = [0.2 0.4 0.8;   % Theta = bluish
            0.9 0.4 0.1;   % Alpha = orange
            0.95 0.8 0.2]; % Beta  = yellow-ish
for bi = 1:numel(bands)
    bh(bi).FaceColor = band_rgb(bi,:);
end

% Add SEM error bars
for bi = 1:numel(bands)
    x = bh(bi).XEndPoints;
    errorbar(x, group_means(:,bi), group_stes(:,bi), ...
        'k','linestyle','none','LineWidth',1, ...
        'HandleVisibility','off'); % keep these out of legend
end

% Label axes
set(gca,'XTickLabel', {sprintf('CN (n=%d)', n_CN), ...
                       sprintf('AD (n=%d)', n_AD), ...
                       sprintf('FTD (n=%d)', n_FTD)});
ylabel('Normalized Bandpower (0–1)');
title('Group-Level Bandpower Comparison');
ylim([0 1]);  % temporary, we may push higher after plotting sig

% Lock in legend *now*, BEFORE we add bracket lines.
legend(bh, band_labels, ...
    'Location','northoutside','Orientation','horizontal');

% Now compute significance vs CN and add brackets
% We'll look inside anova_results.PostHoc for:
%   "C" vs "A"
%   "C" vs "F"
% but posthoc might list ('A','C') instead of ('C','A'), etc.
%
% We'll do per band (Theta,Alpha,Beta) so each bracket sits near that bar's peak.
y_offset = 0.03; % spacing above tallest bar of that band

for bi = 1:numel(bands)
    T_posthoc = anova_results.PostHoc{bi};

    % grab p for CN vs AD
    mask_CA = ( (strcmp(T_posthoc.Group1,'C') & strcmp(T_posthoc.Group2,'A')) | ...
                (strcmp(T_posthoc.Group1,'A') & strcmp(T_posthoc.Group2,'C')) );
    % grab p for CN vs FTD
    mask_CF = ( (strcmp(T_posthoc.Group1,'C') & strcmp(T_posthoc.Group2,'F')) | ...
                (strcmp(T_posthoc.Group1,'F') & strcmp(T_posthoc.Group2,'C')) );

    % baseline height for THIS band = max(mean+SEM across CN,AD,FTD)
    y_base = max(group_means(:,bi) + group_stes(:,bi));
    y_top_for_band = y_base; % track where we last placed a bracket

    compare_list = {
        'C','A', mask_CA;   % CN vs AD
        'C','F', mask_CF;   % CN vs FTD
    };

    for ci = 1:size(compare_list,1)
        gA = compare_list{ci,1};  % 'C'
        gB = compare_list{ci,2};  % 'A' or 'F'
        mask_here = compare_list{ci,3};

        if any(mask_here)
            pval_here = T_posthoc.pValue(mask_here);
            pval_here = pval_here(1); % scalar

            if pval_here < 0.05
                % convert p to stars
                if pval_here < 0.001
                    stars = '***';
                elseif pval_here < 0.01
                    stars = '**';
                else
                    stars = '*';
                end

                % x positions: C -> 1, A -> 2, F -> 3
                xA = find(strcmp(gA, group_codes)); % 1
                xB = find(strcmp(gB, group_codes)); % 2 or 3

                % y position for bracket line
                y_line = y_top_for_band + y_offset;

                % draw bracket line (keep it OUT of legend)
                plot([xA xB],[y_line y_line], 'k-', 'LineWidth',1.2, ...
                    'HandleVisibility','off');

                % star text
                text(mean([xA xB]), y_line + 0.01, stars, ...
                    'HorizontalAlignment','center','FontSize',12, ...
                    'HandleVisibility','off');

                % update top so next bracket for this band goes higher
                y_top_for_band = y_line + 0.05;
            end
        end
    end

    % after adding both possible brackets for this band,
    % we might have extended beyond current ylim. update ylim if needed:
    current_ylim = ylim;
    if y_top_for_band + 0.08 > current_ylim(2)
        ylim([current_ylim(1), y_top_for_band + 0.08]);
    end
end

% save combined summary figure
saveas(fig_sum, fullfile(stats_dir,'group_bandpower_comparison.png'));
exportgraphics(fig_sum, fullfile(stats_dir,'group_bandpower_comparison.pdf'), ...
    'ContentType','vector');
close(fig_sum);


%% ============================================================
%  5. Effect sizes (Cohen's d) for reporting
%     We'll still compute:
%        AD vs CN
%        FTD vs CN
%        FTD vs AD
%     and write to CSV
%% ============================================================
comparisons = {'AD_vs_CN','FTD_vs_CN','FTD_vs_AD'};
all_d = zeros(numel(comparisons), numel(bands));

for bi = 1:numel(bands)
    band = bands{bi};

    % AD vs CN
    x1 = merged{isAD, band}; x2 = merged{isCN, band};
    s1 = std(x1); s2 = std(x2); n1 = length(x1); n2 = length(x2);
    pooled_sd = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2));
    all_d(1,bi) = (mean(x1)-mean(x2))/pooled_sd;

    % FTD vs CN
    x1 = merged{isFTD, band}; x2 = merged{isCN, band};
    s1 = std(x1); s2 = std(x2); n1 = length(x1); n2 = length(x2);
    pooled_sd = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2));
    all_d(2,bi) = (mean(x1)-mean(x2))/pooled_sd;

    % FTD vs AD
    x1 = merged{isFTD, band}; x2 = merged{isAD, band};
    s1 = std(x1); s2 = std(x2); n1 = length(x1); n2 = length(x2);
    pooled_sd = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2));
    all_d(3,bi) = (mean(x1)-mean(x2))/pooled_sd;
end

T_d_all = array2table(all_d, ...
    'VariableNames', bands, ...
    'RowNames', comparisons);
writetable(T_d_all, fullfile(stats_dir,'effect_size_all_pairs.csv'), ...
    'WriteRowNames', true);


%% ============================================================
%  6. Human-readable summary table
%     - ANOVA p with stars
%     - Which pairs are sig (<0.05) for each band
%     - Max |d|
%% ============================================================
summary_rows = {};

for bi = 1:numel(bands)
    raw_p = anova_results.ANOVA_p(bi);
    if     raw_p < 0.001, p_txt = sprintf('%.3g***', raw_p);
    elseif raw_p < 0.01,  p_txt = sprintf('%.3g**',  raw_p);
    elseif raw_p < 0.05,  p_txt = sprintf('%.3g*',   raw_p);
    else,                 p_txt = sprintf('%.3g',    raw_p);
    end

    T_posthoc = anova_results.PostHoc{bi};
    sig_mask  = T_posthoc.pValue < 0.05;

    if any(sig_mask)
        sig_pairs = arrayfun(@(r) ...
            sprintf('%s vs %s', T_posthoc.Group1{r}, T_posthoc.Group2{r}), ...
            find(sig_mask), 'UniformOutput', false);
        posthoc_sig = strjoin(sig_pairs, ', ');
    else
        posthoc_sig = '–';
    end

    max_d = max(abs(all_d(:,bi)));

    summary_rows = [summary_rows;
        {band_labels{bi}, p_txt, posthoc_sig, max_d}];
end

T_summary = cell2table(summary_rows, ...
    'VariableNames', {'Band','ANOVA_p','PostHoc_Significant','Max_Cohens_d'});

writetable(T_summary, fullfile(stats_dir,'bandpower_stats_summary.csv'));


%% ============================================================
%  7. Group topographic maps (Theta/Alpha/Beta scalp averages)
%% ============================================================
subjects = merged.Subject;
groups   = merged.Group;  % 'C','A','F'

band_cols          = {'ThetaPower','AlphaPower','BetaPower'};
band_titles        = {'Theta (4–7 Hz)','Alpha (8–13 Hz)','Beta (13–30 Hz)'};
group_codes_plot   = {'C','A','F'};
group_names_pretty = {'Control','AD','FTD'};

% reference chanlocs from first subject (assumes same montage)
ref_sub = subjects{1};
setfile = fullfile('outputs', ref_sub, [ref_sub '_v2.1_cleaned_final.set']);
EEG = pop_loadset(setfile);
chanlocs = EEG.chanlocs;
expected_n = length(chanlocs);

for gi = 1:numel(group_codes_plot)
    gcode = group_codes_plot{gi};
    these_subs = subjects(strcmp(groups, gcode));

    fig_topo = figure('Position',[100 100 1200 400]);

    for bi = 1:numel(band_cols)
        this_col = band_cols{bi};
        all_vals = [];

        for si = 1:numel(these_subs)
            sid = these_subs{si};
            band_file = fullfile('outputs', sid, [sid '_v2.1_bandpower_matrix.csv']);
            if isfile(band_file)
                Tsub = readtable(band_file);
                if ismember(this_col, Tsub.Properties.VariableNames)
                    vals = Tsub.(this_col);
                    if numel(vals) == expected_n
                        all_vals = [all_vals; vals']; %#ok<AGROW>
                    end
                end
            end
        end

        subplot(1, numel(band_cols), bi);
        if isempty(all_vals)
            text(0.5,0.5,'No Data','HorizontalAlignment','center');
        else
            mean_vals = mean(all_vals,1);
            topoplot(mean_vals, chanlocs, ...
                'maplimits','maxmin', ...
                'electrodes','on');
            title([group_names_pretty{gi} ' - ' band_titles{bi}]);
            colorbar;
        end
    end

    saveas(fig_topo, fullfile(stats_dir, ['topomap_group_' lower(gcode) '.png']));
    close(fig_topo);
end

disp('[INFO] Finished: cleaned single-band plots (no stars), main summary bar fig (only CN comparisons, clean legend), stats tables, and topomaps.');
