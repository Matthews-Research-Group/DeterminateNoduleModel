%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains supporting information of the following publication:
%
% Kinetic model of a determinate nodule reveals plant metabolic characteristics 
% for more efficient nitrogen fixation symbiosis
%
% by Rourou Ji, Joshua A.M. Kaste and Megan L. Matthews
%
% This MATLAB script allows to reproduce supplemental figures
%
% Author: rourouj2@illinois.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Some analysis
clear all

addpath('../kineticmodel/');  
[metabolite_names, reaction_names, x0, tspan, options] = Parameter_settings();
Vmax = Vmaxvalues(); 

param_names = fieldnames(Vmax); 

%% load data
load('../data/alldata.mat');

%% Figure S4
all = Vmax_matrix_initial';
filtered = (newdata.Vmax)';

num_variables = size(param_names,1);

figure('Units', 'pixels', 'Position', [200, 200, 850, 900]);

t = tiledlayout(ceil(num_variables / 8), 8, 'Padding', 'compact', 'TileSpacing', 'compact');  % 4 列布局

for i = 1:num_variables
    ax = nexttile;
    
    all_var = all(:,i);
    filtered_var = filtered(:, i);
    
    % Determine global range for binning
    combined = [all_var; filtered_var];
    min_val = min(combined);
    max_val = max(combined);
    bin_edges = logspace(log10(min_val), log10(1000), 30);

    % histogram
    histogram(ax, all_var, 'Normalization', 'count', 'BinEdges', bin_edges, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.6, 'DisplayName', 'Original Data');
    hold(ax, 'on');
    histogram(ax, filtered_var, 'Normalization', 'count', 'BinEdges', bin_edges, 'FaceColor', [0.9 0.2 0.2], 'FaceAlpha', 0.6, 'DisplayName', 'Filtered Data');

    set(ax, 'XScale', 'log');
    xlim(ax, [0.1 1000]); 
    xticks(ax, [0.1 1 10 100 1000]);
    title(ax, param_names{i});
    hold(ax,'off');
end

%% Figure S5
without_set = load('../data/wo_proteomicdata_result.mat');

with_beforefilter = load('../data/alldata.mat');
with_afterfilter = with_beforefilter.newdata;

%%
with_before_gcgn = with_beforefilter.gCgN;
with_after_gcgn = with_afterfilter.gCgN;

with_before_Nfix = with_beforefilter.Nfix;
with_after_Nfix = with_afterfilter.N_ExportRate;

%% plot histogram together

figure('Units', 'pixels', 'Position', [200, 200, 850, 800]);
t = tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'none');  % Removes space between and around plots

% Figure 1
% subplot(2,2,1);
ax1=nexttile;
hold(ax1, 'on'); 
h1=histogram(ax1, without_set.without_before_gcgn, 'BinWidth', 1, 'FaceColor', [222, 137, 255]/256,...
    'FaceAlpha', 0.4, 'Normalization', 'count'); 
h2=histogram(ax1, with_before_gcgn, 'BinWidth', 1, 'FaceColor', [0.85, 0.6, 0.25],...
    'FaceAlpha', 0.4, 'Normalization', 'count');

legend(ax1, {'w/o proteomic data', 'w/ proteomic data'}, 'FontSize', 12, 'Location', 'northeast');

% Add boundary lines
xline(ax1, 2, '--k', 'LB=2', 'LabelHorizontalAlignment', 'right', 'FontSize', 15, 'HandleVisibility', 'off');
xline(ax1, 8, '--k', 'UB=8', 'LabelHorizontalAlignment', 'right', 'FontSize', 15, 'HandleVisibility', 'off');
xlim(ax1, [1,200]);

% >>> remove the old ylim or set it dynamically <<<
ymax1 = max(h2.Values);
if ~isempty(ymax1) && isfinite(ymax1)
    ylim(ax1, [0, ymax1*1.15]);  % add a little headroom
end

xlabel(ax1, '$\mathbf{Nitrogen \; fixation \; efficiency \; (gC g}^{-1}\mathbf{N)}$', ...
       'FontSize', 16, ...
       'Interpreter', 'latex');
ylabel(ax1, '\textbf{Counts}',...
    'FontWeight','bold',...
    'Interpreter', 'latex');
ax1.FontSize = 16;
ax1.XScale = 'log';

% Figure 2
% subplot(2,2,2);
ax2=nexttile;
hold(ax2, 'on'); 

h3=histogram(ax2, without_set.without_before_Nfix, 'BinWidth', 0.1, 'FaceColor', [222, 137, 255]/256,...
    'FaceAlpha', 0.4, 'Normalization', 'count'); 
h4=histogram(ax2, with_before_Nfix, 'BinWidth', 0.1, 'FaceColor', [0.2549, 0.5961, 0.6745],...
    'FaceAlpha', 0.4, 'Normalization', 'count'); 

legend(ax2, {'w/o proteomic data', 'w/ proteomic data'}, 'FontSize', 12, 'Location', 'northeast');

xlim(ax2, [1e-1,12]);

% >>> remove the old ylim or set it dynamically <<<
ymax2 = max(h4.Values);
if ~isempty(ymax2) && isfinite(ymax2)
    ylim(ax2, [0, ymax2*1.1]);  % add a little headroom
end

xlabel(ax2, '$\textbf{Nitrogen fixation rate (mM s}^{-1}\textbf{)}$', ...
       'FontSize', 16, ...
       'Interpreter', 'latex');
ylabel(ax2, '\textbf{Counts}',...
    'FontWeight','bold',...
    'Interpreter', 'latex');
ax2.FontSize = 16;
ax2.XScale = 'log';


% Figure 3
% subplot(2,2,3);
ax3=nexttile;
hold(ax3, 'on'); 

h5=histogram(ax3, without_set.without_after_gcgn, 'BinWidth', 1, 'FaceColor', [222, 137, 255]/256,...
    'FaceAlpha', 0.4, 'Normalization', 'count');
h6=histogram(ax3, with_after_gcgn, 'BinWidth', 1, 'FaceColor', [0.85, 0.6, 0.25],...
    'FaceAlpha', 0.4, 'Normalization', 'count');

legend(ax3, {'w/o proteomic data', 'w/ proteomic data'}, 'FontSize', 12, 'Location', 'northeast');

% Add boundary lines
xline(ax3, 2, '--k', 'LB=2', 'LabelHorizontalAlignment', 'right', 'FontSize', 15, 'HandleVisibility', 'off');
xline(ax3, 8, '--k', 'UB=8', 'LabelHorizontalAlignment', 'right', 'FontSize', 15, 'HandleVisibility', 'off');
xlim(ax3, [1,200]);

% >>> remove the old ylim or set it dynamically <<<
ymax3 = max(h6.Values);
if ~isempty(ymax3) && isfinite(ymax3)
    ylim(ax3, [0, ymax3*1.15]);  % add a little headroom
end
xlabel(ax3, '$\mathbf{Nitrogen \; fixation \; efficiency \; (gC g}^{-1}\mathbf{N)}$', ...
       'FontSize', 16, ...
       'Interpreter', 'latex');
ylabel(ax3, '\textbf{Counts}',...
    'FontWeight','bold',...
    'Interpreter', 'latex');
ax3.FontSize = 16;
ax3.XScale = 'log';


% Figure 4
% subplot(2,2,4);
ax4=nexttile;
hold(ax4, 'on'); 

h7=histogram(ax4, without_set.without_after_Nfix, 'BinWidth', 0.1, 'FaceColor', [222, 137, 255]/256,...
    'FaceAlpha', 0.4, 'Normalization', 'count');
h8=histogram(ax4, with_after_Nfix, 'BinWidth', 0.1, 'FaceColor', [0.2549, 0.5961, 0.6745],...
    'FaceAlpha', 0.4, 'Normalization', 'count'); 

legend(ax4, {'w/o proteomic data', 'w/ proteomic data'}, 'FontSize', 12, 'Location', 'northeast');

xlim(ax4, [1e-1,12]);

% >>> remove the old ylim or set it dynamically <<<
ymax4 = max(h4.Values);
if ~isempty(ymax4) && isfinite(ymax4)
    ylim(ax4, [0, ymax2*1.1]);  % add a little headroom
end

xlabel(ax4, '$\textbf{Nitrogen fixation rate (mM s}^{-1}\textbf{)}$', ...
       'FontSize', 16, ...
       'Interpreter', 'latex');
ylabel(ax4, '\textbf{Counts}',...
    'FontWeight','bold',...
    'Interpreter', 'latex');
ax4.FontSize = 16;
ax4.XScale = 'log';

labels = {'A.', 'B.', 'C.', 'D.'};
axes_handles = [ax1, ax2, ax3, ax4];
for i = 1:length(axes_handles)
    ax = axes_handles(i);
    xlims = xlim(ax);
    ylims = ylim(ax);
    
    % Position text near top-left, with some padding
    text(ax, xlims(1) + 0.001*range(xlims), ...
             ylims(2) - 0.015*range(ylims), ...
             labels{i}, ...
             'FontWeight', 'bold', ...
             'FontSize', 16, ...
             'HorizontalAlignment', 'left', ...
             'VerticalAlignment', 'top');
end


%% Figure S6
x1 = categorical(metabolite_names, metabolite_names, 'Ordinal', true);
y1 = newdata.x;

% Reshape into column vector
y_all = y1(:);  

% Repeat categorical labels for each sample
num_samples = size(y1, 2); 
x_all = repmat(x1', num_samples, 1);

% Define the index ranges for each run
index_ranges = {
    1:11, ...%glycolysis
    12:18, ...%PP
    19:29, ...%TCA
    [31:33, 41, 49],...%NH4, GS/GOGAT, ALN and ASP
    [35:37,42:45], ...%glycine and serine metabolism
    [34,38:40,46:48, 50:60]%Purine pathway
};

% Define custom titles for each subplot
title_names = {
    'Glycolysis', ...
    'PP pathway', ...
    'TCA', ...
    'GSGOGAT, alanine and aspartate metabolism', ...
    'Glycine and serine metabolism',...
    'De novo purine metabolism'
};

figure('Position', [1, 700, 700, 1000]);
num_plots = length(index_ranges); 
t = tiledlayout(ceil(num_plots / 2), 2, 'Padding', 'compact', 'TileSpacing', 'compact'); 

% Loop over each index range and generate the plots
for i = 1:length(index_ranges)
    indices = index_ranges{i};
    selected_reactions = metabolite_names(indices); 
    
    selected_categories = ismember(x_all, selected_reactions);
    
    x_selected = x_all(selected_categories);
    y_selected = y_all(selected_categories);
    
    x_selected = categorical(x_selected, selected_reactions, 'Ordinal', true);

    ax = nexttile; 
    hold(ax, 'on');

    %%%%%%%%%%% boxchart %%%%%%%%%%%%%%%%%%%%%%
    b = boxchart(ax, x_selected, y_selected);
    b.MarkerStyle = '+';   % hide outlier markers
    b.MarkerSize = 1.3;
    b.BoxFaceAlpha = 0.3;      % a little transparency (optional)
    b.LineWidth = 1.1;
    b.WhiskerLineStyle = '-';
    
    set(gca, 'YScale', 'log');   % make y-axis logarithmic
    
    ylabel('Metabolic concentration (mM)');
    title(ax, title_names{i});

    xticklabels(ax, selected_reactions);  
    set(ax, 'XTickLabelRotation', 90);
    hold(ax, 'off');
end


%% Figure S7
x1 = categorical(reaction_names, reaction_names, 'Ordinal', true);
y1 = newdata.v;

% Reshape into column vector
y_all = y1(:); 

num_samples = size(y1, 2); 
x_all = repmat(x1', num_samples, 1);

index_ranges = {
    1:12, ...%glycolysis
    13:20, ...%PP
    [22:23,25:27,30:32,34,36,28:29], ...%TCA and glyoxylate cycle
    [24,33,35, 61], ...%transports and sink
    [37:40], ... %GSGOGAT
    [42:44, 50:52], ...%glycine
    [48:49, 56], ... %ALN and ASP
    [21, 41, 45:47, 53:55, 57:60, 62:68]%purine
};

title_names = {
    'Glycolysis', ...
    'PP pathway', ...
    'TCA and glyoxylate cycle', ...
    'Transports and sink', ...
    'GSGOGAT', ...
    'Glycine and serine metabolism',...
    'Alanine and aspartate metabolism',...
    'De novo purine metabolism'
};

figure('Position', [1, 700, 700, 1000]);
num_plots = length(index_ranges);  
t = tiledlayout(ceil(num_plots / 2), 2, 'Padding', 'compact', 'TileSpacing', 'compact'); 

for i = 1:length(index_ranges)
    indices = index_ranges{i};
    selected_reactions = reaction_names(indices); 
    
    selected_categories = ismember(x_all, selected_reactions);
    
    x_selected = x_all(selected_categories);
    y_selected = y_all(selected_categories);
    
    x_selected = categorical(x_selected, selected_reactions, 'Ordinal', true);
    
    ax = nexttile; 
    hold(ax, 'on');

    %%%%%%%%%%% boxchart %%%%%%%%%%%%%%%%%%%%%%
    b = boxchart(ax, x_selected, y_selected);
    b.MarkerStyle = '+';   % hide outlier markers
    b.MarkerSize = 1.3;
    b.BoxFaceAlpha = 0.3;      % a little transparency (optional)
    b.LineWidth = 1.1;
    b.WhiskerLineStyle = '-';
    b.BoxFaceColor = [0 0.2 0.4];%[0.4940 0.1840 0.5560];
    b.MarkerColor = [0 0.2 0.4];%[0.4940 0.1840 0.5560]; 
    
    ylabel('Flux (mM/s)');
    title(ax, title_names{i});
    set(ax, 'XTickLabelRotation', 90);  
    hold(ax, 'off');

end

%% Figure S8
load('../data/fitness_comparison.mat');

log_exp = log10(fitness_comparison.E_prop_toPDH_original);
log_model_median = log10(fitness_comparison.E_median_prop_toPDH);
p = polyfit(log_exp, log_model_median, 1);  
y_fit = polyval(p, log_exp);
SS_res = sum((log_model_median - y_fit).^2);
SS_tot = sum((log_model_median - mean(log_model_median)).^2);
R_squared = 1 - SS_res / SS_tot;

figure;
scatter(log_exp, log_model_median, 'k','filled');
hold on;
plot(log_exp, y_fit, 'k-', 'LineWidth', 1.5);
xlabel('Experimental ratio', 'FontSize', 14);
yl=ylabel('Simulated ratio','Rotation', 90, 'FontSize', 14);
yl.Position(1) = yl.Position(1) + 1.9;  
yl.Position(2) = yl.Position(2) + 1.7;  
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
annotation('textbox', [0.7, 0.75, 0.1, 0.05], ...
           'String', sprintf('slope = %.4f\nR^2 = %.4f', p(1), R_squared), ...
           'FitBoxToText', 'on', ...
           'FontSize', 14, ...  
           'BackgroundColor', 'w');

xlim([-2,1.5]);
ylim([-2,1.5]);
axis equal; grid on;


%% Figure S15
[efficient_nod, inefficient_nod] = efficient_identify2(newdata);
energy_ef = datastruct(efficient_nod);
energy_inef = datastruct(inefficient_nod);

data = {energy_ef.NetATP_mM_s, energy_inef.NetATP_mM_s};  % each x# is a 1D vector
labels = {'Efficient nodules', 'Inefficient nodules'};

figure; hold on;
colors = lines(numel(data));  % Distinct colors
for i = 1:numel(data)
    [f, xi] = ksdensity(data{i});  % Kernel density estimation
    plot(xi, f, 'LineWidth', 2, 'Color', colors(i,:));
end
xlabel('NetATP');
ylabel('Density');

title('Net ATP variation for efficient and inefficient nodules');
legend(labels, 'Location', 'northwest');


%% Figure S16 & S17
load('../data/pw_allsets_analysis.mat');
nCombo = numel(comboKeys);

% distribution of % change for selected pairs 
requested = { ...
    'PYK-',  'ZWF-'; ...
    'ZWF-',  'PEPC+'; ...
    'PFK+',  'PYK-'; ...
    'PFK+',  'PEPC+'; ...
    'PYK-',  'TAL-'; ...
    'TAL-',  'PEPC+'; ...
    'PYK-',  'PGL-'; ...
    'PGL+',  'PEPC+'; ...
    'ZWF-',  'PGL-'; ...
    'PFK+',  'ZWF-'; ...
    'ZWF-',  'MALOut+' ...
};
nReq = size(requested, 1);

% locate each requested pair in pairLabels (try both orderings)
reqIdx    = NaN(nReq, 1);
reqLabels = cell(nReq, 1);
for i = 1:nReq
    lab1 = sprintf('%s %s', requested{i,1}, requested{i,2});
    lab2 = sprintf('%s %s', requested{i,2}, requested{i,1});
    hit  = find(strcmp(pairLabels, lab1) | strcmp(pairLabels, lab2), 1);
    if ~isempty(hit)
        reqIdx(i)    = hit;
        reqLabels{i} = pairLabels{hit};
    else
        warning('Pair not found in pairLabels: "%s" / "%s"', lab1, lab2);
    end
end
ok        = ~isnan(reqIdx);
reqIdx    = reqIdx(ok);
reqLabels = reqLabels(ok);
nFound    = numel(reqIdx);

nCols = min(nFound, 3);
nRows = ceil(nFound / nCols);

colG = [212, 127, 102] / 255;
colV = [77, 116, 140] / 255;

% same bin width across every subplot so bars are comparable
binW_gcgn  = 3.5;   % % units for gCgN change
binW_vNout = 3.5;   % % units for vNout change

% font sizes for the distribution subplots (raise these to taste)
axFS    = 18;   % axis tick numbers
labFS   = 24;   % x / y axis labels
titleFS = 25;   % per-tile title

% ---- Figure: gCgN % change distributions ----
figure('Position', [100 100 1000 900], 'Name', 'gCgN % change distribution');
tl1 = tiledlayout(nRows, nCols, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:nFound
    c  = reqIdx(i);
    T  = pairTables.(comboKeys{c});
    dg = T.gcgn_change_pct(~isnan(T.gcgn_change_pct));

    nexttile;
    if ~isempty(dg)
        histogram(dg, 'Normalization', 'pdf', ...
                  'BinWidth', binW_gcgn, ...
                  'FaceColor', colG, 'EdgeColor', 'w');
    end
    xline(0, 'k--', 'LineWidth', 1);
    set(gca, 'FontSize', axFS);
    xlabel('Change in Nitrogen fixation efficiency (%) ', 'FontSize', labFS);
    ylabel('Density', 'FontSize', labFS);
    % ylim([0, 1]);

    title(reqLabels{i}, 'FontSize', titleFS);
    box on; grid on;
end
% title(tl1, 'Change in Nitrogen fixation efficiency % ', ...
%       'FontWeight', 'bold', 'FontSize', 14);

% ---- Figure: vNout % change distributions ----
figure('Position', [150 150 1000 900], 'Name', 'vNout % change distribution');
tl2 = tiledlayout(nRows, nCols, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:nFound
    c  = reqIdx(i);
    T  = pairTables.(comboKeys{c});
    dv = T.vNout_change_pct(~isnan(T.vNout_change_pct));

    nexttile;
    if ~isempty(dv)
        histogram(dv, 'Normalization', 'pdf', ...
                  'BinWidth', binW_vNout, ...
                  'FaceColor', colV, 'EdgeColor', 'w');
    end
    xline(0, 'k--', 'LineWidth', 1);
    set(gca, 'FontSize', axFS);
    xlabel('Change in Nitrogen fixation rate (%) ', 'FontSize', labFS);
    ylabel('Density', 'FontSize', labFS);
    % ylim([0, 1]);
    title(reqLabels{i}, 'FontSize', titleFS);
    box on; grid on;
end
% title(tl2, 'Nitrogen fixation rate % change', ...
%       'FontWeight', 'bold', 'FontSize', 14);


%% Figure S18
load('../data/results_allsets_singleE.mat');
nCombo    = numel(comboKeys);

min_change_pct = 0.1;
nShow = Inf;    % of what survives the filter, keep the nShow biggest |gcgn median|

v1 = gcgn_median_change;    % negative = gcgn drops below baseline -> bar extends left
v2 = vNout_median_change;   % positive = N export above baseline   -> bar extends right

isBig = (abs(v1) >= min_change_pct) | (abs(v2) >= min_change_pct);
isBig = isBig & ~isnan(v1) & ~isnan(v2);

% enzymes to leave out of the figure (exact enzyme name, direction-agnostic:
% both the 'up' and 'down' combo of each listed enzyme are dropped)
exclude = {'HKI'};
enzName = regexprep(string(comboKeys(:)), '(up|down)$', '');   % 'HKIup' -> 'HKI'
isBig   = isBig & ~ismember(enzName, string(exclude));

fprintf('Plotting %d / %d combos (|median change| >= %.2g%% in gcgn or N export; excluding %s).\n', ...
    nnz(isBig), nCombo, min_change_pct, strjoin(exclude, ', '));
if ~any(isBig)
    error('No combo reaches %.2g%% -- nothing to plot.', min_change_pct);
end

% sort so the largest gcgn reduction sits on top, then cut the negligible ones
v1_forSort = v1;
v1_forSort(~isBig) = +Inf;              % park them at the end
[~, order] = sort(v1_forSort, 'ascend');
order = order(isBig(order));

if isfinite(nShow) && nShow < numel(order)
    % biggest movers in either direction, then re-sorted for the plot
    [~, byMag] = sort(abs(v1), 'descend', 'MissingPlacement','last');
    keep  = byMag(1:nShow);
    order = order(ismember(order, keep));
end

v1 = v1(order);
v2 = v2(order);
enznames_ordered = enzLabels(order);
nn = numel(order);

% asymmetric 5th-95th percentile whiskers (currently not drawn -- bars only).
% Uncomment these and the errorbar() calls below to show the percentile band.
% e1neg = v1 - gcgn_p5(order);    e1pos = gcgn_p95(order)  - v1;
% e2neg = v2 - vNout_p5(order);   e2pos = vNout_p95(order) - v2;

% one row per combo -> scale the height with the number of rows
figure('Position',[200 60 900 max(400, 20*nn)]);

ax = gca;
ax.Color = [240 240 240]/255;   % light gray background
ax.FontSize = 15;
box(ax, 'on');
ax.LineWidth = 3;

hold on;

hb1 = barh(1:nn, v1, 0.6, 'FaceColor',[212, 127, 102] / 255);
hb2 = barh(1:nn, v2, 0.6, 'FaceColor',[77, 116, 140] / 255);

% vertical separators
for i = 0.5 : 1 : nn+0.5
    yline(i, 'w-', 'LineWidth',1.6);
end

% 5th-95th percentile band -- commented out: median bars only
% errorbar(v1, 1:nn, e1neg, e1pos, 'horizontal', 'LineStyle','none', ...
%     'Color',[0.15 0.15 0.15], 'LineWidth',1.4, 'CapSize',16);
% errorbar(v2, 1:nn, e2neg, e2pos, 'horizontal', 'LineStyle','none', ...
%     'Color',[0.15 0.15 0.15], 'LineWidth',1.4, 'CapSize',16);

% axis tick fonts: the x-axis numbers and the y-axis enzyme labels share this
% size. It shrinks as more rows are plotted; raise the two bounds for bigger text.
tickFS = max(27, min(40, round(600/nn)));
set(gca,'YTick',1:nn,'YTickLabel',enznames_ordered,'YDir','reverse','FontSize', tickFS);
xline(0,'k','LineWidth',1.5);

% axis fits the median bars (no whiskers drawn), with margin for the value labels
xhi = max([v1; v2], [], 'omitnan');
xlo = min([v1; v2], [], 'omitnan');
if isempty(xhi) || ~isfinite(xhi); xhi =  1; end
if isempty(xlo) || ~isfinite(xlo); xlo = -1; end
sp  = 0.15 * max([abs(xlo), abs(xhi), 1]);
xlim([min(xlo,0) - sp, max(xhi,0) + sp]);

% the median % change printed at the end of each bar
labelFS  = max(7, tickFS - 2);
pad      = 0.012 * diff(xlim);          % gap between the bar tip and the text
for i = 1:nn
    putLabel(v1(i), 0, i, pad, labelFS);
    putLabel(v2(i), 0, i, pad, labelFS);
end

xlabel('% change','FontSize', 35,'FontWeight','bold');

yl = ylabel('Enzyme','FontSize', 35,'FontWeight','bold');
yl.Units = 'normalized';
yl.Position = [-0.06, 0.55, 0];

lgd = legend([hb1 hb2], ...
             {'N-fixation efficiency (gC g^{-1}N)', ...
              'N-fixation rate (mM s^{-1})'}, ...
              'Location', 'southeast', ...
              'FontSize', 28, ...
              'FontWeight', 'bold');
lgd.Box = 'off';
lgd.ItemTokenSize = [18, 8];



%% Figure S19
load('../data/pw_kparam_results.mat');
nCombo    = numel(comboKeys);


% tornado plot, % change in gcgn and vNout per enzyme pair
gsort = gcgn_median_change;
gsort(isnan(gsort)) = +Inf;             % push empty pairs to the very bottom
[~, order] = sort(gsort, 'ascend');

v_g = gcgn_median_change(order);
v_v = vNout_median_change(order);
q_glo = gcgn_p5(order);    q_ghi = gcgn_p95(order);
q_vlo = vNout_p5(order);   q_vhi = vNout_p95(order);
labels_ordered = pairLabels(order);

e_gneg = v_g - q_glo;   e_gpos = q_ghi - v_g;
e_vneg = v_v - q_vlo;   e_vpos = q_vhi - v_v;

n   = nCombo;
yb  = (1:n)';

colG = [212, 127, 102] / 255;   % gcgn  (efficiency)
colV = [ 77, 116, 140] / 255;   % vNout (rate)

figure('Position', [200 60 950 1500], 'Color', 'w');
ax = gca;
ax.Color    = [240 240 240]/255;
ax.FontSize = 15;
box(ax, 'on'); ax.LineWidth = 3;
hold on;

hb1 = barh(yb, v_g, 0.6, 'FaceColor', colG);
hb2 = barh(yb, v_v, 0.6, 'FaceColor', colV);

% white separators between pair rows
for i = 0.5 : 1 : n+0.5
    yline(i, 'w-', 'LineWidth', 1.6);
end

% 5th-95th percentile band (Step 3) drawn on top of the bars
errorbar(v_g, yb, e_gneg, e_gpos, 'horizontal', 'LineStyle','none', ...
    'Color',[0.15 0.15 0.15], 'LineWidth',1.4, 'CapSize',6);
errorbar(v_v, yb, e_vneg, e_vpos, 'horizontal', 'LineStyle','none', ...
    'Color',[0.15 0.15 0.15], 'LineWidth',1.4, 'CapSize',6);

set(gca, 'YTick', 1:n, 'YTickLabel', labels_ordered, 'YDir', 'reverse', ...
    'FontSize', 23);
xline(0, 'k', 'LineWidth', 1.5);

% axis fits the whiskers (5th / 95th percentiles)
xmax = max(abs([q_glo; q_ghi; q_vlo; q_vhi]), [], 'omitnan');
if isempty(xmax) || ~isfinite(xmax) || xmax == 0; xmax = 1; end
xlim([-1.2*xmax, 1.2*xmax]);

xlabel('% change', 'FontSize', 40, 'FontWeight', 'bold');

yl = ylabel('Pairwise enzyme change', 'FontSize', 40, 'FontWeight', 'bold');
yl.Units = 'normalized';
yl.Position = [-0.03, 0.55, 0];

lgd = legend([hb1 hb2], ...
    {'N-fixation efficiency (gC g^{-1}N)', 'N-fixation rate (mM s^{-1})'}, ...
    'Location', 'southeast', 'FontSize', 25, 'FontWeight', 'bold');
lgd.Box = 'off';
lgd.ItemTokenSize = [18, 8];



%% Figure S20
alpha=0.5;
load('../data/Enames.mat');

GS_idx = find(strcmp(param_names, 'GS'));
GOGAT_idx = find(strcmp(param_names, 'GOGAT'));
v_GS=newdata.v(GS_idx,:);
v_GOGAT=newdata.v(GOGAT_idx,:);
v_netGSGOGAT=(v_GS-v_GOGAT)';

Vmaxvalues = newdata.Vmax;
Vmaxvalues([19,25],:)=[];

[prcc, prcc_sign,sign_indices,sign_param,sign_prccvalues] = PRCC_modified(Vmaxvalues', v_netGSGOGAT, Enames, alpha);

% Calculate confidence interval
param_num=length(Enames);
samp_num=size(newdata.gCgN,2);
z_score = norminv(1 - alpha/2);  
ci_lb = zeros(1, param_num);
ci_ub = zeros(1, param_num);

for i = 1:param_num
    r = prcc(1, i);
    
    % Fisher z-transform
    z = 0.5 * log((1 + r) / (1 - r));
    se = 1 / sqrt(samp_num - param_num - 3);
    
    % Confidence interval in z-space
    z_low = z - z_score * se;
    z_up = z + z_score * se;
    
    % Back to r-space
    ci_lb(1, i) = tanh(z_low);
    ci_ub(1, i) = tanh(z_up);
end

ciBCadata=cell(4,param_num+1);
ciBCadata(1,2:end)=Enames;
ciBCadata(2:end,1)={'prcc','lb','ub'};
ciBCadata(2,2:end)=num2cell(prcc);
ciBCadata(3,2:end)=num2cell(ci_lb);
ciBCadata(4,2:end)=num2cell(ci_ub);

n = numel(Enames);
x = 1:n;  % bar locations
y = prcc(:)'; % ensure row vector of length N

% extract CI bounds from your cell array
lowerBound = cell2mat(ciBCadata(3,2:end));  
upperBound = cell2mat(ciBCadata(4,2:end));

% compute error distances
lowerErr = y - lowerBound;            
upperErr = upperBound - y;            

figure('Position',[100 100 1400 1000]);
hold on;


% 1) histogram-style bars
h = bar(x, y, 0.7, ...
    'FaceColor', [152, 180, 206]/255, ...
    'EdgeColor','none');

% 2) vertical errorbars for 95% CI
errorbar(x, y, lowerErr, upperErr, ...
    'k', ...             % black whiskers
    'LineStyle','none', ...
    'LineWidth',1.5, ...
    'CapSize',8);

% 3) filled dot at the point estimate
scatter(x, y, 20, 	[0,0,0], 'filled');

% 4) build star labels
sig = strings(1,n);
sig(prcc_sign < 0.001) = "***";
sig(prcc_sign < 0.01  & prcc_sign >= 0.001) = "**";
sig(prcc_sign < 0.05  & prcc_sign >= 0.01) = "*";
% (cells with P>=0.05 remain "" and won’t be shown)

% 5) compute where to put them
ypos=nan(1,n);
for i = 1:n
    if lowerBound(i) > 0
        % CI is entirely > 0 → put stars above the upper bound
        ypos(i) = upperBound(i) + 0.001;
    else
        % CI crosses or sits below zero → put stars below the lower bound
        ypos(i) = lowerBound(i) - 0.08;
    end
end

% 6) annotate
for i = 1:n
    if sig(i) ~= ""
        text(x(i), ypos(i), sig(i), ...
             'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom', ...
             'FontSize',15, ...
             'Color','k');
    end
end


% 4) polish axes
ylim([-1, 1]);
xticks(x);
xticklabels(Enames);

ax = gca; % get current axes (parent of heatmap)

annotation('textbox', [0.11 0.15 0.5 0.05], ...   % [x y width height]
    'String', 'PRCC of Vmax to net flux in GS/GOGAT', ...
    'HorizontalAlignment', 'center', ...
    'EdgeColor', 'none', ...
    'FontSize', 21, 'Rotation', 90,...
    'FontWeight', 'bold');

ax.XAxis.FontSize = 15;   % set x-axis tick label font size

box on;
grid off;
hold off;






%% Figure S21
val_range=linspace(0, 5, 30); 

load('../data/changeFBAratio_result.mat');

figure('Units', 'pixels', 'Position', [200, 200, 900, 400]);
t = tiledlayout(1, 2, ...
    'TileSpacing', 'loose', ...
    'Padding', 'loose');  % Removes space between and around plots

ax1=nexttile;
p1=plot(ax1, val_range,gCgN, 'k-', 'LineWidth', 1.5);
xlim(ax1, [0,6]);
ylim(ax1, [0,20]);
xline(border1,'--k','LineWidth',1.5);  % horizontal line at y=0
xline(border2,'--k','LineWidth',1.5);  % horizontal line at y=0

xlabel(ax1, '$V_{MAL} / V_{NH4}$', 'Interpreter', 'latex', 'FontSize', 14, 'Fontweight','bold');
ylabel(ax1, '$gCgN^{-1}$', 'Interpreter', 'latex', 'FontSize', 14, 'Fontweight','bold');
hold on;  % Allow overlaying another plot
plot(transportratio, gCgN_median, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');  % Black filled circle
text(transportratio, gCgN_median, '   Control (1.36, 4.17)', 'FontSize', 14);  % Optional label slightly above the point
hold off;

ax2=nexttile;
p2=plot(ax2, val_range,Nfixrate, 'k-', 'LineWidth', 1.5);
xlim(ax2, [0,6]);
ylim(ax2, [0,2.5]);
xline(border1,'--k','LineWidth',1.5);  % horizontal line at y=0
xline(border2,'--k','LineWidth',1.5);  % horizontal line at y=0

xlabel('$V_{MAL} / V_{NH4}$', 'Interpreter', 'latex', 'FontSize', 14, 'Fontweight','bold');
ylabel('N-fixation rate(mM/s)', 'Interpreter', 'latex', 'FontSize', 14, 'Fontweight','bold');
hold on;  % Allow overlaying another plot
plot(transportratio, Nfix_median, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');  % Black filled circle
text(transportratio, Nfix_median, '   Control (1.36, 1.14)', 'FontSize', 14);  % Optional label slightly above the point
hold off;

labels = {'A.', 'B.'};
axes_handles = [ax1, ax2];
for i = 1:length(axes_handles)
    ax = axes_handles(i);
    xlims = xlim(ax);
    ylims = ylim(ax);
    
    % Position text near top-left, with some padding
    text(ax, xlims(1) - 0.15*range(xlims), ...
             ylims(2) - 0.015*range(ylims), ...
             labels{i}, ...
             'FontWeight', 'bold', ...
             'FontSize', 12, ...
             'Fontweight','bold',...
             'HorizontalAlignment', 'left', ...
             'VerticalAlignment', 'top');
end

width = 9;  
height = 4; 
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width height]);           
set(gcf, 'PaperPosition', [0 0 width height]);    
set(gcf, 'PaperPositionMode', 'manual');           







%% Functions
function z=gcgn(v_matrix)
    %gCgN
    CO2fluxout = v_matrix(15,:) + v_matrix(12,:) + v_matrix(30,:) + v_matrix(31,:) + ...
                 v_matrix(24,:).*4 + v_matrix(33,:).*4 + v_matrix(35,:).*4 + v_matrix(61,:).*2 + v_matrix(67,:);
    CO2fluxin = v_matrix(54,:) + v_matrix(22,:);
    CO2efflux = CO2fluxout - CO2fluxin;
    NH4fluxout = v_matrix(68,:);
    z = 12 * CO2efflux ./ (14 * 4 * NH4fluxout);
end

function [efficient_nod, inefficient_nod] = efficient_identify2(newdata)
    efficient_indices = find(newdata.gCgN(1,:) < 4);
    inefficient_indices = find(newdata.gCgN(1,:) > 6);
    
    efficient_nod=struct();
    efficient_nod.x=newdata.x(:,efficient_indices);
    efficient_nod.v=newdata.v(:,efficient_indices);
    efficient_nod.Vmax=newdata.Vmax(:,efficient_indices);
    efficient_nod.gCgN=newdata.gCgN(:,efficient_indices);
    
    inefficient_nod=struct();
    inefficient_nod.x=newdata.x(:, inefficient_indices);
    inefficient_nod.v=newdata.v(:, inefficient_indices);
    inefficient_nod.Vmax=newdata.Vmax(:, inefficient_indices);
    inefficient_nod.gCgN=newdata.gCgN(:, inefficient_indices);
end

function data = datastruct(data0)
    data = struct();
    data.x0_1_ATPcons = data0.v([1, 3, 32, 38, 21, 45, 47, 53, 55],:);
    data.x0_1_ATPsyn = data0.v([11, 8], :);
    data.x0_1_NADPHcons = data0.v([23, 39], :);
    data.x0_1_NADPHsyn = data0.v([7, 12, 13, 15, 30, 31, 40, 60, 64, 42, 52, 50], :);
    data.x0_1_QH2syn = data0.v([34], :);
    [x1_NetATPcons, x1_NetNADPH_reductant, x1_NetNADPHcons, x1_QH2_total, x1_NetATP_mM_s, x1_NetATP_mmol_day] = EnergyCal(data.x0_1_ATPcons(:,2:end),...
        data.x0_1_ATPsyn(:,2:end), data.x0_1_NADPHcons(:,2:end), ...
        data.x0_1_NADPHsyn(:,2:end), data.x0_1_QH2syn(:,2:end));
    data.NetATPcons = x1_NetATPcons;
    data.NetATP_reductant = x1_NetNADPH_reductant;
    data.NetATP_mM_s = x1_NetATP_mM_s;
    data.NetNADPHcons = x1_NetNADPHcons;
    data.QH2_total = x1_QH2_total;
    data.NetATP_mmol_day = x1_NetATP_mmol_day;
end

function [NetATPcons, NetNADPH_reductant, NetNADPHcons, QH2_total, NetATP_mM_s, NetATP_mmol_day] = EnergyCal(ATPcons, ATPsyn, NADPHcons, NADPHsyn, QH2syn)
    if iscell(ATPcons),   ATPcons   = cell2mat(ATPcons);   end
    if iscell(ATPsyn),    ATPsyn    = cell2mat(ATPsyn);    end
    if iscell(NADPHcons), NADPHcons = cell2mat(NADPHcons); end
    if iscell(NADPHsyn),  NADPHsyn  = cell2mat(NADPHsyn);  end
    if iscell(QH2syn),    QH2syn    = cell2mat(QH2syn);    end

    ATPcons_total = sum(ATPcons);
    ATPsyn_total = sum(ATPsyn);
    NetATPcons = ATPsyn_total - ATPcons_total;

    NADPHcons_total = sum(NADPHcons);
    NADPHsyn_total = sum(NADPHsyn);
    NetNADPHcons = NADPHsyn_total - NADPHcons_total;

    if isvector(QH2syn)
        QH2_total = QH2syn;
    else
        QH2_total = sum(QH2syn);
    end

    conversion_efficiency_NADHtoATP = 2.5;

    NetNADPH_reductant = NetNADPHcons + 2 .* QH2_total;
    NetATP_mM_s = NetNADPH_reductant .* conversion_efficiency_NADHtoATP + NetATPcons;%mM/s

    cellvolume = 12.4 * (10^-14)*1000;%L
    
    NetATP = NetATP_mM_s .* cellvolume;%mmol/s
    NetATP_mmol_day = NetATP .* 3600 .* 24;%mmol/day
end

function putLabel(val, err, row, pad, fs)
    % print val to 2 decimals just beyond the outer end of its error bar, on the
    % side the bar points to, so the text never sits on top of the bar
    if isnan(val)
        return
    end
    if isnan(err)
        err = 0;
    end
    if val < 0
        x     = val - err - pad;
        align = 'right';
    else
        x     = val + err + pad;
        align = 'left';
    end
    text(x, row, sprintf('%.2f', val), ...
        'HorizontalAlignment', align, 'VerticalAlignment', 'middle', ...
        'FontSize', fs, 'Color', [0.15 0.15 0.15]);
end

function t = tcrit95(n)
    % two-sided 95% t critical value; falls back to the normal
    % approximation when the Statistics Toolbox is unavailable
    if n < 2
        t = NaN;
    elseif exist('tinv','file') == 2
        t = tinv(0.975, n-1);
    else
        t = 1.96;
    end
end

function y = pctl(x, p)
    x = sort(x(~isnan(x)));
    n = numel(x);
    if n == 0
        y = NaN;
    elseif n == 1
        y = x(1);
    else
        pos = 100 * ((1:n) - 0.5) / n;
        y = interp1(pos, x, p, 'linear');
        if p <= pos(1),   y = x(1);   end
        if p >= pos(end), y = x(end); end
    end
end

function lab = enzShort(vmaxField)
    lab = erase(vmaxField, '_Vmax');
    if strcmp(lab, 'SourceMAL') 
        lab = 'MALout';
    end
end

