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

%% Figure S3
all = Vmax_matrix';
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

%% Figure S4
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


%% Figure S5
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


%% Figure S6
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

%% Figure S7
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


%% Figure S14
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

%% Figure S15
val_range=linspace(0, 5, 30); 

load('../data/changeFBAratio_result.mat');

%%
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


%% Figure S16
load('../data/result_Vmaxdown.mat');
load('../data/result_Vmaxup.mat');

fieldNames = fieldnames(Vmax);   
met_names = newVmax_stats_table{:, 1};        
x1_Vmax_repre = newVmax_stats_table{:,4};
Vmax1=Vmax;

for i = 1:length(met_names)
    idx = find(strcmp(met_names{i}, fieldNames)); 
    if ~isempty(idx)
        Vmax1.(fieldNames{idx}) = x1_Vmax_repre(i);   
    end
end

[ T1, X1, FLUXss1 ] = nodrunSimulation(Vmax1,tspan,x0,options);

gCgN_x1_original = gcgn(FLUXss1');
Nfix_x1_original = FLUXss1(end) .* 4;

n_E=numel(result_Vmaxup.E);
gcgn_change_up=NaN(n_E,1);
gcgn_change_down=NaN(n_E,1);

for i=1:n_E
    gcgn_change_up(i) = (result_Vmaxup.gcgn(i) - gCgN_x1_original) ./ gCgN_x1_original .* 100;
    gcgn_change_down(i) = (result_Vmaxdown.gcgn(i) - gCgN_x1_original) ./ gCgN_x1_original .* 100;
end

Nfix_change_up=NaN(n_E,1);
Nfix_change_down=NaN(n_E,1);

for i=1:n_E
    Nfix_change_up(i) = (result_Vmaxup.Nfix(i)-Nfix_x1_original) ./ Nfix_x1_original .* 100;
    Nfix_change_down(i) = (result_Vmaxdown.Nfix(i)-Nfix_x1_original) ./ Nfix_x1_original .* 100;

end

% plot
threshold = 0.001;  % Define cutoff for sensitivity

gcgnchange_idx = (abs(gcgn_change_up) >= threshold | abs(gcgn_change_down) >= threshold);  
gcgnchange = gcgn_change_up(gcgnchange_idx, :);  % Keep only rows where values >= threshold

gcgnchange_up=gcgn_change_up(gcgnchange_idx,:);
gcgnchange_down=gcgn_change_down(gcgnchange_idx,:);
enzyme_left=result_Vmaxdown.E(gcgnchange_idx,:);

n=numel(gcgnchange);

color1 = [
    0.85, 0.6, 0.25
];

color2 = color1 + (1 - color1) .* 0.5;

y_pos = 1:n;  % Y positions for each parameter
bar_width = 0.35;  % Width of each bar (adjust as needed)
kups_vals=gcgnchange_up;
kdowns_vals=gcgnchange_down;

figure('Position', [100, 100, 600, 900]); 
% Set gray background
ax = gca;
ax.Color = [0.9 0.9 0.9];   % light gray background

hold on;
% Plot increase bars (right side)
bar1=barh(y_pos, kups_vals, bar_width, 'FaceColor', 'flat', 'EdgeColor', 'none');
% Plot decrease bars (left side)
bar2=barh(y_pos, kdowns_vals, bar_width, 'FaceColor', 'flat', 'EdgeColor', 'none');

% White horizontal separator lines
for i = 0.5 : 1 : n+0.5
    yline(i, 'w-', 'LineWidth', 0.5);
end

bar_colors = zeros(length(y_pos), 3);  % initialize RGB array
bar1.CData = color1;
bar2.CData = color2;

ynames = enzyme_left;

xline(0, 'k-', 'LineWidth', 0.5);
set(gca, 'Ydir','reverse');
set(gca, 'YTick', y_pos, 'YTickLabel', ynames,'FontSize', 15);

% annotation('textbox', [0.37 0.03 0.4 0.05], 'String', '\textbf{Percent Change in $\mathbf{gCg}^{-1}\mathbf{N}$ \%}', ...
%     'Interpreter', 'latex', 'HorizontalAlignment', 'center',...
%     'EdgeColor', 'none', 'FontSize', 16, 'FontWeight', 'bold');
annotation('textbox', [0.37 0.03 0.4 0.05], ...
    'String', '\textbf{Percent Change in $\mathbf{gC\,g^{-1}\,N}$ (\%)}', ...
    'Interpreter', 'latex', ...
    'HorizontalAlignment', 'center', ...
    'EdgeColor', 'none', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

color_all=[color1; color2];
hold on
h(1) = bar(nan, nan, 'FaceColor', color_all(1,:), 'EdgeColor', 'none');
h(2) = bar(nan, nan, 'FaceColor', color_all(2,:), 'EdgeColor', 'none');

legend_handle = legend({'increase 3-fold', 'decrease 3-fold'});
set(legend_handle, 'Location', 'northeast','fontsize',15); 


%% Figure S17
threshold = 0.001;  % Define cutoff for sensitivity

Nfix_change_idx = (abs(Nfix_change_up) >= threshold | abs(Nfix_change_down) >= threshold);  
Nfixchange = Nfix_change_up(Nfix_change_idx, :);  % Keep only rows where values >= threshold

Nfixchange_up=Nfix_change_up(Nfix_change_idx,:);
Nfixchange_down=Nfix_change_down(Nfix_change_idx,:);
enzyme_left=result_Vmaxdown.E(Nfix_change_idx,:);

n=numel(Nfixchange);

color1 = [
    0.2549, 0.5961, 0.6745
];

color2 = color1 + (1 - color1) .* 0.5;

y_pos = 2:n;  % Y positions for each parameter
bar_width = 0.35;  % Width of each bar (adjust as needed)
kups_vals=Nfixchange_up(2:end);
kdowns_vals=Nfixchange_down(2:end);


figure('Position', [100, 100, 600, 900]);  
% Set gray background
ax = gca;
ax.Color = [0.9 0.9 0.9];   % light gray background

hold on;
% Plot increase bars (right side)
bar1=barh(y_pos, kups_vals, bar_width, 'FaceColor', 'flat', 'EdgeColor', 'none');
% Plot decrease bars (left side)
bar2=barh(y_pos, kdowns_vals, bar_width, 'FaceColor', 'flat', 'EdgeColor', 'none');

% White horizontal separator lines
for i = 0.5 : 1 : n+0.5
    yline(i, 'w-', 'LineWidth', 0.5);
end

bar_colors = zeros(length(y_pos), 3);  % initialize RGB array
bar1.CData = color1;
bar2.CData = color2;

ynames = enzyme_left(2:end);


xline(0, 'k-', 'LineWidth', 0.5);
set(gca, 'Ydir','reverse');
set(gca, 'YTick', y_pos, 'YTickLabel', ynames,'FontSize', 15);

annotation('textbox', [0.37 0.03 0.4 0.05], 'String', '\textbf{Percent Change in N-fixation rate (\%)}', ...
    'Interpreter', 'latex', 'HorizontalAlignment', 'center',...
    'EdgeColor', 'none', 'FontSize', 16, 'FontWeight', 'bold');


color_all=[color1; color2];
hold on
h(1) = bar(nan, nan, 'FaceColor', color_all(1,:), 'EdgeColor', 'none');
h(2) = bar(nan, nan, 'FaceColor', color_all(2,:), 'EdgeColor', 'none');

legend_handle = legend({'increase 3-fold', 'decrease 3-fold'});
set(legend_handle, 'Location', 'northwest','fontsize',15);  


%% Figure S18
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

