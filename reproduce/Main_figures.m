%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains supporting information of the following publication:
%
% Mechanistic modeling of the determinate nodule metabolism reveals enzymatic 
% influences on improving nitrogen fixation efficiency
%
% by Rourou Ji, Joshua A.M. Kaste and Megan L. Matthews
%
% This MATLAB script allows to reproduce Figure 2, 4, 6, and 7
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

%% Plot Figure 2
figure('Units', 'pixels', 'Position', [200, 200, 850, 800]);
t = tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'none');  % Removes space between and around plots

% Figure 1
% subplot(2,2,1);
ax1=nexttile;
h1=histogram(ax1, gCgN, 'BinWidth', 1, 'FaceColor', [0.85, 0.6, 0.25], 'Normalization', 'count'); % blue-green
% Add boundary lines
xline(ax1, 2, '--k', 'LB=2', 'LabelHorizontalAlignment', 'right', 'FontSize', 15);
xline(ax1, 8, '--k', 'UB=8', 'LabelHorizontalAlignment', 'right', 'FontSize', 15);
xlim(ax1, [1,200]);

% >>> remove the old ylim or set it dynamically <<<
ymax1 = max(h1.Values);
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
h2=histogram(ax2, Nfix, 'BinWidth', 0.1, 'FaceColor', [0.2549, 0.5961, 0.6745], 'Normalization', 'count'); %orange
xlim(ax2, [1e-1,12]);

% >>> remove the old ylim or set it dynamically <<<
ymax2 = max(h2.Values);
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
h3=histogram(ax3, newdata.gCgN, 'BinWidth', 1, 'FaceColor', [0.85, 0.6, 0.25], 'Normalization', 'count'); 
% Add boundary lines
xline(ax3, 2, '--k', 'LB=2', 'LabelHorizontalAlignment', 'right', 'FontSize', 15);
xline(ax3, 8, '--k', 'UB=8', 'LabelHorizontalAlignment', 'right', 'FontSize', 15);
xlim(ax3, [1,200]);

% >>> remove the old ylim or set it dynamically <<<
ymax3 = max(h3.Values);
if ~isempty(ymax3) && isfinite(ymax3)
    ylim(ax3, [0, ymax1*1.15]);  % add a little headroom
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
h4=histogram(ax4, newdata.N_ExportRate, 'BinWidth', 0.1, 'FaceColor', [0.2549, 0.5961, 0.6745], 'Normalization', 'count'); 
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


%% Figure 4A
% load data
load('../data/LSA_ALTNout_50%.mat');
load('../data/LSA_gcgn_50%.mat');

LSAparam=LSA_ALTNout_data.Var1;

% change other param names
idx_MALexport=find(strcmp(LSAparam, 'SourceMAL_KmMAL'));
idx_SUCexport=find(strcmp(LSAparam, 'SourceSUC_KmSUC'));
idx_FUMexport=find(strcmp(LSAparam, 'SourceFUM_KmFUM'));

LSAparam(idx_MALexport)={'MAL export_KmMAL'};
LSAparam(idx_SUCexport)={'SUC export_KmSUC'};
LSAparam(idx_FUMexport)={'FUM export_KmFUM'};


LSA_ALTNout_data.Var1=LSAparam;
LSA_gcgn_data.Var1=LSAparam;


% plot
threshold = 0.001;  % Define cutoff for sensitivity
Kup_gcgn = cell2mat(LSA_gcgn_data.Kupchange);  
Kdown_gcgn = cell2mat(LSA_gcgn_data.Kdownchange);  
keep_idx_gcgn = (abs(Kup_gcgn) >= threshold | abs(Kdown_gcgn) >= threshold);  
filtered_gcgn_SA = LSA_gcgn_data(keep_idx_gcgn, :);  % Keep only rows where values >= threshold

filtered_gcgn_SA.Var1 = replace(filtered_gcgn_SA.Var1, '_', ' ');

n=numel(filtered_gcgn_SA.Var1);

color1 = [
     0.85, 0.6, 0.25;%orange
];

color2 = color1 + (1 - color1) .* 0.5;

y_pos = 1:length(filtered_gcgn_SA.Var1);  % Y positions for each parameter
bar_width = 0.8;%0.35;  % Width of each bar (adjust as needed)
kups_vals=cell2mat(filtered_gcgn_SA.Kupchange).*100;
kdowns_vals=cell2mat(filtered_gcgn_SA.Kdownchange).*100;

figure('Position', [100, 100, 600, 900]); 

% Set gray background
ax = gca;
ax.Color = [240 240 240]/255;   % light gray background
hold on;

% Plot increase bars (right side)
bar1=barh(y_pos, kups_vals, bar_width, 'FaceColor', 'flat', 'EdgeColor', 'none');
% Plot decrease bars (left side)
bar2=barh(y_pos, kdowns_vals, bar_width, 'FaceColor', 'flat', 'EdgeColor', 'none');

ylim([0.5, n+0.5])

% White horizontal separator lines
for i = 0.5 : 1 : n+0.5
    yline(i,'Color', [170 170 170]/255,'LineWidth', 0.5);%'w','-',
end

% Apply colors to both increase and decrease bars
bar1.CData = color1;
bar2.CData = color2;

yline(3.5, 'k-', 'LineWidth', 0.8);  
y_min = ylim(); 
text(3.5, 2.5, 'Transport', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 16,'FontWeight','bold');

yline(12.5, 'k-', 'LineWidth', 0.8);  
y_min = ylim();  
text(3.5, 8.5, 'Glycolysis', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 16,'FontWeight','bold');

yline(28.5, 'k-', 'LineWidth', 0.8);
y_min = ylim(); 
text(3.5, 20.5, 'PP Pathway', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 16,'FontWeight','bold');

yline(34.5, 'k-', 'LineWidth', 0.8); 
y_min = ylim(); 
text(3.5, 32.5, 'TCA cycle', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 16,'FontWeight','bold');

y_min = ylim();
text(3.5, 35.5, 'Glyoxylate cycle', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 16,'FontWeight','bold');

ynames_gcgn = filtered_gcgn_SA.Var1;
unknownKindex=[1,2,3,7,11,12,16,20,31,33,34,35];

for i=1:length(ynames_gcgn)
    ynames_gcgn{i} = strrep(ynames_gcgn{i}, '__', ' ');
    if ismember(i, unknownKindex)
        ynames_gcgn{i} = ['\bf', ynames_gcgn{i}];
    end
end

xline(0, 'k-', 'LineWidth', 0.5);
set(gca, 'YTick', y_pos, 'YTickLabel', ynames_gcgn,'FontSize', 15);
xlim([-5, 5]); % Adjust x-axis limits based on data range

annotation('textbox', [0.37 0.03 0.4 0.05], 'String', '\textbf{Change in nitrogen fixation efficiency (\%)}', ...
    'Interpreter', 'latex', 'HorizontalAlignment', 'center',...
    'EdgeColor', 'none', 'FontSize', 18, 'FontWeight', 'bold');


color_all=[color1; color2];
hold on
h(1) = bar(nan, nan, 'FaceColor', color_all(1,:), 'EdgeColor', 'none');
h(2) = bar(nan, nan, 'FaceColor', color_all(2,:), 'EdgeColor', 'none');

legend_handle = legend({'+50%', '-50%'});
set(legend_handle, 'Position', [0.78, 0.6, 0.1, 0.04],'fontsize',17);  % Bottom-right corner


%% Figure 4B
threshold = 0.001;  % Define cutoff for sensitivity

Kup_vALTNout = cell2mat(LSA_ALTNout_data.Kupchange);  
Kdown_vALTNout = cell2mat(LSA_ALTNout_data.Kdownchange);  
keep_idx_vALTNout = (abs(Kup_vALTNout) >= threshold | abs(Kdown_vALTNout) >= threshold);  
filtered_vALTNout_SA = LSA_ALTNout_data(keep_idx_vALTNout, :);  % Keep only rows where values >= threshold

filtered_vALTNout_SA.Var1 = replace(filtered_vALTNout_SA.Var1, '_', '__');

color1 = [
    0.2549, 0.5961, 0.6745
];

color2 = color1 + (1 - color1) .* 0.5;

y_pos = 1:length(filtered_vALTNout_SA.Var1);  % Y positions for each parameter
bar_width = 0.8;%0.35;  % Width of each bar (adjust as needed)
kups_vals=cell2mat(filtered_vALTNout_SA.Kupchange).*100;
kdowns_vals=cell2mat(filtered_vALTNout_SA.Kdownchange).*100;

figure('Position', [100, 100, 600, 900]); 
% Set gray background
ax = gca;
ax.Color = [240 240 240]/255;   % light gray background

hold on;
% Plot increase bars (right side)
bar1=barh(y_pos, kups_vals, bar_width, 'FaceColor', 'flat', 'EdgeColor', 'none');
% Plot decrease bars (left side)
bar2=barh(y_pos, kdowns_vals, bar_width, 'FaceColor', 'flat', 'EdgeColor', 'none');

n=numel(y_pos);
ylim([0.5, n+0.5]);% border of the y axis

% White horizontal separator lines
for i = 0.5 : 1 : n+0.5
    yline(i, 'Color', [170 170 170]/255,  'LineWidth', 0.5);
end

bar_colors = zeros(length(y_pos), 3);  % initialize RGB array
bar1.CData = color1;
bar2.CData = color2;

yline(3.5, 'k-', 'LineWidth', 0.8);  
text(3.5, 2.5, 'Transport', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 16,'FontWeight','bold');

yline(12.5, 'k-', 'LineWidth', 0.8);  
text(3.5, 8.5, 'Glycolysis', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 16,'FontWeight','bold');

yline(27.5, 'k-', 'LineWidth', 0.8);
text(3.5, 18.5, 'PP Pathway', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 16,'FontWeight','bold');

text(3.5, 30.5, 'TCA cycle', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 16,'FontWeight','bold');

ynames_vALTN = filtered_vALTNout_SA.Var1;
unknownKindex=[1,2,3,8,12,16,20,29,31,32];

for i=1:length(ynames_vALTN)
    ynames_vALTN{i} = strrep(ynames_vALTN{i}, '__', ' ');
    if ismember(i, unknownKindex)
        ynames_vALTN{i} = ['\bf', ynames_vALTN{i}];
    end
end

xline(0, 'k-', 'LineWidth', 0.5);
set(gca, 'YTick', y_pos, 'YTickLabel', ynames_vALTN,'FontSize', 15);

annotation('textbox', [0.37 0.03 0.4 0.05], 'String', '\textbf{Change in nitrogen fixation rate (\%)}', ...
    'Interpreter', 'latex', 'HorizontalAlignment', 'center',...
    'EdgeColor', 'none', 'FontSize', 18, 'FontWeight', 'bold');


legend('+50%', '-50%', 'Location', 'best');
xlim([-5, 5]); % Adjust x-axis limits based on data range

color_all=[color1; color2];
hold on
h(1) = bar(nan, nan, 'FaceColor', color_all(1,:), 'EdgeColor', 'none');
h(2) = bar(nan, nan, 'FaceColor', color_all(2,:), 'EdgeColor', 'none');

legend_handle = legend({'+50%', '-50%'});
set(legend_handle, 'Position', [0.78, 0.6, 0.1, 0.04],'fontsize',17);  % Bottom-right corner


%% Figure 6
load('../data/Enames.mat');

newdata_Vmax=newdata.Vmax;
newdata_Vmax(19,:)=[];

alpha=0.5;
[prcc, prcc_sign,sign_indices,sign_param,sign_prccvalues] = PRCC_modified(newdata_Vmax', newdata.v', Enames, alpha);


%% extract sig rxns
rxn_sig_pval=prcc_sign(sign_indices);
rxn_sig = [sign_param',num2cell(sign_prccvalues),num2cell(rxn_sig_pval)];

%% heatmap
xvalues = reaction_names;%x axis of heatmap
yvalues = Enames;%y axis of heatmap
%check if size is same
[m, n] = size(prcc');
if length(xvalues) ~= n
    error('The number of metabolites (xvalues) must equal the number of columns in prcc (%d columns).', n);
end

if length(yvalues) ~= m
    error('The number of reactions (yvalues) must equal the number of rows in prcc (%d rows).', m);
end

n = 256;  % number of colors in the colormap (higher = smoother)
% Color anchors
blue  = [64, 128, 148] / 255;
white = [1, 1, 1];
red   = [198, 91, 63] / 255;
% Allocate steps
n1 = round(n * 0.46);  % blue → white
n2 = round(n * 0.08);  % white plateau (repeated white)
n3 = round(n * 0.46);     % white → red
% Interpolate
cmap_bluetowhite = [ ...
    linspace(red(1), white(1), n1)', ...
    linspace(red(2), white(2), n1)', ...
    linspace(red(3), white(3), n1)'];
cmap_white = repmat(white, n2, 1);  % plateau region
cmap_whitetored = [ ...
    linspace(white(1), blue(1), n3)', ...
    linspace(white(2), blue(2), n3)', ...
    linspace(white(3), blue(3), n3)'];
% Combine
cmap = [cmap_bluetowhite; cmap_white; cmap_whitetored];

data = prcc';

figure();
h=heatmap(xvalues,yvalues,data);
h.Colormap = cmap;

ax = gca; % get current axes (parent of heatmap)
annotation('textbox', [0.45 0.015 0.13 0.025], 'String', 'Flux (mM s^{-1})', ...
    'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'FontSize', 20, 'FontWeight','bold');

annotation('textbox', [0.07 0.4 0.2 0.05], 'String', 'Vmax (mM s^{-1})', ...
    'HorizontalAlignment', 'center', 'EdgeColor', 'none', ...
    'FontSize', 20, 'Rotation', 90, 'FontWeight','bold');


% Improve readability
h.GridVisible = 'on';%grid of heatmap
h.ColorbarVisible = 'on';%colorbar
% Display numeric values clearly inside cells:
h.CellLabelFormat = '%.3f';%digits in heatmap cells
h.FontSize = 14;%for labels and titles

annotation('textbox',[0.55 0.9 0.7 0.05], ...
    'String','PRCC', ...
    'HorizontalAlignment','center', 'EdgeColor','none', ...
    'FontSize',17);

grid off;



%% Figure 7
%% Figure 7A
alpha=0.05;

[prcc2, prcc_sign2,sign_indices2,sign_param2,sign_prccvalues2] = PRCC_modified(newdata_Vmax', newdata.gCgN', Enames, alpha);

%% Calculate confidence interval
%Fisher z-transformation
param_num=length(Enames);
samp_num=size(newdata.x,2);
z_score = norminv(1 - alpha/2);  % e.g., 1.96 for 95% CI
ci_lb2 = zeros(1, param_num);
ci_ub2 = zeros(1, param_num);

for i = 1:param_num
    r = prcc2(1, i);
    
    % Fisher z-transform
    z = 0.5 * log((1 + r) / (1 - r));
    se = 1 / sqrt(samp_num - param_num - 3);
    
    % Confidence interval in z-space
    z_low = z - z_score * se;
    z_up = z + z_score * se;
    
    % Back to r-space
    ci_lb2(1, i) = tanh(z_low);
    ci_ub2(1, i) = tanh(z_up);
end

ciBCadata2=cell(4,param_num+1);
ciBCadata2(1,2:end)=Enames;
ciBCadata2(2:end,1)={'prcc','lb','ub'};
ciBCadata2(2,2:end)=num2cell(prcc2);
ciBCadata2(3,2:end)=num2cell(ci_lb2);
ciBCadata2(4,2:end)=num2cell(ci_ub2);


%% gcgn PRCC with CI values
n = numel(Enames);
x = 1:n;  % bar locations
y = prcc2(:)'; % ensure row vector of length N

% extract CI bounds from your cell array
lowerBound = cell2mat(ciBCadata2(3,2:end));  
upperBound = cell2mat(ciBCadata2(4,2:end));

% compute error distances
lowerErr = y - lowerBound;            
upperErr = upperBound - y;            

figure('Position',[100 100 1400 1000]);
hold on;

% light gray background
ax = gca;
ax.Color = [240 240 240]/255;   % light gray background

% vertical separators
for i = 0.5 : 1 : n+0.5
    xline(i, 'w-', 'LineWidth',1.6);
end

% 1) histogram-style bars
hb = bar(x, y, 0.7, ...
    'FaceColor',  [238, 192, 163]/255, ...
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
sig(prcc_sign2 < 0.001) = "***";
sig(prcc_sign2 < 0.01  & prcc_sign2 >= 0.001) = "**";
sig(prcc_sign2 < 0.05  & prcc_sign2 >= 0.01) = "*";
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
             'FontSize',12, ...
             'Color','k');
    end
end


% 4) polish axes
ylim([-1, 1]);

xticks(x);
xticklabels(Enames);
ax = gca; % get current axes (parent of heatmap)

annotation('textbox', [0.11 0.15 0.5 0.05], ...   % [x y width height]
    'String', 'PRCC of Vmax to nitrogen fixation efficiency (gC g^{-1}N)', ...
    'HorizontalAlignment', 'center', ...
    'EdgeColor', 'none', ...
    'FontSize', 21, 'Rotation', 90,...
    'FontWeight', 'bold');

ax.XAxis.FontSize = 15;   % set x-axis tick label font size
ax.YAxis.FontSize = 15;   % set x-axis tick label font size

% grid off;
box on;
hold off;

%% Figure 7B
alpha=0.05;
newdata_vALTNout=newdata.v(end,:) .* 4;

[prcc, prcc_sign,sign_indices,sign_param,sign_prccvalues] = PRCC_modified(newdata_Vmax', newdata_vALTNout', Enames, alpha);


%% Calculate confidence interval

%Fisher z-transformation
param_num=length(Enames);
samp_num=size(newdata.x,2);
z_score = norminv(1 - alpha/2);  % e.g., 1.96 for 95% CI
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


%% Nfix PRCC with CI values
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

% light gray background
ax = gca;
ax.Color = [240 240 240]/255;   % light gray background

% vertical separators
for i = 0.5 : 1 : n+0.5
    xline(i, 'w-', 'LineWidth',1.6);
end

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
    'String', 'PRCC of Vmax to nitrogen fixation rate (mM s^{-1})', ...
    'HorizontalAlignment', 'center', ...
    'EdgeColor', 'none', ...
    'FontSize', 21, 'Rotation', 90,...
    'FontWeight', 'bold');

ax.XAxis.FontSize = 15;   % set x-axis tick label font size
ax.YAxis.FontSize = 15;   % set x-axis tick label font size

box on;
grid off;
hold off;


%% Figure 8
%% load pairwise data
load('../data/gcgn_pairwise.mat');
load('../data/Nfix_pairwise.mat');

%% control results
x1_Vmax_repre=newVmax_stats_table{:,4};

% initial values and load parameter settings
fieldNames = fieldnames(Vmax);        % Get all field names from struct A
met_names = newVmax_stats_table{:, 1};                    % Column 1 of B: field names to check
Vmax1=Vmax;

for i = 1:length(met_names)
    idx = find(strcmp(met_names{i}, fieldNames)); % Find index of matching field
    if ~isempty(idx)
        Vmax1.(fieldNames{idx}) = x1_Vmax_repre(i);     % Assign new value from B to A
    end
end

[ T1, X1, FLUXss1 ] = nodrunSimulation(Vmax1,tspan,x0,options);

gCgN_x1_original = gcgn(FLUXss1');
Nfix_x1_original = FLUXss1(end) .* 4;

%% All enzymes affect gcgn & Nfix across all significant level
scaler = [0.33, 0.5, 0.75, 1, 1.5, 2, 3];

Kcat_TKT1=131.2;
Kcat_TKT2=69.05;
E_TKT=Vmax1.TKT1 ./ Kcat_TKT1;

% All enzymes affect gcgn & Nfix across all significant level
both_PGI_range=Vmax1.PGI .* scaler;
both_PFK_range=Vmax1.PFK .* scaler;
both_PBP_range=Vmax1.FBP .* scaler;
both_GDH_range=Vmax1.GDH .* scaler;
both_PYK_range=Vmax1.PYK .* scaler;
both_ZWF_range=Vmax1.ZWF .* scaler;
both_PGL_range=Vmax1.PGL .* scaler;
both_TKT_range=E_TKT .* scaler;
both_TAL_range=Vmax1.TAL .* scaler;
both_PEPC_range=Vmax1.PEPC .* scaler;
both_LPD_range=Vmax1.LPD .* scaler;
both_SDH_range=Vmax1.SDH .* scaler;
both_FUMA_range=Vmax1.FUMA .* scaler;
both_MALout_range=Vmax1.MALOut .* scaler;
both_SUCout_range=Vmax1.SUCOut .* scaler;
both_FUMout_range=Vmax1.FUMOut .* scaler;
both_FGAMS_range=Vmax1.FGAMS .* scaler;
both_IMPCH_range=Vmax1.IMPCH .* scaler;
both_IMPDH_range=Vmax1.IMPDH .* scaler;


enzNames  = {'PGI','PFK','FBP','GDH','PYK',...
    'ZWF','PGL','TKT','TAL','PEPC',...
    'LPD','SDH','FUMA','MALOut','SUCout','FUMout',...
    'FGAMS','IMPCH','IMPDH'};
enzRanges = {both_PGI_range, both_PFK_range,both_PBP_range,both_GDH_range, both_PYK_range,...
    both_ZWF_range,both_PGL_range, both_TKT_range, both_TAL_range, both_PEPC_range,...
    both_LPD_range, both_SDH_range, both_FUMA_range, ...
    both_MALout_range, both_SUCout_range, both_FUMout_range,both_FGAMS_range,...
    both_IMPCH_range, both_IMPDH_range};

%%
result_ranges1=struct();
n = numel(enzNames);
for p = 1:n-1
    for q=p+1:n
        nameP = enzNames{p};
        nameQ = enzNames{q};
        rangeP = enzRanges{p};
        rangeQ = enzRanges{q};

        if strcmp(nameP, 'TKT') || strcmp(nameQ, 'TKT')
            continue
        else
            pairKey = sprintf('%s__%s', nameP, nameQ);

            rangecell=cell(2,1);
            rangecell{1}=rangeP;
            rangecell{2}=rangeQ;
    
            result_ranges1.(pairKey)=rangecell;
        end
    end
end

result_ranges2=struct();
TKTidx=find(strcmp(enzNames,'TKT'));
for p=1:(length(enzNames))
    nameTKT='TKT';
    nameP=enzNames{p};
    rangeTKT=enzRanges{TKTidx};
    rangeP=enzRanges{p};
    if strcmp(nameP, 'TKT')
        continue
    else
        pairKey = sprintf('%s__%s', nameTKT, nameP);
        
        rangecell=cell(2,1);
        rangecell{1}=rangeTKT;
        rangecell{2}=rangeP;
        
        result_ranges2.(pairKey)=rangecell; 
    end   
end

result_ranges=result_ranges1;
fieldsA = fieldnames(result_ranges2);
for k = 1:numel(fieldsA)
    result_ranges.(fieldsA{k}) = result_ranges2.(fieldsA{k});
end


%% change the order of results

% gcgn
struct_gcgn=result_gcgn;
% Flip all matrices row-wise
fn = fieldnames(struct_gcgn);
for k = 1:numel(fn)
    struct_gcgn.(fn{k}) = flipud(struct_gcgn.(fn{k}));   % or S.(fn{k})(end:-1:1,:)
end

% Nfix
struct_Nfix=result_Nfix;
% Flip all matrices row-wise
fn = fieldnames(struct_Nfix);
for k = 1:numel(fn)
    struct_Nfix.(fn{k}) = flipud(struct_Nfix.(fn{k}));   % or S.(fn{k})(end:-1:1,:)
end

%%
C_gcgnall = struct2cell(result_gcgn);  %convert to cell
minVals_gcgnall = cellfun(@(M) min(M(:)), C_gcgnall); 

C_Nfixall = struct2cell(result_Nfix);  %convert to cell
maxVals_Nfixall = cellfun(@(M) max(M(:)), C_Nfixall); 

gcgnchanges=(gCgN_x1_original - minVals_gcgnall) ./ gCgN_x1_original .* 100;
Nfixchanges=(maxVals_Nfixall - Nfix_x1_original) ./ Nfix_x1_original .* 100;

%% convert gcgn and nfix to table
fngcgn = fieldnames(result_gcgn);
valsgcgn = struct2cell(result_gcgn);
Tgcgn = table(fngcgn, valsgcgn, minVals_gcgnall, 'VariableNames', {'FieldName', 'gcgn','min'});

fnNfix = fieldnames(result_Nfix);
valsNfix = struct2cell(result_Nfix);
TNfix = table(fnNfix, valsNfix, maxVals_Nfixall, 'VariableNames', {'FieldName', 'Nfix','min'});

%% enzymes for gcgn and Nfix 

% name for gcgn
enzyme_all={'FBP- FUMA+','FBP- GDH+','FBP- LPD-','FBP- MALout+','FBP- PEPC+','FBP- PGL-','FBP- PYK-','FBP- SDH+','FBP- TAL-','FBP- ZWF-',...
    'GDH+ FUMA+','GDH+ LPD-','GDH+ MALout+','GDH+ PEPC+','GDH+ PGL-','GDH+ PYK-','GDH+ SDH+','GDH+ TAL-','GDH+ ZWF-',...
    'PGI+ FBP+','PGI+ FUMA+','PGI+ GDH+','PGI+ LPD-','PGI+ MALout+','PGI+ PEPC+','PGI+ PFK+','PGI+ PGL-','PGI+ PYK-','PGI+ SDH+','PGI+ TAL-','PGI+ ZWF-',...
    'ZWF- FUMA+','ZWF- LPD-','ZWF- MALout+','ZWF- PEPC+','ZWF- PGL-','ZWF- SDH+','ZWF- TAL-',...
    'FUMA+ MALout+',...
    'LPD- FUMA+','LPD- MALout+','LPD- SDH+',...
    'PGL+ FUMA+','PGL+ LPD-','PGL+ MALout+','PGL+ PEPC+','PGL+ SDH+','PGL+ TAL-',...
    'PYK- FUMA+','PYK- LPD-','PYK- MALout+','PYK- PEPC+','PYK- PGL-','PYK- SDH+','PYK- TAL-','PYK- ZWF-',...
    'PEPC+ FUMA+','PEPC+ LPD-','PEPC+ MALout+','PEPC+ SDH+',...
    'PFK+ FBP-','PFK+ FUMA+','PFK+ GDH+','PFK+ LPD-','PFK+ MALout+','PFK+ PEPC+','PFK+ PGL-','PFK+ PYK-','PFK+ SDH+','PFK+ TAL-','PFK+ ZWF-',...
    'SDH+ FUMA+','SDH+ MALout+',...
    'TAL- LPD-','TAL- MALout+','TAL- PEPC+','TAL- SDH+','TAL- FUMA+',...
    'TKT- FBP-','TKT- FUMA+','TKT- GDH+','TKT- LPD-','TKT- MALout+','TKT- PEPC+','TKT- PFK+','TKT- PGI+','TKT- PGL-','TKT- PYK-','TKT- SDH+','TKT- TAL-','TKT- ZWF-'
    };

idx_pair_all  = [2,4,7,8,9,10,11,12,14,15,...
    28,32,33,34,35,36,37,39,40,...
    80,82,84,87,88,89,90,91,92,93,95,96,...
    161,165,166,167,168,169,171,...
    22,...
    43,44,48,...
    98,102,103,104,105,107,...
    109,113,114, 115,116,117,119,120,...
    56,60,61,62,...
    64,66,68,71,72,73,74,75,76,78,79,...
    122, 126,...
    137,138,139,140,133,...
    142,144,146,149,150,151,152,153,154,155,156,158,159
    ];

allidx=1:numel(fieldnames(result_ranges));

temp_gcgn=cell(numel(idx_pair_all),1);
temp_Nfix=cell(numel(idx_pair_all),1);


[tf, idx] = ismember(idx_pair_all, allidx);
temp_gcgn(:, 1) = num2cell(gcgnchanges(idx, 1));

[tf, idx] = ismember(idx_pair_all, allidx);
temp_Nfix(:, 1) = num2cell(Nfixchanges(idx, 1));

alldata=cell(numel(enzyme_all),3);
alldata(:,1)=enzyme_all;
alldata(:,2)=temp_gcgn;
alldata(:,3)=temp_Nfix;

%% check their location to find different patterns in pairs

% Assume you already have struct A and struct B in your workspace

% Assume A and B are structs where each field is a 7x7 double matrix
A=result_gcgn;
B=result_Nfix;

fieldsA = fieldnames(A);
fieldsB = fieldnames(B);
numFields = numel(fieldsA);

A_results = cell(numFields, 4);
B_results = cell(numFields, 4);

for i = 1:numFields
    % ---------- Struct A ----------
    fieldNameA = fieldsA{i};
    dataA = round(A.(fieldNameA), 4);        % round to 3 decimal places
    A.(fieldNameA) = dataA;

    % find min and its coordinates
    [minVal, idxMin] = min(dataA(:));
    [rowMin, colMin] = ind2sub(size(dataA), idxMin);

    % check if any row or column is constant
    rowConstant = any(all(diff(dataA, 1, 2) == 0, 2));  % true if any row constant
    colConstant = any(all(diff(dataA, 1, 1) == 0, 1));  % true if any column constant

    % if constant row/col found, replace coordinate with 'NA'
    if all(diff(dataA(rowMin, :)) == 0)  % row is constant
        rowCoord = "NA";
    else
        rowCoord = rowMin;
    end

    if all(diff(dataA(:, colMin)) == 0)  % column is constant
        colCoord = "NA";
    else
        colCoord = colMin;
    end

    % store result
    A_results{i, 1} = fieldNameA;
    A_results{i, 2} = rowCoord;
    A_results{i, 3} = colCoord;
    A_results{i, 4} = minVal;

    % ---------- Struct B ----------
    fieldNameB = fieldsB{i};
    dataB = round(B.(fieldNameB), 4);        % round to 3 decimal places
    B.(fieldNameB) = dataB;

    % find max and its coordinates
    [maxVal, idxMax] = max(dataB(:));
    [rowMax, colMax] = ind2sub(size(dataB), idxMax);

    % check for constant rows/columns
    if all(diff(dataB(rowMax, :)) == 0)
        rowCoordB = "NA";
    else
        rowCoordB = rowMax;
    end

    if all(diff(dataB(:, colMax)) == 0)
        colCoordB = "NA";
    else
        colCoordB = colMax;
    end

    % store result
    B_results{i, 1} = fieldNameB;
    B_results{i, 2} = rowCoordB;
    B_results{i, 3} = colCoordB;
    B_results{i, 4} = maxVal;
end


% === Filtering results without "NA" in coordinates ===

% For A_table
validA = ~strcmp(string(A_results(:,2)), "NA") & ~strcmp(string(A_results(:,3)), "NA");
A_table_filtered = cell2table(A_results(validA, :), ...
    'VariableNames', {'FieldName', 'Row_num', 'Column_num', 'MinValue'});


% For B_table
validB = ~strcmp(string(B_results(:,2)), "NA") & ~strcmp(string(B_results(:,3)), "NA");
B_table_filtered = cell2table(B_results(validB, :), ...
    'VariableNames', {'FieldName', 'Row_num', 'Column_num', 'MinValue'});


%%
% --- Assume A_table_filtered and B_table_filtered are tables ---
% Columns: FieldName | Row | Column | MinValue/MaxValue

% 1. Sort both tables by FieldName (to align)
A_sorted = sortrows(A_table_filtered, 'FieldName');
B_sorted = sortrows(B_table_filtered, 'FieldName');

% 2. Find matching field names
[commonFields, ia, ib] = intersect(A_sorted.FieldName, B_sorted.FieldName, 'stable');

% 3. Extract matching rows
A_common = A_sorted(ia, :);
B_common = B_sorted(ib, :);

% 4. Identify rows where Row or Column differ
diff_idx = (A_common.Row_num ~= B_common.Row_num) | (A_common.Column_num ~= B_common.Column_num);

% 5. Extract those differing rows
A_diff = A_common(diff_idx, :);
B_diff = B_common(diff_idx, :);

% 6. Combine results into one comparison table
diff_table = table( ...
    A_diff.FieldName, ...
    A_diff.Row_num, A_diff.Column_num, ...
    B_diff.Row_num, B_diff.Column_num, ...
    'VariableNames', {'FieldName', 'A_Row', 'A_Column', 'B_Row', 'B_Column'} );



%%
v1 = alldata(:,2);   % first vector
v2 = alldata(:,3);   % second vector

v1 = cell2mat(v1);
v2 = cell2mat(v2);

nn=numel(v1);

% mark rows where either vector has NaN
nanMask = isnan(v1) | isnan(v2);

% pick sorting variable (here v2, but could be v1)
v1_forSort = v1;
v1_forSort(nanMask) = -Inf;   % force NaN rows to bottom

[~, order] = sort(v1_forSort, 'descend');

% reorder everything
v1 = v1(order);
v2 = v2(order);

pairnames_ordered = enzyme_all(order);

% Make horizontal bar chart
figure('Position',[200 200 800 400]);

% light gray background
ax = gca;
ax.Color = [240 240 240]/255;   % light gray background
ax.FontSize=14;

hold on;

% Plot left bars (negative values so they extend left)
barh(1:nn, -v1, 0.6, 'FaceColor',[238, 192, 163]/255); %

hold on;

% Plot right bars (positive values so they extend right)
barh(1:nn,  v2, 0.6, 'FaceColor',[152, 180, 206]/255);%[0.3 0.6 0.8]


% vertical separators
for i = 0.5 : 1 : nn+0.5
    yline(i, 'w-', 'LineWidth',1.6);
end

% Y-axis labels
set(gca,'YTick',1:nn,'YTickLabel',pairnames_ordered,'YDir','reverse','FontSize', 10.5);

% Zero line in the middle
xline(0,'k','LineWidth',1.5);

xlim([-8, 8]);
xlabel('% change','FontSize', 16,'FontWeight','bold');

ylabel('Pairwise enzyme change','FontSize', 16,'FontWeight','bold');
legend({'N-fixation efficiency (gC g^{-1}N)','N-fixation rate (mM s^{-1})'},'Location','best','FontSize', 13,'FontWeight','bold');

%%
rmpath('../kineticmodel/'); 


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

function C = mergeStructs(A, B)
    C = A;
    fields = fieldnames(B);
    for i = 1:numel(fields)
        C.(fields{i}) = B.(fields{i});
    end
end

