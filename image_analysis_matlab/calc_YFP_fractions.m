% This script calls the function time_YFP_hist to generate histogram of
% YFP intensities. It then plots fraction of cells that have reactivated
% YFP over time in the specified positions, both for each position
% separately and for all overlaid. It also plots cell count over time for
% each position. 

close all
clear all

dirname = '/data/abadiek/Imaging/Data_and_Processsing/20210628_CD8_D3_sort/20210628_CD8_D3_sort_raw_ANALYZED_20210810_21-38-02/';
% all postions used for Fig 4 summary line plot: 
positions = ["pos 28", "pos 30", "pos 32", "pos 34", "pos 35", "pos 36", "pos 37", "pos 39", "pos 40", "pos 41", "pos 42", "pos 44", "pos 45", "pos 47", "pos 48", "pos 49", "pos 50", "pos 51", "pos 53", "pos 54", "pos 56", "pos 57", "pos 59", "pos 64", "pos 65","pos 66", "pos 67", "pos 68", "pos 69", "pos 70", "pos 71", "pos 72", "pos 74", "pos 75","pos 76"]
% pos  used for Fig 4 snapshots:
% positions = ["pos 30", "pos 45", "pos 76"] 
filename = 'segprop.mat'
outfolder = '/data/abadiek/Imaging/Data_and_Processsing/20210628_CD8_D3_sort/20210628_CD8_D3_sort_raw_ANALYZED_20210810_21-38-02/YFP_frac/';
window = 5; % time window, hours
early_window = 25; % early time window, hours
n_std = 2; % standard deviations of the mean to draw YFP+/- gate
pos_thresh = 0.075;

% initialze vector to store starting cell numbers
clones_vec = NaN * ones(length(positions),1);
f_pos_final = NaN * ones(length(positions), 1);
t_react = NaN * ones(length(positions), 1);
tp_uni_pos = {}; % list
yy_pos = {}; % list
f_pos_vec_pos = {}; % list
tp_pos = {}; % list
yy_counts_pos={}; %list

% go through each position
for i = 1:length(positions)
    pos = positions(i)
    
    % run function to generate timepoint by YFP ave intensity vectors
    [tp, val] = time_YFP_hist(dirname, pos, filename);
    
    pos = convertStringsToChars(pos); % convert pos to char for figure naming 
    
    % remove nan (non fluor timepoints)
    tp_with_bf = tp; % store original timepoint vector just to count objects in each timepoint
    tp_with_bf = (tp_with_bf*5)/60; % convert to hours
    tp = tp(~isnan(val));
    val = val(~isnan(val));

    % remove negative values 
    tp = tp(val>0);
    val = val(val>0);

    % convert to hours
    tp = (tp*5)/60;
    
    % take log of yfp vals
    val_log = log10(val);
    
    % Get the distribution of YFP log vals from the first hours of the movie
    early_val = val_log(tp < early_window);
    % set yfp gate
    yfp_gate = mean(early_val) + n_std*std(early_val);
    
    % store starting cell number using full tp vector (including bf)
    early_tp = tp_with_bf(tp_with_bf < 5);
    counts = groupcounts(early_tp');
    clones = mode(counts);
    clones_vec(i) = clones;
    
    % get vector of unique timepoints
    tp_uni = unique(tp);
    
    % go through each timepoint and calculate YFP+ fraction 
    % first, intialize vectors to store f_pos
    f_pos_vec = NaN * ones(length(tp_uni),1);
    for j = 1:length(tp_uni)
        
        t = tp_uni(j);
        
        % get object yfp values corresponding to unique timepoint t
        n_obj = sum(tp == t);
        val_t = val_log(tp == t);
        n_pos = sum(val_t > yfp_gate);
        n_neg = sum(val_t <= yfp_gate);
        f_pos = n_pos / n_obj;
        f_pos_vec(j) = f_pos;

    end
    
    % store final fraction pos for each position
    f_pos_final(i) = mean(f_pos_vec(tp_uni > 70)); % take final as mean of all tp > 70 hr

    % smoothed frac over time
    xx = linspace(0,max(tp_uni),1000);
    yy = smooth(f_pos_vec, 20);
    
    % store each pos data
    tp_pos{i} = tp;
    tp_uni_pos{i} = tp_uni;
    f_pos_vec_pos{i} = f_pos_vec;
    yy_pos{i} = yy;
    
    % using smoothed data, quantify reactivation time 
    if f_pos_final(i) > pos_thresh
        tp_uni_test = tp_uni(tp_uni>40);
        yy_test = yy(tp_uni>40);   
        pass_thresh = tp_uni_test(yy_test > pos_thresh);
        t_react(i) = pass_thresh(1);
    else
        t_react(i) = NaN;
    end
    

   
    % for each position, plot f_pos over time
    figure(i)
    scatter(tp_uni, f_pos_vec, '.', 'k')
    hold on
    plot(tp_uni, yy, 'k', 'LineWidth', 2)
    hold off
    ylim([-.05,1])
    xlim([0,90])
    xticks([0, 20, 40, 60, 80])
    yticks([0, 0.2, 0.4, 0.6, 0.8])
    ax = gca; ax.FontSize = 32; ax.LineWidth=2;
    %set(gcf, 'Position',  [10, 10, 300, 600])

    % print density plot
    figname = strcat(outfolder, pos(~isspace(pos)), '_YFP_pos_frac','.pdf');
    print(figure(i), figname, '-dpdf')
    
    % for each position, plot object number over time
    figure(i+200)
    counts = groupcounts(tp')
    yy_counts = smooth(counts, 20);
    yy_counts_pos{i} = yy_counts;
    scatter(tp_uni,counts, '.', 'k'); hold on
    plot(tp_uni, yy_counts, 'k')
    hold off
    xticks([0, 20, 40, 60, 80])
    ax = gca; ax.FontSize = 32;
    
    % print  plot
    figname = strcat(outfolder, pos(~isspace(pos)), '_pos_counts','.pdf');
    print(figure(i+200), figname, '-dpdf')
    
    
          
end

% mean and std react time
print('mean and std deviation of reactivation time')
t_react_mean = mean(t_react(~isnan(t_react)))
t_react_std = std(t_react(~isnan(t_react)))
x1 = t_react_mean - t_react_std;
x2 = t_react_mean + t_react_std;

% make all position overlay plot
figure(600)
xline(t_react_mean, 'k', 'LineWidth',1, 'LineStyle','--')
hold on
area([x1 x2], [1 1], 'FaceColor','r', 'FaceAlpha', 0.2, 'LineStyle', 'none');
hold on 

for i = 1:length(positions)
    tp_uni = tp_uni_pos{i};
    yy = yy_pos{i};
    if f_pos_final(i) > pos_thresh
        plot(tp_uni,yy, 'r', 'LineWidth', 2); hold on
    else
        plot(tp_uni, yy, 'k', 'LineWidth', 2); hold on
    end
    
end
hold off 
ylim([-.05,1])
xlim([0,90])
xticks([0, 20, 40, 60, 80])
yticks([0, 0.2, 0.4, 0.6, 0.8])
ax = gca; ax.FontSize = 32; ax.LineWidth=2;

% print
figname = strcat(outfolder, 'YFP_pos_all','.pdf');
print(figure(600),figname, '-dpdf')

% make all position overlay plot for object number
figure(601)

for i = 1:length(positions)
    tp_uni = tp_uni_pos{i};
    yy_counts = yy_counts_pos{i};
    if f_pos_final(i) > pos_thresh
        plot(tp_uni,yy_counts, 'r'); hold on
    else
        plot(tp_uni, yy_counts, 'k'); hold on
    end
     
end
hold off 
xlim([0,90])
xticks([0, 20, 40, 60, 80])
ax = gca; ax.FontSize = 32;

% print
figname = strcat(outfolder, 'YFP_pos_all_counts','.pdf');
print(figure(601),figname, '-dpdf')


% for all positions, plot distribution of activation times

figure(602)
X = t_react';
histogram(X, 'BinWidth', 15)
ax = gca; ax.FontSize = 32;
% print
figname = strcat(outfolder, 'react_time','.pdf');
print(figure(602),figname, '-dpdf')

% store reactivation info
react = f_pos_final > pos_thresh;
pos_react = [positions', clones_vec, f_pos_final, react, t_react];
save([outfolder 'react_pos_info_matrix'], 'pos_react')
writematrix(pos_react,[outfolder 'react_pos_info_matrix.csv']) 

% react_pos
react_pos_only = positions(react);
save([outfolder 'react_pos_only'], 'react_pos_only'); 


