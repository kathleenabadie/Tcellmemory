% This script calls the function time_YFP_hist to generate histogram of YFP
% intensities for a specified time window and plots


dirname = '/data/abadiek/Imaging/Data_and_Processsing/20210628_CD8_D3_sort/20210628_CD8_D3_sort_raw_ANALYZED_20210810_21-38-02/';
positions = ["pos 30", "pos 45", "pos 76"] % pos  used for Fig 4 

time_slices = [10, 30, 50, 70];
filename = 'segprop.mat'
outfolder = '/data/abadiek/Imaging/Data_and_Processsing/20210628_CD8_D3_sort/20210628_CD8_D3_sort_raw_ANALYZED_20210810_21-38-02/YFP_hist_1D/';
bw = 0.1;
window = 5; % time window, hours
early_window = 25; % early time window, hours
n_std = 2; % standard deviations of the mean to draw YFP+/- gate

for i = 1:length(positions)
    pos = positions(i)
    
    % run function to generate timepoint by YFP ave intensity vectors
    [tp, val] = time_YFP_hist(dirname, pos, filename);
    
    % remove nan (non fluor timepoints)
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
  
    for j = 1:length(time_slices)
        
        pos = convertStringsToChars(pos);
        
        slice = time_slices(j)
        
%         % histogram of only a single time point
%         [min_diff, min_loc] =  min(abs(tp - time_slice));
%         slice_real = tp(min_loc);
%         slice_val = val_log(tp == slice_real);
        
        % alternatively, take multiple timepoints within a window
        locs = tp > (slice - window) & tp < (slice + window);
        slice_val = val_log(locs);
        X = [slice_val'];
      
        % density plot
        h = histogram(X, 'BinWidth', bw, 'normalization', 'pdf');
        cnts = h.BinEdges + h.BinWidth / 2;
        cnts = cnts(1:end-1);
        h_val = h.Values;
        if ~isnan(h_val)
            [f,xi] = ksdensity(X);

            figure(2)
            scatter(h_val, cnts, 38, 'k', 'filled')
            %scatter(h_val, cnts, '.', 'k')
            hold on
            plot(f,xi,'k', 'LineWidth', 2)
            %plot(f,xi,'k')
            hold on
            yfp_gate_vec = yfp_gate*ones(100,1);
            y_vec = linspace(0,max(max(f), max(h_val)), 100);
            plot(y_vec, yfp_gate_vec, 'r', 'LineWidth', 2)
            hold off
            ylim([0,3.5])
            set(gcf, 'Position',  [10, 10, 300, 600])
            ax = gca; ax.FontSize = 36; ax.LineWidth=2;
            if slice ~= 10
                ax.YTickLabel=[];
            end

            % print density plot
            figname = strcat(outfolder, pos(~isspace(pos)), '_t', string(slice), '_density_1D','.pdf');
            figure(2)
            print(figname, '-dpdf')
        end
              
    end
end