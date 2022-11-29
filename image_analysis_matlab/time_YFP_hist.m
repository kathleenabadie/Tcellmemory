function [tp, val] = time_YFP_hist(dirname, pos, filename)
% This function calculates histograms for the following:
% TCF1 average intensity 

cd(strcat(dirname,pos))
load(filename)


% intialize time and value vectors
tp = []; % store timepoint
val = [];
cell_counter = 1;
for t = 1:length(objects)
    cells = objects(t).obj;
    for c = 1:length(cells)
        value = objects(t).obj(c).data.TCF1_YFP_cor;
        tp(cell_counter) = t;
        val(cell_counter)=value;
        cell_counter = cell_counter + 1;
    end
end





        

        
    


