function [cell_atts] = calculate_cell_atts(cell_type)
% [cell_atts] = calculate_cell_atts(cell_type)
%   Cell_type refers to data set. Takes HP raw current clamp data in hd5
%   structure and returns an array with n rows and 12 columns
%   representing 12 extracted variables:
%    attVarNames = {'spkthr',...
%             'spkmax',...
%             'spkhlf',...
%             'rheo',...
%             'ahp',...
%             'vm',...
%             'ir',...
%             'sag',...
%             'tau',...
%             'resf',...
%             'resmag',...
%             'dvloc'};
%
%   Last Edited by Oliver Shipston-Sharman 10/06/18

switch nargin
    case 0
        cell_type = 0;
end

switch cell_type
    case 0 % Stellate cells only 
        filename = {'../raw_data/stellate_cell_recordings.h5'};
        dv_loc_idx = 11;
    case 2 % Both
        filename = {'path_to_raw_data/calbindin_cell_recordings.h5'};
        dv_loc_idx = 7;
    case 2 % Both
        filename = [{'path_to_raw_data/stellate_cell_recordings.h5'},...
                    {'path_to_raw_data/calbindin_cell_recordings.h5'}];
        dv_loc_idx = [11 7];
end
         
for file = 1:length(filename)
    % Collect hd5 file structure...
    info = h5info(filename{file});
    g1 = info.Groups;
    for m = 1:length(g1) % For each mouse...
        g2 = info.Groups(m).Groups;
        for I_input = 1:length(g2) % For each category of experiment...
            path = g2(I_input).Name;
            clear v_dat I_dat
            for c = 1:length(g2(I_input).Datasets)  % For each cell...
                datasetname = info.Groups(m).Groups(I_input).Datasets(c).Name;
                fs = g2(I_input).Datasets(c).Attributes(2).Value;
                data = h5read(filename{file},[path '/' datasetname]);
                if I_input == 1 % Ramp
                    v_dat = {data(:,1)};
                    I_dat = {data(:,2)};
                    cell_dat(c,1:5) = spikemeasure(fs, v_dat, I_dat);
                    cell_label{1,c} = g2(I_input).Datasets(c).Name(1:end-4); % Cell Label
                    cell_label{2,c} = g2(I_input).Datasets(c).Attributes(3).Value; % Experimenter
                    cell_label{3,c} = g2(I_input).Datasets(c).Attributes(5).Value; % Animal Age
                    cell_label{4,c} = g2(I_input).Datasets(c).Attributes(6).Value; % Hemisphere
                    cell_label{5,c} = g2(I_input).Datasets(c).Attributes(7).Value; % Total Cells
                    cell_label{6,c} = g2(I_input).Datasets(c).Attributes(8).Value; % Housing
                    cell_label{7,c} = g2(I_input).Datasets(c).Attributes(9).Value; % Patch Direction
                    cell_label{8,c} = g2(I_input).Datasets(c).Attributes(10).Value; % Slice Position
                    cell_label{9,c} = g2(I_input).Datasets(c).Attributes(12).Value; % Recording Time
                elseif I_input == 2 % Subthresh
                    for ep = 1:5
                        v_dat{ep} = data(:,ep);
                        I_dat{ep} = data(:,5+ep);
                    end
                    cell_dat(c,6:9) = subthreshmeasure(fs, v_dat, I_dat);
                elseif I_input == 3 % zap
                    v_dat = {data(:,1)};
                    I_dat = {data(:,2)};
                    cell_dat(c,10:11) = zapmeasure(fs, v_dat, I_dat);
                end
                cell_dat(c,12) = g2(I_input).Datasets(c).Attributes(dv_loc_idx(file)).Value;
            end
        end
        cell_atts(m).mouse = cell_dat; % Store Cell att data
        cell_atts(m).labels = {cell_label};
        disp(['Mouse:' num2str(m) '/27']);
    end
end
end

