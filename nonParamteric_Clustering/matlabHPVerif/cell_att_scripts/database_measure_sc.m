% Script performs three main tasks:
% (i) Reads data and metadat from an HDF5 database
% (ii) Analyses the imported data
% (iii) Saves data analysis output and metadata into ascii (plaintext) file
% appropriate for statistical analysis in R

clear

% Database info
dataFolder = '../../raw_data/';
databaseName = 'stellate_cell_recordings.h5';
databasePath = [dataFolder databaseName];
databaseInfo = h5info(databasePath);

% Output info
saveFolder = '../../raw_data/';
saveFileName = 'SC_HPdatatable.txt';

% Initialise storage structures for data and metadata
allData = []; % Measured properties
allMetadata.location = []; % Dorsal-ventral location
allMetadata.position = []; % Medial-lateral slice position
allMetadata.hemisphere = []; % Recorded hemisphere
allMetadata.age = []; % Animal age
allMetadata.housing = [];
allMetadata.totalCells = [];
allMetadata.id = []; % Animal ID
allMetadata.experimenter = []; % Experimenter ID
allMetadata.patchDirection = []; % General patch direction (not systematic)
allMetadata.recordingTimes = []; % Recording times (converted from Excel format)

% Each top level group in the database contains all data for a mouse
for nAnimal = 1:length(databaseInfo.Groups);
    
    animalName = databaseInfo.Groups(nAnimal).Name;
    
    % -------------------- Analyse subthreshold data ------------------------
    
    protocolName = 'subthresh';
    nSubthreshProperties = 4;
    nWaveforms = 5;
    
    % Indices of voltage and current command traces in each dataset
    indsVoltageTrace = 1:nWaveforms;
    indsCurrentTrace = indsVoltageTrace + nWaveforms;
    
    % List of recordings and associated information
    datasetsInfo = h5info(databasePath, [animalName '/' protocolName]);
    
    % Storage matrix for subthreshold stats
    subthreshData = zeros(length(datasetsInfo.Datasets), nSubthreshProperties);
    
    
    for nDataset = 1:length(datasetsInfo.Datasets)
        
        % Determine path to data
        datasetName = datasetsInfo.Datasets(nDataset).Name;
        datasetPath = [animalName '/' protocolName '/' datasetName];
        
        % Load the data and associated attribute(s)
        data = h5read(databasePath, datasetPath);
        sampleRate = h5readatt(databasePath, datasetPath, 'Sample rate');
        
        % Extract recorded voltage and current commands separately
        voltageData = double(data(:, indsVoltageTrace));
        currentData = double(data(:, indsCurrentTrace));
        
        % Measure subthreshold properties
        subthreshData(nDataset, :) = subthreshmeasure(sampleRate, voltageData, currentData);
        
    end
    
    % % -------------------- Analyse resonance data ------------------------
    
    protocolName = 'zap';
    nResonanceProperties = 2;
    nWaveforms = 1;
    
    % Indices of voltage and current command traces in each dataset
    indsVoltageTrace = 1:nWaveforms;
    indsCurrentTrace = indsVoltageTrace + nWaveforms;
    
    % List of recordings and associated information
    datasetsInfo = h5info(databasePath, [animalName '/' protocolName]);
    
    % Storage matrix for resonance stats
    zapData = zeros(length(datasetsInfo.Datasets), nResonanceProperties);
    
    
    for nDataset = 1:length(datasetsInfo.Datasets)
        
        % Determine path to data
        datasetName = datasetsInfo.Datasets(nDataset).Name;
        datasetPath = [animalName '/' protocolName '/' datasetName];
        
        % Load the data and associated attribute(s)
        data = h5read(databasePath, datasetPath);
        sampleRate = h5readatt(databasePath, datasetPath, 'Sample rate');
        
        % Extract recorded voltage and current commands separately
        voltageData = double(data(:, indsVoltageTrace));
        currentData = double(data(:, indsCurrentTrace));
        
        % Measure resonance properties
        zapData(nDataset, :)...
            = zapmeasure(sampleRate, voltageData, currentData);
        
    end
    
    % -------------------- Analyse spiking data ------------------------
    
    protocolName = 'ramp';
    nSpikeProperties = 5;
    nWaveforms = 1;
    
    % Indices of voltage and current command traces in each dataset
    indsVoltageTrace = 1:nWaveforms;
    indsCurrentTrace = indsVoltageTrace + nWaveforms;
    
    % List of recordings and associated information
    datasetsInfo = h5info(databasePath, [animalName '/' protocolName]);
    
    % Storage matrix for spiking stats
    spikeData = zeros(length(datasetsInfo.Datasets), nSpikeProperties);
    
    
    for nDataset = 1:length(datasetsInfo.Datasets)
        
        % Determine path to data
        datasetName = datasetsInfo.Datasets(nDataset).Name;
        datasetPath = [animalName '/' protocolName '/' datasetName];
        
        % Load the data and associated attribute(s)
        data = h5read(databasePath, datasetPath);
        sampleRate = h5readatt(databasePath, datasetPath, 'Sample rate');
        
        % Extract recorded voltage and current commands separately
        voltageData = double(data(:, indsVoltageTrace));
        currentData = double(data(:, indsCurrentTrace));
        
        % Measure spiking properties
        spikeData(nDataset, :)...
            = spikemeasure(sampleRate, voltageData, currentData);
        
    end
    
    % -------------------- Read all metadata --------------------
    
    % Anatomical information
    animalMetadata.location = zeros(length(datasetsInfo.Datasets), 1);
    animalMetadata.recordingTimes = zeros(length(datasetsInfo.Datasets), 1);
    animalMetadata.position = cell(length(datasetsInfo.Datasets), 1);
    animalMetadata.hemisphere = cell(length(datasetsInfo.Datasets), 1);
    % Animal information
    animalMetadata.age = cell(length(datasetsInfo.Datasets), 1);
    animalMetadata.id = cell(length(datasetsInfo.Datasets), 1);
    animalMetadata.housing = cell(length(datasetsInfo.Datasets), 1);
    % Experiment information
    animalMetadata.experimenter = cell(length(datasetsInfo.Datasets), 1);
    animalMetadata.patchDirection = cell(length(datasetsInfo.Datasets), 1);
    animalMetadata.totalCells = cell(length(datasetsInfo.Datasets), 1);
    
    for nDataset = 1:length(datasetsInfo.Datasets)
        
        % Path to dataset metadata
        datasetName = datasetsInfo.Datasets(nDataset).Name;
        datasetPath = [animalName '/' protocolName '/' datasetName];
        
        animalMetadata.location(nDataset) = ...
            h5readatt(databasePath, datasetPath, 'DV location');
        animalMetadata.position{nDataset} = ...
            h5readatt(databasePath, datasetPath, 'Slice position');
        animalMetadata.hemisphere{nDataset} = ...
            h5readatt(databasePath, datasetPath, 'Recorded hemisphere');
        animalMetadata.age{nDataset} = ...
            h5readatt(databasePath, datasetPath, 'Animal age');
        animalMetadata.id{nDataset} = animalName(2:end); % Remove filesep character
        animalMetadata.housing{nDataset} = ...
            h5readatt(databasePath, datasetPath, 'Housing');
        animalMetadata.experimenter{nDataset} = ...
            h5readatt(databasePath, datasetPath, 'Experimenter');
        animalMetadata.patchDirection{nDataset} = ...
            h5readatt(databasePath, datasetPath, 'Patch direction');
        animalMetadata.totalCells{nDataset} = ...
            h5readatt(databasePath, datasetPath, 'Total cells recorded');
        animalMetadata.recordingTimes(nDataset) = ...
            h5readatt(databasePath, datasetPath, 'Recording time');
        
        
    end
    
    % -------------------- Collate data ------------------------
    
    animalData = [subthreshData zapData spikeData];
    
    allData = [allData; animalData];
    allMetadata.location = [allMetadata.location; animalMetadata.location];
    allMetadata.position = [allMetadata.position; animalMetadata.position];
    allMetadata.hemisphere = [allMetadata.hemisphere; animalMetadata.hemisphere];
    allMetadata.age = [allMetadata.age; animalMetadata.age];
    allMetadata.id = [allMetadata.id; animalMetadata.id];
    allMetadata.housing = [allMetadata.housing; animalMetadata.housing];
    allMetadata.experimenter = [allMetadata.experimenter; animalMetadata.experimenter];
    allMetadata.patchDirection = [allMetadata.patchDirection; animalMetadata.patchDirection];
    allMetadata.totalCells = [allMetadata.totalCells; animalMetadata.totalCells];
    allMetadata.recordingTimes = [allMetadata.recordingTimes; animalMetadata.recordingTimes];
end

% -------------------- Format and save data ------------------------

% Configure all data as a single cell array
cellData = [num2cell(allData) num2cell(allMetadata.location)...
    allMetadata.position allMetadata.hemisphere...
    allMetadata.age allMetadata.id...
    allMetadata.housing allMetadata.experimenter...
     allMetadata.patchDirection allMetadata.totalCells...
     num2cell(allMetadata.recordingTimes)];

% Table headers
headers = {'vm' 'ir' 'sag' 'tau' 'resf' 'resmag' 'spkthr' 'spkmax' 'spkhlf'...
    'rheo' 'ahp' 'dvloc' 'mlpos' 'hemi' 'age' 'id' 'housing' 'expr' 'patchdir' 'totalcells' 'rectime'};
% 'ahp' 'fi'
% Convert cell array to a table
tableData = cell2table(cellData, 'VariableNames', headers);

% Write table to text file
currentFolder = cd(saveFolder);
writetable(tableData, saveFileName, 'Delimiter', 'tab')
cd(currentFolder)












