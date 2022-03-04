%%  Load data of interest - under /data/ folder

dir = '/home/janaki/Dropbox/project_with_melanie/DataForJanaki/Across_many_datasets/som/';
if ~exist(dir)
    mkdir(dir); 
    cd(dir); 
else
    cd(dir); 
end

load('DetectabilityVsDiscriminability_som.mat');

%% counting total cells in each category

numberOfDetectableHighCells = 0;
numberOfDetectableMidCells = 0;
numberOfDetectableNoCells = 0;
numberOfDiscriminableHighCells = 0;
numberOfDiscriminableMidCells = 0;
numberOfDiscriminableNoCells = 0;

for ii = 1:11
    numberOfDetectableHighCells = numberOfDetectableHighCells + length(DetectableCellsHighLight{ii});
    numberOfDetectableMidCells = numberOfDetectableMidCells + length(DetectableCellsMidLight{ii});
    numberOfDetectableNoCells = numberOfDetectableNoCells + length(DetectableCellsNoLight{ii});    
    numberOfDiscriminableHighCells = numberOfDiscriminableHighCells + length(DiscriminableCellsHighLight{ii});
    numberOfDiscriminableMidCells = numberOfDiscriminableMidCells + length(DiscriminableCellsMidLight{ii});
    numberOfDiscriminableNoCells = numberOfDiscriminableNoCells + length(DiscriminableCellsNoLight{ii});    
end