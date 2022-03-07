function [raster_ordered, mean_Trace, std_mean_Trace]=makeOrdRasterAndAveTrace_MT(raster,stimindex,order,repeats)

% Takes a raster, orders it by stimuli conditions  and outputs the ordered
% raster, and the average and std of each stimulus condition

% stimindex=stimInfo.index;
% repeats=stimInfo.repeats;
% order=stimInfo.order;

% Inputs : 

% raster : a matrix of dimension trials x time x ncells ; the ordering will
%          be along the first dimension
% stimindex : stimulus conditions
% order : order in which the stimuli conditions are presented
% repeats : number of repeats for each stim condition

% Outputs

% raster_ordered : raster with trials organized by stim condition
% mean_Trace : nanmean trace for each stim condition across the repeats
% std_mean_Trace : nanstd trace for each stim condition across the repeats


for ii=1:length(stimindex)
    TrialOrder(ii,:) = find(order==ii);
    raster_ordered(repeats*(ii-1)+1:repeats*ii,:,:)=raster(TrialOrder(ii,:),:,:);
end

% average the trials for each test frequency

mean_Trace=zeros(length(stimindex),size(raster,2),size(raster,3));

for ii=1:length(stimindex)  
    mean_Trace(ii,:,:) = nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,:,:),1);
    std_mean_Trace(ii,:,:) = nanstd(raster_ordered(repeats*(ii-1)+1:repeats*ii,:,:),1);
end

end