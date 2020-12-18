
%% plotting image with cells and neuropils
figure;
imagesc(brighten(brighten(dat.mimg(:,:,2),1),0.5));
hold on;

%spatialInfo.ROIs gives the ROIs only for the cells and not neuropils
for cell = 1:205
    plot(spatialInfo.ROIs{1,cell}(:,1), spatialInfo.ROIs{1,cell}(:,2), '.', 'Color', 'k')
end

%% events
event_length = events.eventsOff(2:end) - events.eventsOn;
duration_bet_events = events.eventsOn(2:end) - events.eventsOff(2:end-1);

% logical variable for cells
iscell = false(674,1);
for i = 1:674
    iscell(i) = dat.stat(i).iscell;
end

%% check if calcium npilTraces are selected from dat.FcellNeu and calcium rawTraces are selected from dat.Fcell

% found that dat.ca is the "smoothed" version of dat.Fcell

cnt = 1; 
check_neuropil = zeros(674,1);
check_cell = zeros(674,1);
for i = 1:length(iscell)
    if iscell(i)
        check_neuropil(i) = sum(calcium.npilTraces(cnt,:) == dat.FcellNeu{1,1}(i,:));
        check_cell(i) = sum(calcium.rawTraces(cnt,:) == dat.Fcell{1,1}(i,:));
        cnt = cnt + 1;
    end
end

%% check if spikes.raster are selected from dat.sp
cnt = 1;
check_spikes = zeros(674,1);
for i = 1:length(iscell)
    if iscell(i)
        check_spikes(i) = sum(spikes.raster(cnt,:) == dat.sp{1,1}(i,:));
        cnt = cnt + 1;
    end
end

%% run deconvolution on Ca traces

% found that sp is equal to spikes.raster in melanie's data
% found that ca is equal to dat.ca 
% found that coefs is equal to calcium.npilCoeffs

dat.ops.sensorTau = 1; % 0.7 for GCaMP6f, 1.0 for GCaMP6m, 1.25-1.5 for GCaMP6s
dat.ops.deconvType = 'OASIS';
dat.ops.recomputeKernel = true;
dat.ops.maxNeurop = 0.7;
dat.ops.estimateNeuropil = true;

[new_sp_runb, new_ca_runb, new_coefs_runb, new_B_runb, new_sd_runb, new_ops_runb, new_baselines_runb] = wrapperDECONV(dat.ops, calcium.rawTraces', calcium.npilTraces');



