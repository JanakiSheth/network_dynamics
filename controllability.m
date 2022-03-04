%%  Load data of interest - under /data/ folder

clear variables

load('MT060_20190816_FRA_OptoStim_03.mat');
load('MT060_20190816_RedCells_SlopeRemoval_Movie.mat');
load('F_MT060_20190816MT060_tifStacks_plane1_proc.mat');

dir = '/home/janaki/Dropbox/project_with_melanie/DataForJanaki/MT060_20190816';
if ~exist(dir)
    mkdir(dir); 
    cd(dir); 
else
    cd(dir); 
end

%% 

for j = 1:length(dat.Fcell)
    Fcell{j}    = dat.Fcell{j};
    FcellNeu{j} = dat.FcellNeu{j};
end
    
Ff = [];
Fneu = [];
for j = 1:numel(Fcell)
    Ff   = cat(1, Ff, Fcell{j}');
    Fneu = cat(1,Fneu, FcellNeu{j}');
end

%% Ca and firing rates from suite2p with new parameter values

% found that sp is equal to spikes.raster in melanie's data
% found that calcium is equal to dat.ca 
% found that coefs is equal to calcium.npilCoeffs

dat.ops.sensorTau = 1; % 0.7 for GCaMP6f, 1.0 for GCaMP6m, 1.25-1.5 for GCaMP6s
dat.ops.deconvType = 'OASIS';
dat.ops.recomputeKernel = true;
dat.ops.maxNeurop = 0.7;
dat.ops.estimateNeuropil = false;
dat.ops.runningBaseline = true;

[new_sp_runb, new_ca_runb, new_coefs_runb, new_B_runb, new_sd_runb, new_ops_runb, new_baselines_runb] = wrapperDECONV(dat.ops, calcium.rawTraces', calcium.npilTraces');

% storing old variables since we will be changing them
calcium_raw_0 = calcium;
spikes_raw_0 = spikes;
calcium.npilCoeffs = new_coefs_runb';
calcium.npilSubTraces = calcium.rawTraces - calcium.npilCoeffs.*calcium.npilTraces;

spikes.raster = new_sp_runb';

%% calculate where laser bands are 
 % Janaki: when the laser shines, it makes bands of higher intensity in the
 % image, so this section is to find out where these bands are, which cells
 % are within the bands at that time frame (changes because of registration/movement), 
 % and replace the fluorescence at that time frame by a NaN. 

 % I will modify calcium directly, so make a copy in calcium_raw before
 % modifications
 
calcium_raw = calcium;

fr=exptInfo.fr; % frame rate

% find in original recording where the laser bands are

fs=stimInfo.fs;

index_laserevents=find(stimInfo.index(:,2)~=0);

% find laser band indices in terms of fs
indfs_laser=[];
for i=1:length(stimInfo.order)
    if ismember(stimInfo.order(i),index_laserevents)
        begin_indfs=exptInfo.preStimSilence*fs+(i-1)*fs*(stimInfo.ISI/1000+stimInfo.tDur_opto);
        % at the beginning of the recording, there is 20 seconds of silence before the stimulus file starts being presented. So the first stim arrives at frame ~600
        for j=1:stimInfo.nPulses % number of laser pulses  
            interval=stimInfo.evOn_opto(j):stimInfo.evOff_opto(j);
            indfs_laser=[indfs_laser begin_indfs+interval];
        end
        
    end    
end
   

% convert the indices in terms of frame lines where fr is camera frame rate
% and dat.ops.Ly is number of lines in a frame
indlinefr_laser=indfs_laser*(fr*dat.ops.Ly)/fs;  
indlinefr_laser=unique(round(indlinefr_laser));   

% go to the frame and lines where there is laser
indfr_laser=[];
indline_laser=[];
indline_reg_laser=[];
indcell_reg=zeros(size(calcium.npilSubTraces)); % index replaced by 1 if there is laser band

% find extremes in registered cell boundaries
for i=1:length(spatialInfo.ROIs)
    cell_reg_boundaries1(i,:)=[min(spatialInfo.ROIs{i}(:,1)) max(spatialInfo.ROIs{i}(:,1))];
end

for i=1:length(indlinefr_laser)
    indfr_laser(i)=ceil(indlinefr_laser(i)/dat.ops.Ly);
    indline_laser(i)=mod(indlinefr_laser(i),dat.ops.Ly);
    % line of registered frames
    indline_reg_laser(i)=indline_laser(i)+dat.ops.DS(indfr_laser(i),1);
    for j=1:length(spatialInfo.ROIs)
        if indline_reg_laser(i)>=cell_reg_boundaries1(j,1)&indline_reg_laser(i)<=cell_reg_boundaries1(j,2)
            calcium.npilSubTraces(j,indfr_laser(i)) = NaN;
            indcell_reg(j,indfr_laser(i))=1;
        end
    end
end



%% extract parameters

% Janaki: this uses the function makeCaRaster_NaN_MT.m

fr=exptInfo.fr; % frame rate
% these relevant_info strings are used later in titles for plotting - not
% relevant to making the raster
relevant_info=['t_{sound} = ',num2str(stimInfo.nClicks/stimInfo.clickRate),'s ; ','t_{opto} = ', num2str(stimInfo.tDur_opto),'s ; '];
relevant_info2 = ['max Amp_{sound} = ', num2str(max(stimInfo.intensity)),'dB; ','max Amp_{opto} = ', num2str(max(stimInfo.amplitude_opto*100)),'%'];

%  make the raster, eg cut up time series in each trial

preOnsetTime=1; %start raster 1 second before stim onset
postOnsetTime=round((events.eventsOn(2)-events.eventsOn(1))/fr)-1; % end raster x seconds after stim onset; here duration between two stimuli = 5 seconds; 
doDFoverF0=1; % =1 means it does the z-scoring -> (F-mean(Fbaseline))/std(Fbaseline)
raster = makeCaRaster_NaN_JS(calcium.npilSubTraces,events.eventsOn,round(preOnsetTime*fr),round(postOnsetTime*fr),doDFoverF0);

% average over all events, shape = time*no. of cells
av_raster = squeeze(nanmean(raster,1));

%% sort depending on components of index

% Janaki: uses function makeOrdRasterAndAveTrace_MT.m
% this rearranges the raster by stim combination and calculates average and
% std for each stim combination

repeats=stimInfo.repeats;

[raster_ordered,mean_Trace,std_mean_Trace]=makeOrdRasterAndAveTrace_MT(raster,stimInfo.index,stimInfo.order,stimInfo.repeats);

% number of components along each dimension of index
for i=1:size(stimInfo.index,2)
dim_index(i)=length(unique(stimInfo.index(:,i)));
end

dim_relevant=find(dim_index~=1);

% ordered in a  matrix, on several figures if necessary    

if length(dim_relevant)==2
    figure();
    hold on
    for ii=1:length(stimInfo.index)
        subplot(dim_index(dim_relevant(1)),dim_index(dim_relevant(2)),ii)
        imagesc([1:size(raster_ordered,2)]/fr,[1:size(raster_ordered,3)],squeeze(mean_Trace(ii,:,:))')
        colormap gray
        set(gca,'TickDir','out','FontSize',14)
        xlabel('Time (s)')
        ylabel('Cells')
        title({['index 1 = ',num2str(stimInfo.index(ii,dim_relevant(1)))],['index 2 = ',num2str(stimInfo.index(ii,dim_relevant(2)))]})
        if exist('RedCells')==1
            Ind_cells_red=find(RedCells.isCell_red);
            for i=1:length(Ind_cells_red)
                line([0 0.5],[Ind_cells_red(i) Ind_cells_red(i)],'Color','r','LineWidth',2)
            end
        end
    end
end

% figure with just the max parameter values in sound and laser

for i=1:length(dim_relevant)
    % find indices where dim_relevant is max
    b=find(stimInfo.index(:,dim_relevant(i))==max(stimInfo.index(:,dim_relevant(i))));
    % find indices where dim non relevant is zero
    if dim_relevant==1 % condition added for only laser changing and sound always on
        ind_max_dim_relevant(i)=b(find(stimInfo.index(b)));
    else
        ind_max_dim_relevant(i)=b(find(stimInfo.index(b,setdiff(dim_relevant,dim_relevant(i)))==0));
    end
end


if length(dim_relevant) > 1
    figure
    a=[1 fliplr(ind_max_dim_relevant) length(stimInfo.index)];
        hold on
    for ii=1:4
        subplot(2,2,ii)
        imagesc([1:size(raster_ordered,2)]/fr,[1:size(raster_ordered,3)],squeeze(mean_Trace(a(ii),:,:))')
        colormap gray
        set(gca,'TickDir','out','FontSize',14)
        xlabel('Time (s)')
        ylabel('Cells')
        title({['index 1 = ',num2str(stimInfo.index(a(ii),dim_relevant(1)))],['index 2 = ',num2str(stimInfo.index(a(ii),dim_relevant(2)))]})
        hold on
        if exist('RedCells')==1
            Ind_cells_red=find(RedCells.isCell_red);
            for i=1:length(Ind_cells_red)
                line([0 0.5],[Ind_cells_red(i) Ind_cells_red(i)],'Color','r','LineWidth',2)
            end
        end
    end
end

%close all;

%% sort for every condition depending on sound response and red intensity and show rasters

% Janaki: this is the section which looks at which cells are
% stim-responsive. I tried a few options out, the ones I'm not using
% anymore more are commented out, but I leave them there so you see what I
% tried. 


%%%% Option 2 : varying interval then take min of pstat t-test; adapated from
%%%% above; no correction to t-test because I then take only one window.

Tdelay_test=[0.1:0.03:1.5];
Tdelay_windur=1; % window duration of 300ms for averaging wih firing rates.
alpha=0.01;
 
Tdelay=[];
tstats_time=[];
pstats_time=[];
zscore_time=[];
raster_ordered_poststim = zeros(length(stimInfo.order),size(raster,3));
av_raster_poststim = zeros(length(stimInfo.index),size(raster,3));
sem_raster_poststim = zeros(length(stimInfo.index),size(raster,3));

for jj=1:size(raster,3)
    for ii=1:length(stimInfo.index)  
        % caulculate stats with all Tdelay_test
        for kk=1:length(Tdelay_test)
            [tstats_time(ii,jj,kk) pstats_time(ii,jj,kk)]=ttest(nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,1:round(fr*preOnsetTime),jj),2),nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,round(fr*(preOnsetTime+Tdelay_test(kk))):round(fr*(preOnsetTime+Tdelay_test(kk)+Tdelay_windur)),jj),2),'Alpha',alpha);
            zscore(ii,jj,kk)=(mean(nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,1:round(fr*preOnsetTime),jj),2))-mean(nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,round(fr*(preOnsetTime+Tdelay_test(kk))):round(fr*(preOnsetTime+Tdelay_test(kk)+Tdelay_windur)),jj),2)))/sqrt(0.5*(std(nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,1:round(fr*preOnsetTime),jj),2))^2+std(nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,round(fr*(preOnsetTime+Tdelay_test(kk))):round(fr*(preOnsetTime+Tdelay_test(kk)+Tdelay_windur)),jj),2))^2)^2);
        end
        Tdelay(ii,jj)=min(Tdelay_test(find(pstats_time(ii,jj,:)==min(pstats_time(ii,jj,:)))));

        % compute stats, average and std with chosen Tdelay
        [tstats(ii,jj) pstats(ii,jj)]=ttest(nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,1:round(fr*preOnsetTime),jj),2),nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,round(fr*(preOnsetTime+Tdelay(ii,jj))):round(fr*(preOnsetTime+Tdelay(ii,jj)+Tdelay_windur)),jj),2),'Alpha',alpha);
        b2=squeeze(nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,round(fr*(preOnsetTime+Tdelay(ii,jj))):round(fr*(preOnsetTime+Tdelay(ii,jj)+Tdelay_windur)),jj),2));
        raster_ordered_poststim(repeats*(ii-1)+1:repeats*ii,jj) = (b2-nanmean(b2))/abs(nanstd(b2));
        av_raster_poststim(ii,jj)=nanmean(b2);
        sem_raster_poststim(ii,jj)=nanstd(b2)/sqrt(repeats);
    end
end

for i=1:length(dim_relevant)
    [a Ind_sort_poststim{i}]=sort(av_raster_poststim(ind_max_dim_relevant(i),:),'descend');
    mean_Trace_ord{i}=mean_Trace(:,:,Ind_sort_poststim{i});
end


 %%%%% mark which cells show an activity change, not ordered
Cells_ActivityChange=cell(length(stimInfo.index),1);
Cells_ActivityChange2=zeros(length(stimInfo.index),size(raster,3));
for kk=1:length(stimInfo.index)
    Cells_ActivityChange{kk}=zeros(1,size(raster,3));
    for jj=1:size(raster,3)
        if pstats(kk,jj)<alpha && av_raster_poststim(kk,jj) >0
            Cells_ActivityChange{kk}(jj)=1;
            Cells_ActivityChange2(kk,jj)=1;
        elseif pstats(kk,jj)<alpha && av_raster_poststim(kk,jj) <0
            Cells_ActivityChange{kk}(jj)=-1;
            Cells_ActivityChange2(kk,jj)=-1;
        end
    end
end

%% do we have amplitude tuning?

for ii = 1:size(raster,3)
    figure(); hold on;
    %errorbar(av_raster_poststim([1,4,7,10,13,16,19], ii), sem_raster_poststim([1,4,7,10,13,16,19], ii));    
    errorbar(av_raster_poststim([2,5,8,11,14,17,20], ii), sem_raster_poststim([2,5,8,11,14,17,20], ii)); 
    %errorbar(av_raster_poststim([3,6,9,12,15,18,21], ii), sem_raster_poststim([3,6,9,12,15,18,21], ii)); 
end

cellsIncreasing = [176,172,169,155,148,131,111,109,104,103,97,95,93,92,91,82,80,78,74,67,64,61,60,57,55,48,40,39,33,25,23,19,17,15,8,7];
cellsDecreasing = [165,160,146,139,137,135,123,114,116,101,100,73,59,52,49,18,12,1];
cellsTunedUp = [154,145,142,134,130,125,124,107,99,89,88,85,84,81,79,75,69,68,62,41,38,31,30,27,26,21,16,14,6,5];
cellsTunedDown = [166,163,158,153,147,144,140,132,126,112,108,90,87,72,50,47,34,28,11,4];

for ii = 1:length(cellsIncreasing)
    h = patch(spatialInfo.ROIs{cellsIncreasing(ii)}(:,2),spatialInfo.ROIs{cellsIncreasing(ii)}(:,1),'y','facealpha',0.8,'edgealpha',0,...
        'LineWidth',1,'facecolor','r');
end

for ii = 1:length(cellsDecreasing)
    h = patch(spatialInfo.ROIs{cellsDecreasing(ii)}(:,2),spatialInfo.ROIs{cellsDecreasing(ii)}(:,1),'y','facealpha',0.8,'edgealpha',0,...
        'LineWidth',1,'facecolor','y');
end

for ii = 1:length(cellsTunedUp)
    h = patch(spatialInfo.ROIs{cellsTunedUp(ii)}(:,2),spatialInfo.ROIs{cellsTunedUp(ii)}(:,1),'y','facealpha',0.8,'edgealpha',0,...
        'LineWidth',1,'facecolor','k');
end

for ii = 1:length(cellsTunedDown)
    h = patch(spatialInfo.ROIs{cellsTunedDown(ii)}(:,2),spatialInfo.ROIs{cellsTunedDown(ii)}(:,1),'y','facealpha',0.8,'edgealpha',0,...
        'LineWidth',1,'facecolor','magenta');
end

%% extended differential covariance 

% build a matrix with rows time points and stim conds and columns cells
% cov(dX/dt, X), cov(X,X) gives A.

X = zeros(prod(size(mean_Trace,[1,2])),size(mean_Trace,3));
X = zeros(prod(size(mean_Trace,[1,2])),size(mean_Trace,3));

for ii=1:size(mean_Trace,1)
    X(((ii-1)+1:ii)*size(mean_Trace,2),:) = mean_Trace(ii,:,:);
    diffX
end

%% euclidean distance between states (an idea from the grant)
close all;

dist = zeros(21, size(raster,2));
for cc = 2:21
    for ii = 1:size(raster,2)
        dist(cc,ii) = vecnorm(mean_Trace(cc,ii,:)-mean_Trace(1,ii,:),2,3);
    end
end

figure(); hold on;
subplot(3,3,1); hold on;
plot(dist(2,:));
plot(dist(3,:));

subplot(3,3,2); hold on;
plot(dist(4,:))
plot(dist(5,:))
plot(dist(6,:))

subplot(3,3,3); hold on;
plot(dist(7,:))
plot(dist(8,:))
plot(dist(9,:))

subplot(3,3,4); hold on;
plot(dist(10,:))
plot(dist(11,:))
plot(dist(12,:))

subplot(3,3,5); hold on;
plot(dist(13,:))
plot(dist(14,:))
plot(dist(15,:))

subplot(3,3,6); hold on;
plot(dist(16,:))
plot(dist(17,:))
plot(dist(18,:))

subplot(3,3,7); hold on;
plot(dist(19,:))
plot(dist(20,:))
plot(dist(21,:))

figure(); hold on;
for ii = 1:3
    % subplot(3,1,ii); hold on;
    if ii==1
        c = 'r';
    elseif ii==2
        c = 'b';
    else
        c = 'k';
    end
    plot(dist(ii,:), 'Color', c);
    plot(dist(ii+3,:), 'Color', c);
    plot(dist(ii+6,:), 'Color', c);
    plot(dist(ii+9,:), 'Color', c);
    plot(dist(ii+12,:), 'Color', c);
    plot(dist(ii+15,:), 'Color', c);
    plot(dist(ii+18,:), 'Color', c);
end
