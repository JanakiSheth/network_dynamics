%% Noise click + optogenetics Analysis : 

% OUTLINE
% - load the data
% - get fr from spikes
% - remove laser bands from calcium trace
% - make raster triggered on events
% - do stats on each condition to see what cells respond, and from that
% get which cells are sound-responsive (from max sound no opto condition)
% or laser-responsive
% - plot curve for sound-responsive cells of response versus sound
% amplitude for different laser levels (look at indiv cells and average)
% - use 99% criterion and 1-s moving window between [0 1] or [1 2] seconds after sound
% onset 


%%  Load data of interest - under /data/ folder

clear variables
close all

load('ImagingDatasetsAndCorrespondingInfo_DirectlyFromExpts/MT063_20190822_FRA_OptoStim_03.mat');
load('ImagingDatasetsAndCorrespondingInfo_DirectlyFromExpts/MT063_20190822_RedCells_SlopeRemoval_Movie.mat');
load('ImagingDatasetsAndCorrespondingInfo_DirectlyFromExpts/F_MT063_20190822MT063_tifStacks_plane1_proc.mat');

dir = '/home/janaki/Dropbox/project_with_melanie/DataForJanaki/Across_many_datasets/som';
if ~exist(dir)
    mkdir(dir); 
    cd(dir); 
else
    cd(dir); 
end

changeNumberOfRepeats = 0;
if length(calcium.rawTraces)<39100
    oldCalciumRawTraces = calcium.rawTraces;
    calcium.rawTraces = zeros(size(calcium.rawTraces,1),39100);
    calcium.rawTraces(:,1:38100) = oldCalciumRawTraces;
    calcium.rawTraces(:,38101:39100) = nan;
    oldCalciumNPilTraces = calcium.npilTraces;
    calcium.npilTraces = zeros(size(calcium.npilTraces,1),39100);
    calcium.npilTraces(:,1:38100) = oldCalciumNPilTraces;
    calcium.npilTraces(:,38101:39100) = nan;
    oldCalciumNPilSubtraces = calcium.npilSubTraces;
    calcium.npilSubTraces = zeros(size(calcium.npilSubTraces,1),39100);
    calcium.npilSubTraces(:,1:38100) = oldCalciumNPilSubtraces;
    calcium.npilSubTraces(:,38101:39100) = nan;
    changeNumberOfRepeats = 1;
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

%% Get FR from spikes

spikes_var = spikes.raster;
firing_rate = zeros(size(spikes_var));

sigma = 4; % sigma = 3 is chosen since the uncertainty in firing rate is ~100 ms (Paper:Computational processing of calcium imaging data)
x = [-12:1:12]; % 3*sigma = 9, hence sampling the gaussian dist over a longer range than 9
k = exp(-(x/sigma).^2/2)/(sigma*sqrt(2*pi)); % Gaussian kernel
for i_spike = 1:size(spikes_var,1)
    temp = conv(spikes_var(i_spike,:),k);
    firing_rate(i_spike,:) = temp(ceil(length(k)/2):end-floor(length(k)/2));
end
f = figure();
plot(spikes_var(1,:));
hold on
h=plot(firing_rate(1,:), '-r');
h.Color(4) = 0.5;  % 50% transparent

f = figure();
plot(spikes_var(10,:));
hold on
h=plot(firing_rate(10,:), '-r');
h.Color(4) = 0.5;  % 50% transparent

f = figure();
plot(spikes_var(20,:));
hold on
h=plot(firing_rate(20,:), '-r');
h.Color(4) = 0.5;  % 50% transparent

legend('Spike Train', 'Rate');

%% calculate when laser bands are 
 % Janaki: when the laser shines, it makes bands of higher intensity in the
 % image, so this section is to find out where these bands are, which cells
 % are within the bands at that time frame (changes because of registration/movement), 
 % and replace the fluorescence at that time frame by a NaN. 

 % I will modify calcium directly, so make a copy in calcium_raw before
 % modifications
 
FR_raw=firing_rate;
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
            firing_rate(j,indfr_laser(i)) = NaN;
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
% raster = makeCaRaster_NaN_JS(firing_rate,events.eventsOn,round(preOnsetTime*fr),round(postOnsetTime*fr),doDFoverF0);
raster = makeCaRaster_NaN_JS(calcium.npilSubTraces,events.eventsOn,round(preOnsetTime*fr),round(postOnsetTime*fr),doDFoverF0);

% average over all events, shape = time*no. of cells
av_raster = squeeze(nanmean(raster,1));
% av_raster_ca = squeeze(nanmean(raster_ca,1));

%% sort depending on components of index

% Janaki: uses function makeOrdRasterAndAveTrace_MT.m
% this rearranges the raster by stim combination and calculates average and
% std for each stim combination

repeats=stimInfo.repeats;

[raster_ordered,mean_Trace,std_mean_Trace]=makeOrdRasterAndAveTrace_MT(raster,stimInfo.index,stimInfo.order,stimInfo.repeats);
% [raster_ca_ordered,mean_ca_Trace,std_mean_ca_Trace]=makeOrdRasterAndAveTrace_MT(raster_ca,stimInfo.index,stimInfo.order,stimInfo.repeats);

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


%% sort for every condition depending on sound response and red intensity and show rasters

% Janaki: this is the section which looks at which cells are
% stim-responsive. I tried a few options out, the ones I'm not using
% anymore more are commented out, but I leave them there so you see what I
% tried. 


% Option 2 : varying interval then take min of pstat t-test; adapated from
% above; no correction to t-test because I then take only one window.

Tdelay_windur=1; % window duration of 1s for averaging.
Tdelay_test=[0:0.03:postOnsetTime-Tdelay_windur];
alpha=0.01;
if changeNumberOfRepeats
    repeats=9;
end

Tdelay=[];
tstats_time=[];
pstats_time=[];
zscore_time=[];

% calculate stats for all cells, all stims, all windows
for jj=1:size(raster,3)
for ii=1:length(stimInfo.index)  
    % calculate stats with all Tdelay_test
    for kk=1:length(Tdelay_test)
     [tstats_time(ii,jj,kk) pstats_time(ii,jj,kk)]=ttest(nanmean(raster_ordered((ii-1)*10+1:10*(ii-1)+repeats,1:round(fr*preOnsetTime),jj),2),...
                                                        nanmean(raster_ordered((ii-1)*10+1:10*(ii-1)+repeats,round(fr*(preOnsetTime+Tdelay_test(kk))):round(fr*(preOnsetTime+Tdelay_test(kk)+Tdelay_windur)),jj),2),'Alpha',alpha);
     zscore(ii,jj,kk)=(mean(nanmean(raster_ordered((ii-1)*10+1:10*(ii-1)+repeats,1:round(fr*preOnsetTime),jj),2))...
                        -mean(nanmean(raster_ordered((ii-1)*10+1:10*(ii-1)+repeats,round(fr*(preOnsetTime+Tdelay_test(kk))):round(fr*(preOnsetTime+Tdelay_test(kk)+Tdelay_windur)),jj),2)))...
                        /sqrt(0.5*(std(nanmean(raster_ordered((ii-1)*10+1:10*(ii-1)+repeats,1:round(fr*preOnsetTime),jj),2))^2+...
                        std(nanmean(raster_ordered((ii-1)*10+1:10*(ii-1)+repeats,round(fr*(preOnsetTime+Tdelay_test(kk))):round(fr*(preOnsetTime+Tdelay_test(kk)+Tdelay_windur)),jj),2))^2)^2);
    % average pop response for each time window  
    b2=squeeze(nanmean(raster_ordered((ii-1)*10+1:10*(ii-1)+repeats,round(fr*(preOnsetTime+Tdelay_test(kk))):round(fr*(preOnsetTime+Tdelay_test(kk)+Tdelay_windur)),jj),2));
    av_raster_poststim_time(ii,jj,kk)=nanmean(b2);
    sem_raster_poststim_time(ii,jj,kk)=nanstd(b2)/sqrt(repeats);
    end
end
end

% Based on conversation with Melanie, we either consider having the same
% time delay for all the cells for a given stimulus condition determined 
% using population analysis or we have varying time delays for each cell
% per stimulus condition as captured in the variable Tdelay_test.

for kk=1:length(Tdelay_test)
    Nsignifcells_allstims_time(kk)=length(find(sum(tstats_time(:,:,kk))>=1));
end
Tdelay_allstims=Tdelay_test(find(Nsignifcells_allstims_time(:)==max(Nsignifcells_allstims_time(:)),1,'first'));

cm=[0 0 0;
    1 0.7 0.7;
    1 0 0];

set(groot,'defaultAxesColorOrder',cm)

figure, 
hold on, 
plot(Tdelay_test,Nsignifcells_allstims_time,'b','LineWidth',2)
plot(Tdelay_allstims,max(Nsignifcells_allstims_time),'bo','LineWidth',2,'MarkerSize',12)
xlabel('Onset of response window (s)')
ylabel('Number of significantly responsive cells')
title({['mouse, - ',num2str(date),' - ',num2str(size(raster,3)),' cells'],...
    ['Number of signif cells vs onset of time window'],...
     ['to at least one stim condition stim conditions (blue)'],...
    ['Window begin = ',num2str(Tdelay_allstims),'s, Number of signif cells = ',num2str(max(Nsignifcells_allstims_time))],...
    ['Duration window = ',num2str(Tdelay_windur),'s, alpha = ',num2str(alpha)]})

for ii=1:length(stimInfo.index)
    % compute stats, average and std with chosen Tdelay
    tstats(ii,:)=tstats_time(ii,:,find(Nsignifcells_allstims_time(:)==max(Nsignifcells_allstims_time(:)),1,'first'));
    pstats(ii,:)=pstats_time(ii,:,find(Nsignifcells_allstims_time(:)==max(Nsignifcells_allstims_time(:)),1,'first'));
    av_raster_poststim(ii,:)=av_raster_poststim_time(ii,:,find(Nsignifcells_allstims_time(:)==max(Nsignifcells_allstims_time(:)),1,'first'));
    sem_raster_poststim(ii,:)=sem_raster_poststim_time(ii,:,find(Nsignifcells_allstims_time(:)==max(Nsignifcells_allstims_time(:)),1,'first'));
end

for i=1:length(dim_relevant)
    [a Ind_sort_poststim{i}]=sort(av_raster_poststim(ind_max_dim_relevant(i),:),'descend');
    mean_Trace_ord{i}=mean_Trace(:,:,Ind_sort_poststim{i});
end


 % mark which cells show an activity change, not ordered
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

%% Testing the idea of vips being intensity tuned (from Functional response properties of VIP-expressing inhibitory neurons in mouse visual and auditory cortex)

for iRed = 1:length(Ind_cells_red)
    figure();
    plot(av_raster_poststim([1,4,7,10,13,16,19],Ind_cells_red(iRed)), 'blue');
    hold on; 
    plot(av_raster_poststim([2,5,8,11,14,17,20],Ind_cells_red(iRed)), 'red');
    plot(av_raster_poststim([3,6,9,12,15,18,21],Ind_cells_red(iRed)), 'black');
end

%% Sort cells depending on different criteria

% modified 2019/09/26 : make cells be only Sound+ or only Sound-
criterion=cell(9,1);
criterion_short=cell(7,1);
Ind_CellSelection=cell(9,1);

%%%% 1) All cells %%%%
    Ind_CellSelection{1}=[1:size(av_raster,2)];
    criterion{1} = 'All cells';
    criterion_short{1} = 'All';
    
% %%%% 2) All red cells %%%%
if exist('RedCells')==1
    Ind_CellSelection{2}=Ind_cells_red;
    criterion{2} = 'All red cells';    
    criterion_short{2} = 'Red';
end

%%%% 3) All sound-responsive cells, increasing %%%%
    b=find(stimInfo.index(:,2)==0 & stimInfo.index(:,1)~=0);
    Ind_CellSelection{3}=[];
    Ind_CellSelection{4}=[];
    for jj=1:size(av_raster,2)
        if sum(abs(Cells_ActivityChange2(b,jj)))>0
            b1=find(Cells_ActivityChange2(b,jj)~=0);
            c=find(abs(av_raster_poststim(:,jj))==max(abs(av_raster_poststim(b(b1),jj))));
            if Cells_ActivityChange2(c,jj)==1
                Ind_CellSelection{3}=union(Ind_CellSelection{3},jj);
            elseif Cells_ActivityChange2(c,jj)==-1
                Ind_CellSelection{4}=union(Ind_CellSelection{4},jj);
            end
        end
    end
    Ind_CellSelection{3}=reshape(Ind_CellSelection{3},[],1);
    Ind_CellSelection{4}=reshape(Ind_CellSelection{4},[],1);
    
    criterion{3} = 'All sound-responsive cells, increasing';
    criterion_short{3} = 'Sound+';
    criterion{4} = 'All sound-responsive cells, decreasing'; 
    criterion_short{4} = 'Sound-';
    
%%%% 5) All laser-responsive cells, increasing %%%%
    b=find(stimInfo.index(:,1)==0 & stimInfo.index(:,2)~=0);
    Ind_CellSelection{5}=[];

    for jj=1:size(av_raster,2)
        if sum(abs(Cells_ActivityChange2(b,jj)))>0
            b1=find(Cells_ActivityChange2(b,jj)~=0);
            c=find(abs(av_raster_poststim(:,jj))==max(abs(av_raster_poststim(b(b1),jj))));
            if Cells_ActivityChange2(c,jj)==1
            Ind_CellSelection{5}=union(Ind_CellSelection{5},jj);
            elseif Cells_ActivityChange2(c,jj)==-1
            Ind_CellSelection{6}=union(Ind_CellSelection{6},jj);
            end
        end
    end
    Ind_CellSelection{5}=reshape(Ind_CellSelection{5},[],1);
    Ind_CellSelection{6}=reshape(Ind_CellSelection{6},[],1);    
    
    criterion{5} = 'All laser-responsive cells, increasing';
    criterion_short{5} = 'Laser+';
%%%% 6) All laser-responsive cells, decreasing %%%%
    criterion{6} = 'All laser-responsive cells, decreasing';
    criterion_short{6} = 'Laser-';
    
%%%% 7) Cells that don't respond and are not red %%%%
    Ind_CellSelection{7}=[1:size(av_raster,2)];
    for q=2:6
        Acommon = intersect(Ind_CellSelection{7},Ind_CellSelection{q});
        Ind_CellSelection{7} = setxor(Ind_CellSelection{7},Acommon);
      % Ind_CellSelection{7}(ismember(Ind_CellSelection{q},Ind_CellSelection{7}))=[]; 
    end
    criterion{7} = 'Cells not repsonsive to laser or sound and not red';
    criterion_short{7} = 'none';
    
%%%% 8) Cells responsive to max sound amplitude with laser off %%%%
    Ind_CellSelection{8}=find(Cells_ActivityChange{ind_max_dim_relevant(1)}==1); %dim index 1 = sound
    criterion{8} = 'Cells responsive to max sound amplitude';
    
%%%% 9) Particular cell %%%%
    Ind_CellSelection{9}=61;%
    criterion{9} = 'Cell 61';
    
    disp(size(Ind_CellSelection{5})/size(Ind_CellSelection{6}));
%% Plot curve cell response vs sound amplitude for stim- responsive-cells
% figure with overlap individual sound-responsive cell curves, different
% laser for different panels
% figure with average for the different laser amplitudes overlapped

for q=[1:6]
    fig(q+1)=figure(q+10);
    if ~isempty(Ind_CellSelection{q})
        % plot the response curves for all individual cells
        for ii=1:dim_index(2) % dim index 2 = laser
            a=[ii:dim_index(2):length(stimInfo.index)];
            subplot(1,dim_index(2)+1,ii)
            hold  on
            if length(Ind_CellSelection{q})==1
                errorbar(repmat(stimInfo.index(a,1),1,length(Ind_CellSelection{q})),av_raster_poststim(a,Ind_CellSelection{q}),sem_raster_poststim(a,Ind_CellSelection{q}),'+-')
            else
                plot(stimInfo.index(a,1),av_raster_poststim(a,Ind_CellSelection{q}),'+-')   
            end
            xlabel('Sound amplitude (dB)')
            ylabel('Cell Response')
            %title({[mouse,' - ',num2str(date),' - ',regexprep(Filename(16:end-4), '_','-')],criterion{q},['index opto = ',num2str(stimInfo.index(ii,dim_relevant(2)))]})
        end
        
        % plot mean cell response over all cells
        subplot(1,dim_index(2)+1,dim_index(2)+1)
        
        for ii=1:dim_index(2) % dim index 2 = laser
            a=[ii:dim_index(2):length(stimInfo.index)];
            hold  on
            if length(Ind_CellSelection{q})>1
                errorbar(stimInfo.index(a,1),mean(av_raster_poststim(a,Ind_CellSelection{q}),2),std(av_raster_poststim(a,Ind_CellSelection{q})')/sqrt(length(Ind_CellSelection{q})),'Color',(ii-1)/dim_index(2)*[1 1 1])
            else
                errorbar(stimInfo.index(a,1),av_raster_poststim(a,Ind_CellSelection{q}),sem_raster_poststim(a,Ind_CellSelection{q}),'Color',(ii-1)/dim_index(2)*[1 1 1])
            end
            xlabel('Sound amplitude (dB)')
            ylabel('Cell Response')
            title({'Average over all selected cells',['Number of cells = ',num2str(length(Ind_CellSelection{q})),'/',num2str(size(av_raster_poststim,2))],'black = laser OFF'})
        end
         
    end
end


%% plot cells with their stats in outline

fig(8)=figure;
 micronsPerPixel=exptInfo.recInfo.micronsPerPixel(1);
%       suptitle({[exptInfo.recDate,' - ', exptInfo.mouse, ' - ', regexprep(Filename(16:end-4), '_','-')], ...
%          relevant_info}) 
 imagesc(brighten(brighten(dat.mimg(:,:,2),1),0.5))
colormap gray
hold on
% c = colormap('jet');
 %   c = c(1:floor(size(c,1)/length(uModules{jj})):floor(size(c,1)/length(uModules{jj}))*length(uModules{jj}),:);
 %   colormap gray


 for kk = 1:size(spatialInfo.ROIs,2) %all cells grey outline
     h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',1,'edgecolor',[0.7 0.7 0.7]);
 end
 for kk=1:length(Ind_CellSelection{5}) %laser+, red inside
     h = patch(spatialInfo.ROIs{Ind_CellSelection{5}(kk)}(:,2),spatialInfo.ROIs{Ind_CellSelection{5}(kk)}(:,1),'y','facealpha',0.8,'edgealpha',0,'LineWidth',1,'facecolor',[1 0.6 0.6]);
 end
 for kk=1:length(Ind_CellSelection{6}) %laser-, blue inside
     h = patch(spatialInfo.ROIs{Ind_CellSelection{6}(kk)}(:,2),spatialInfo.ROIs{Ind_CellSelection{6}(kk)}(:,1),'y','facealpha',0.8,'edgealpha',0,'LineWidth',1,'facecolor',[0 0 1]);
 end
 for kk=1:length(Ind_CellSelection{3}) %sound+, yellow linewidth 3
     h = patch(spatialInfo.ROIs{Ind_CellSelection{3}(kk)}(:,2),spatialInfo.ROIs{Ind_CellSelection{3}(kk)}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',3,'edgecolor',[1 1 0]);
 end
 for kk=1:length(Ind_CellSelection{4}) %sound-, cyan linewidth 3
     h = patch(spatialInfo.ROIs{Ind_CellSelection{4}(kk)}(:,2),spatialInfo.ROIs{Ind_CellSelection{4}(kk)}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',3,'edgecolor',[0 1 1]);
 end
  if exist('RedCells')==1 % red cells, red outline
     for kk=1:length(Ind_CellSelection{2})
         h = patch(spatialInfo.ROIs{Ind_CellSelection{2}(kk)}(:,2),spatialInfo.ROIs{Ind_CellSelection{2}(kk)}(:,1),'r','facealpha',0,'edgealpha',1,'LineWidth',1,'edgecolor',[1 0 0]);
     end
 end
 
axis tight equal
set(gca,'TickDir','out','FontSize',14)
xlabel('Dorsal-ventral (pixels)')
ylabel('Rostral-caudal (pixels)')

%title({[num2str(date),'_',mouse,'_',num2str(Filename(16:end-4))], ...
        % ['n = ',num2str(size(raster,3)),' total cells'], ...
        % ['Sound+ = yellow : n = ',num2str(length(Ind_CellSelection{3})), '/', num2str(size(raster,3))], ...
        % ['Sound- = cyan : n = ',num2str(length(Ind_CellSelection{4})), '/', num2str(size(raster,3))], ...
        % ['Laser+ = red interior : n = ',num2str(length(Ind_CellSelection{5})), '/', num2str(size(raster,3))], ...
        % ['Laser- = blue interior : n = ',num2str(length(Ind_CellSelection{6})), '/', num2str(size(raster,3))], ...
        % ['Red cells = red : n = ',num2str(length(Ind_CellSelection{2})), '/', num2str(size(raster,3))] })

%% Sound responsive cell selection based on diff in means of baseline and stim evoked response
% Melanie's comment -  this is basically what the t-test does, so it is not
% that different.

close all;

Tdelay_test=[0:0.03:rasterInfo.postOnsetTime-Tdelay_windur];
 
soundMeanTrace_diffs_time=[];

for jj=1:size(raster,3)
    % caulculate stats with all Tdelay_test
    for kk=1:length(Tdelay_test)
        soundMeanTrace_diffs_time(1:7,kk)=abs(nanmean(mean_Trace([1,4,7,10,13,16,19],...
            round(fr*(preOnsetTime+Tdelay_test(kk))):round(fr*(preOnsetTime+Tdelay_test(kk)))+29,jj),2)-...
            nanmean(mean_Trace([1,4,7,10,13,16,19],1:round(fr*preOnsetTime),jj),2));
    end
    diff_mean_responsive_cells(jj)= max(max(soundMeanTrace_diffs_time,[],2))>0.1;
end

cellsSoundResponsive = find(diff_mean_responsive_cells==1);

meanResponse_fromMeanTrace = zeros(21,length(cellsSoundResponsive));
semResponse_fromMeanTrace = zeros(21,length(cellsSoundResponsive));
for jj=1:length(cellsSoundResponsive)
    for kk=1:length(Tdelay_test)
        meanResponse_Tdelay_fromMeanTrace(:,kk)=nanmean(mean_Trace(:,round(fr*(preOnsetTime+Tdelay_test(kk))):round(fr*(preOnsetTime+Tdelay_test(kk)))+29,...
            cellsSoundResponsive(jj)),2);
        semResponse_Tdelay_fromMeanTrace(:,kk)=nanstd(mean_Trace(:,round(fr*(preOnsetTime+Tdelay_test(kk))):round(fr*(preOnsetTime+Tdelay_test(kk)))+29,...
            cellsSoundResponsive(jj)),[],2)/sqrt(29);
    end
    [~,index_Tdelay] = max(meanResponse_Tdelay_fromMeanTrace,[],2);
    for ii=1:21  
        meanResponse_fromMeanTrace(ii,jj)= meanResponse_Tdelay_fromMeanTrace(ii,index_Tdelay(ii));
        semResponse_fromMeanTrace(ii,jj)= semResponse_Tdelay_fromMeanTrace(ii,index_Tdelay(ii));
    end
end

% After determining sound responsive cells plot response functions using raster_poststim variables

for i=1:length(Ind_CellSelection{3,1})
    figure();
    errorbar([0,30,50,60,70,80,90],av_raster_poststim([1,4,7,10,13,16,19],Ind_CellSelection{3,1}(i)),sem_raster_poststim([1,4,7,10,13,16,19],Ind_CellSelection{3,1}(i)));
    hold on;
    errorbar([0,30,50,60,70,80,90],av_raster_poststim([2,5,8,11,14,17,20],Ind_CellSelection{3,1}(i)),sem_raster_poststim([2,5,8,11,14,17,20],Ind_CellSelection{3,1}(i)));
    errorbar([0,30,50,60,70,80,90],av_raster_poststim([3,6,9,12,15,18,21],Ind_CellSelection{3,1}(i)),sem_raster_poststim([3,6,9,12,15,18,21],Ind_CellSelection{3,1}(i)));
end

% After determining sound responsive cells plot response functions using
% mean_Trace variables

% for i=1:length(cellsSoundResponsive)
% figure();
% errorbar([0,30,50,60,70,80,90],meanResponse_fromMeanTrace([1,4,7,10,13,16,19],i),semResponse_fromMeanTrace([1,4,7,10,13,16,19],i));
% hold on;
% errorbar([0,30,50,60,70,80,90],meanResponse_fromMeanTrace([2,5,8,11,14,17,20],i),semResponse_fromMeanTrace([2,5,8,11,14,17,20],i));
% errorbar([0,30,50,60,70,80,90],meanResponse_fromMeanTrace([3,6,9,12,15,18,21],i),semResponse_fromMeanTrace([3,6,9,12,15,18,21],i));
% end

% save('MeanTrace_threshold0.6_responseFunctions.mat','cellsSoundResponsive','av_raster_poststim','sem_raster_poststim','meanResponse_fromMeanTrace','semResponse_fromMeanTrace'); 


%% Test fitting of sigmoid to one of response functions

x = [0;30;50;60;70;80;90];
y = av_raster_poststim([1,4,7,10,13,16,19],cellsSoundResponsive(8));

ft = fittype('sigmoid(x,y0,ymax,x0,delta_x)');
f = fit(x,y,ft,'StartPoint',[y(1),y(end),60,(y(end)-y(1))/50]);

figure();
plot(f,x,y);
        
%% Prepping mat file for pca using \delta F/F0 rasters
%{
% Removing nans by interpolating

raster_noNaN = raster;
for ii = 1:size(raster_noNaN,1)
    for kk = 1:size(raster_noNaN,3)
        nanidxs = find(isnan(raster_noNaN(ii,:,kk)));
        if ~isempty(nanidxs)
            for jj = 1:length(nanidxs) 
                raster_noNaN(ii,nanidxs(jj),kk) = nanmean(raster_noNaN(ii,nanidxs(jj)-4:nanidxs(jj)+4,kk));
            end
        end
    end
end

raster_ordered_noNaN = raster_ordered;
for ii = 1:size(raster_ordered_noNaN,1)
    for kk = 1:size(raster_ordered_noNaN,3)
        nanidxs = find(isnan(raster_ordered_noNaN(ii,:,kk)));
        if ~isempty(nanidxs)
            for jj = 1:length(nanidxs) 
                raster_ordered_noNaN(ii,nanidxs(jj),kk) = nanmean(raster_ordered_noNaN(ii,nanidxs(jj)-4:nanidxs(jj)+4,kk));
            end
        end
    end
end

smoothed_raster = zeros(size(raster));
smoothed_raster_ordered = zeros(size(raster_ordered));
for ii = 1:size(raster,1)
    for kk = 1:size(raster,3)
        smoothed_raster(ii,:,kk) = smooth(raster_noNaN(ii,:,kk),5);
        smoothed_raster_ordered(ii,:,kk) = smooth(raster_ordered_noNaN(ii,:,kk),5);
    end
end

raster_smoothed_n_by_ct = zeros(210*size(raster,2),size(raster,3));
raster_ordered_smoothed_n_by_ct = zeros(210*size(raster,2),size(raster,3));
for i = 1:210
    max_vector = abs(max(smoothed_raster(i,:,:),[],2,'omitnan'));
    % max_vector(max_vector < 1) = 1; % using soft normalization like the paper neural population dynamics during reaching
    raster_smoothed_n_by_ct((i-1)*size(raster,2)+1:i*size(raster,2),:) = squeeze(smoothed_raster(i,:,:));%./reshape(max_vector,[1,size(raster,3)]);
end
raster_smoothed_n_by_ct = raster_smoothed_n_by_ct-nanmean(raster_smoothed_n_by_ct,1);
% raster_smoothed_n_by_ct = raster_smoothed_n_by_ct./abs(nanstd(raster_smoothed_n_by_ct,1));

for i = 1:210
    max_vector = abs(max(smoothed_raster_ordered(i,:,:),[],2,'omitnan'));
    % max_vector(max_vector < 1) = 1; % using soft normalization like the paper neural population dynamics during reaching
    temp = squeeze(smoothed_raster_ordered(i,:,:))./reshape(max_vector,[1,size(raster,3)]);
    raster_ordered_smoothed_n_by_ct((i-1)*size(raster,2)+1:i*size(raster,2),:) = squeeze(smoothed_raster_ordered(i,:,:));%./reshape(max_vector,[1,size(raster,3)]);
end
raster_ordered_smoothed_n_by_ct = raster_ordered_smoothed_n_by_ct-nanmean(raster_ordered_smoothed_n_by_ct,1);
% raster_ordered_smoothed_n_by_ct = raster_ordered_smoothed_n_by_ct./abs(nanstd(raster_ordered_smoothed_n_by_ct,1));

%% Prepping mat file for pca using averaged \delta F/F0 rasters

raster_ordered_mean_n_by_ct = zeros(prod(size(mean_Trace,1:2)),size(mean_Trace,3));

for i = 1:length(stimInfo.index)
    max_vector = abs(max(mean_Trace(i,:,:),[],2,'omitnan'));
    % max_vector(max_vector < 1) = 1; % using soft normalization like the paper neural population dynamics during reaching
    raster_ordered_mean_n_by_ct((i-1)*size(raster,2)+1:i*size(raster,2),:) = squeeze(mean_Trace(i,:,:));%./reshape(max_vector,[1,size(raster,3)]);
end
raster_ordered_mean_n_by_ct = raster_ordered_mean_n_by_ct-nanmean(raster_ordered_mean_n_by_ct,1);
% raster_ordered_mean_n_by_ct = raster_ordered_mean_n_by_ct./abs(nanstd(raster_ordered_mean_n_by_ct,1));

smoothed_mean_Trace = zeros(size(mean_Trace));
raster_ordered_smoothed_mean_n_by_ct = zeros(prod(size(smoothed_mean_Trace,1:2)),size(smoothed_mean_Trace,3));

for ii = 1:size(mean_Trace,1)
    for kk = 1:size(mean_Trace,3)
        % by default smoothing is over 5 points
        smoothed_mean_Trace(ii,:,kk) = smooth(mean_Trace(ii,:,kk),5); 
    end
end

for i = 1:length(stimInfo.index)
    max_vector = abs(max(smoothed_mean_Trace(i,:,:),[],2,'omitnan'));
    % max_vector(max_vector < 1) = 1; % using soft normalization like the paper neural population dynamics during reaching
    raster_ordered_smoothed_mean_n_by_ct((i-1)*size(raster,2)+1:i*size(raster,2),:) = squeeze(smoothed_mean_Trace(i,:,:));%./reshape(max_vector,[1,size(raster,3)]);
end
% raster_ordered_smoothed_mean_n_by_ct = raster_ordered_smoothed_mean_n_by_ct-nanmean(raster_ordered_smoothed_mean_n_by_ct,1);
% raster_ordered_smoothed_mean_n_by_ct = raster_ordered_smoothed_mean_n_by_ct./abs(nanstd(raster_ordered_smoothed_mean_n_by_ct,1));

%% Cell selection using Melanie's test on av_raster

Tdelay_test=[0.1:0.03:1.5];
Tdelay_windur=1; % window duration of 300ms for averaging wih firing rates.
alpha=0.01;
 
mean_Trace_Tdelay=[];
mean_Trace_tstats_time=[];
mean_Trace_pstats_time=[];

for jj=1:size(raster,3)
    % caulculate stats with all Tdelay_test
    for kk=1:length(Tdelay_test)
        [mean_Trace_tstats_time(jj,kk), mean_Trace_pstats_time(jj,kk)]=ttest(av_raster(1:round(fr*preOnsetTime),jj), av_raster(round(fr*(preOnsetTime+Tdelay_test(kk))):round(fr*(preOnsetTime+Tdelay_test(kk)))+29,jj),'Alpha',alpha);
    end
    mean_Trace_Tdelay(jj)=min(Tdelay_test(find(mean_Trace_pstats_time(jj,:)==min(mean_Trace_pstats_time(jj,:)))));

    % compute stats, average and std with chosen Tdelay
    [mean_Trace_tstats(jj), mean_Trace_pstats(jj)]=ttest(av_raster(1:round(fr*preOnsetTime),jj),av_raster(round(fr*(preOnsetTime+mean_Trace_Tdelay(jj))):round(fr*(preOnsetTime+mean_Trace_Tdelay(jj)))+29,jj),'Alpha',alpha);

end

%% Alternative cell selection based on diff in means of baseline and stim evoked response

Tdelay_test=[0.1:0.03:1.5];
 
mean_Trace_diffs_Tdelay=[];
mean_Trace_diffs_time=[];

for jj=1:size(raster,3)
    % caulculate stats with all Tdelay_test
    for kk=1:length(Tdelay_test)
        mean_Trace_diffs_time(jj,kk)=abs(nanmean(av_raster(round(fr*(preOnsetTime+Tdelay_test(kk))):round(fr*(preOnsetTime+Tdelay_test(kk)))+29,jj))-nanmean(av_raster(1:round(fr*preOnsetTime),jj)));
    end
    diff_mean_responsive_cells(jj)= max(mean_Trace_diffs_time(jj,:))>0.1;
end
%}