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

%clear variables

load('/Users/janaki/Dropbox/project_with_melanie/DataForJanaki/ImagingDatasetsAndCorrespondingInfo_DirectlyFromExpts/MT065_20190905_FRA_OptoStim_03.mat');
load('/Users/janaki/Dropbox/project_with_melanie/DataForJanaki/ImagingDatasetsAndCorrespondingInfo_DirectlyFromExpts/MT065_20190905_RedCells_SlopeRemoval_Movie.mat');
load('/Users/janaki/Dropbox/project_with_melanie/DataForJanaki/ImagingDatasetsAndCorrespondingInfo_DirectlyFromExpts/F_MT065_20190905MT065_tifStacks_plane1_proc.mat');

dir = '/Users/janaki/Dropbox/project_with_melanie/DataForJanaki/Across_many_datasets/vip/';
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

close all;

%% sort for every condition depending on sound response and red intensity and show rasters

% Janaki: this is the section which looks at which cells are
% stim-responsive. I tried a few options out, the ones I'm not using
% anymore more are commented out, but I leave them there so you see what I
% tried. 


%%%% Option 2 : varying interval then take min of pstat t-test; adapated from
%%%% above; no correction to t-test because I then take only one window.
%%%% (ask melanie: what is a t-test correction?)

Tdelay_test=[0.1:0.03:1.5];
Tdelay_windur=1; % window duration of 1s for averaging.
alpha=0.01;
 
Tdelay=[];
tstats_time=[];
pstats_time=[];
zscore_time=[];

% calculate stats for all cells, all stims, all windows
for jj=1:size(raster,3)
for ii=1:length(stimInfo.index)  
    % calculate stats with all Tdelay_test
    for kk=1:length(Tdelay_test)
     [tstats_time(ii,jj,kk) pstats_time(ii,jj,kk)]=ttest(nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,1:round(fr*preOnsetTime),jj),2),nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,round(fr*(preOnsetTime+Tdelay_test(kk))):round(fr*(preOnsetTime+Tdelay_test(kk)+Tdelay_windur)),jj),2),'Alpha',alpha);
     zscore(ii,jj,kk)=(mean(nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,1:round(fr*preOnsetTime),jj),2))-mean(nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,round(fr*(preOnsetTime+Tdelay_test(kk))):round(fr*(preOnsetTime+Tdelay_test(kk)+Tdelay_windur)),jj),2)))/sqrt(0.5*(std(nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,1:round(fr*preOnsetTime),jj),2))^2+std(nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,round(fr*(preOnsetTime+Tdelay_test(kk))):round(fr*(preOnsetTime+Tdelay_test(kk)+Tdelay_windur)),jj),2))^2)^2);
    % average pop response for each time window  
    b2=squeeze(nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,round(fr*(preOnsetTime+Tdelay_test(kk))):round(fr*(preOnsetTime+Tdelay_test(kk)+Tdelay_windur)),jj),2));
    av_raster_ordered_time(repeats*(ii-1)+1:repeats*ii,jj,kk) = b2;
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
    av_raster_ordered(repeats*(ii-1)+1:repeats*ii,:)=av_raster_ordered_time(repeats*(ii-1)+1:repeats*ii,:,...
        find(Nsignifcells_allstims_time(:)==max(Nsignifcells_allstims_time(:)),1,'first'));
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

%% Which cells are responsive to sound? 
 % first pass: count all cells with 1 / -1 in no light 
 cellsActiveInNoLight = find(sum(abs(Cells_ActivityChange2([4,7,10,13,16,19],:)),1)>0);
 cellsActiveInMidLight = find(sum(abs(Cells_ActivityChange2([4,7,10,13,16,19]+1,:)),1)>0);
 cellsActiveInHighLight = find(sum(abs(Cells_ActivityChange2([4,7,10,13,16,19]+2,:)),1)>0);
 
 % discriminability: how many cells are responsive to only a selective
 % number of sound stimuli. by selective we mean 1/2. do these cells
 % exhibit non-monotonic behavior? 
 % (caveat: finding overlap between matlab and python cells is not very obvious)
 discriminableCellsNoLight = [find(sum(abs(Cells_ActivityChange2([4,7,10,13,16,19],:)),1)==1)';
 find(sum(abs(Cells_ActivityChange2([4,7,10,13,16,19],:)),1)==2)'];
 discriminableCellsMidLight = [find(sum(abs(Cells_ActivityChange2([4,7,10,13,16,19]+1,:)),1)==1)';
 find(sum(abs(Cells_ActivityChange2([4,7,10,13,16,19]+1,:)),1)==2)'];
 discriminableCellsHighLight = [find(sum(abs(Cells_ActivityChange2([4,7,10,13,16,19]+2,:)),1)==1)';
 find(sum(abs(Cells_ActivityChange2([4,7,10,13,16,19]+2,:)),1)==2)'];

 % detectability: how many cells are responsive to a majority of sound
 % stimuli. by majority we mean 5/6. do these cells exhibit monotonic
 % behavior?
 detectableCellsNoLight = find(sum(abs(Cells_ActivityChange2([4,7,10,13,16,19],:)),1)>4);
 detectableCellsMidLight = find(sum(abs(Cells_ActivityChange2([4,7,10,13,16,19]+1,:)),1)>4);
 detectableCellsHighLight = find(sum(abs(Cells_ActivityChange2([4,7,10,13,16,19]+2,:)),1)>4);
 
load('DetectabilityVsDiscriminability_vip.mat');
dataset = 15;
DetectableCellsHighLight{dataset} = detectableCellsHighLight;
DetectableCellsMidLight{dataset} = detectableCellsMidLight;
DetectableCellsNoLight{dataset} = detectableCellsNoLight;
DiscriminableCellsNoLight{dataset} = discriminableCellsNoLight;
DiscriminableCellsMidLight{dataset} = discriminableCellsMidLight;
DiscriminableCellsHighLight{dataset} = discriminableCellsHighLight;
save('DetectabilityVsDiscriminability_vip','DetectableCellsHighLight','DetectableCellsMidLight','DetectableCellsNoLight',...
    'DiscriminableCellsHighLight','DiscriminableCellsMidLight','DiscriminableCellsNoLight');

%% find clusters of cells based on proximity to red cells

% check which red cells respond to light only (no sound), we study each 
% light condition separately for reponse functions and study both the
% conditions together for the correlation functions

responsiveRedcells = [];
responsiveRedcellstoMidLight = [];
responsiveRedcellstoHighLight = [];
for i = 1:length(Ind_cells_red)
    if any(abs(Cells_ActivityChange2([2,3],Ind_cells_red(i))))
        responsiveRedcells = [responsiveRedcells, ...
            Ind_cells_red(i)*any(abs(Cells_ActivityChange2([2,3],Ind_cells_red(i))))];
        responsiveRedcellstoMidLight = [responsiveRedcellstoMidLight, ...
            Ind_cells_red(i)*any(abs(Cells_ActivityChange2(2,Ind_cells_red(i))))];
        responsiveRedcellstoHighLight = [responsiveRedcellstoHighLight, ...
            Ind_cells_red(i)*any(abs(Cells_ActivityChange2(3,Ind_cells_red(i))))];
    end
end
responsiveRedcellstoMidLight = nonzeros(responsiveRedcellstoMidLight);
responsiveRedcellstoHighLight = nonzeros(responsiveRedcellstoHighLight);

% find ROIs of active red cells
Redcells_ROI = zeros(length(responsiveRedcells),2);
for i = 1:size(Redcells_ROI,1)
    Redcells_ROI(i,:) = mean(spatialInfo.ROIs{responsiveRedcells(i)});
end

%find closest red cell to all cells
Allcells_ROI = zeros(size(raster,3),2);
for i = 1:size(Allcells_ROI,1)
    Allcells_ROI(i,:) = mean(spatialInfo.ROIs{i});
end

DistfromRedcells = zeros(length(Allcells_ROI),length(Redcells_ROI));
for i = 1:length(Allcells_ROI)
    DistfromRedcells(i,:) = vecnorm(Redcells_ROI-Allcells_ROI(i,:),2,2);
end

[closestRedcell(:,1),closestRedcell(:,2)] = min(DistfromRedcells,[],2);

%find distance between red cells
DistbetRedcells = zeros(length(responsiveRedcells),length(responsiveRedcells));
for i = 1:size(DistbetRedcells,1)
    DistbetRedcells(i,:) = vecnorm(Redcells_ROI-Redcells_ROI(i,:),2,2);
end

% plotting red cells and other cells that are within its associated
% cluster

fig1 = figure();
micronsPerPixel=exptInfo.recInfo.micronsPerPixel(1);
imagesc(brighten(brighten(dat.mimg(:,:,2),1),0.5))
colormap gray
hold on

for ii = 1:length(closestRedcell)
    h = patch(spatialInfo.ROIs{ii}(:,2),spatialInfo.ROIs{ii}(:,1),'y','facealpha',0.8,'edgealpha',0,...
        'LineWidth',1,'facecolor',[closestRedcell(ii,2)/length(responsiveRedcells) closestRedcell(ii,2)/length(responsiveRedcells) 1]);
end

for ii = 1:length(responsiveRedcells)
    h = patch(spatialInfo.ROIs{responsiveRedcells(ii)}(:,2),spatialInfo.ROIs{responsiveRedcells(ii)}(:,1),'y',...
        'facealpha',0.8,'edgealpha',0,'LineWidth',1,...
        'facecolor',[1 (length(responsiveRedcells)-ii)/length(responsiveRedcells) (length(responsiveRedcells)-ii)/length(responsiveRedcells)]);
    text(spatialInfo.ROIs{responsiveRedcells(ii)}(:,2), spatialInfo.ROIs{responsiveRedcells(ii)}(:,1)+ 0.3, ...
        num2str(responsiveRedcells(ii)), 'Color', 'g', 'FontSize', 5);
end
%savefig(fig1, 'Cell layout');
%saveas(fig1,'Cell layout.png')
%close all;

%% average max distance for the red cell regions 

% this is to document how large the red cell regions are to compare them
% across vip and som datasets

cellIdxs = 1:length(closestRedcell);
%fileID = fopen('maxDistance1.txt','a+');

for ii = 1:length(responsiveRedcells)
   closestCells = cellIdxs(closestRedcell(:,2)==ii);
    [distclosestCells, sortDistIdx] = sort(closestRedcell(closestCells,1));
    max_dist = max(distclosestCells);
    %fprintf(fileID, strcat(string(max_dist), '\n'));
end

%fclose(fileID);

%% create influence tables for distance and sound

% influence = (activity in light + sound - activity in only sound)/(red cell activity in light without sound)

influencemidSoundDist = {};
influencehighSoundDist = {};%cell(length(responsiveRedcells),1);

for ii = 1:length(responsiveRedcells)
    fig2 = figure(); ax2=axes(fig2); 
    closestCells = cellIdxs(closestRedcell(:,2)==ii);
    closestCells(closestCells==responsiveRedcells(ii)) = [];
    [distclosestCells, sortDistIdx] = sort(closestRedcell(closestCells,1));
    influencemidSoundDist{ii} = (av_raster_poststim([2,5,8,11,14,17,20],closestCells)-av_raster_poststim([1,4,7,10,13,16,19],closestCells));
    influencemidSoundDist{ii} = influencemidSoundDist{ii}/(av_raster_poststim(2,responsiveRedcells(ii))-av_raster_poststim(1,responsiveRedcells(ii)));
    for i=1:dim_index(1)
        plot(ax2, distclosestCells, influencemidSoundDist{ii}(i,sortDistIdx),...
            'color',[(dim_index(1)-i)/dim_index(1),(dim_index(1)-i)/dim_index(1),(dim_index(1)-i)/dim_index(1)]);
        hold on;
    end    
    title(join(['Redcell ',num2str(responsiveRedcells(ii))]));
    xlabel('Distance'); ylabel('Light = 0.34 influence');
    %savefig(fig2, join(['Mid light response ',num2str(responsiveRedcells(ii))]));
    %saveas(fig2, strcat(join(['Mid light response ',num2str(responsiveRedcells(ii))]),'.png'))
    close all;
end

for ii = 1:length(responsiveRedcells)
    fig3 = figure(); ax3=axes(fig3);
    closestCells = cellIdxs(closestRedcell(:,2)==ii);
    closestCells(closestCells==responsiveRedcells(ii)) = [];
    [distclosestCells, sortDistIdx] = sort(closestRedcell(closestCells,1));
    influencehighSoundDist{ii} = (av_raster_poststim([3,6,9,12,15,18,21],closestCells)-av_raster_poststim([1,4,7,10,13,16,19],closestCells));
    influencehighSoundDist{ii} = influencehighSoundDist{ii}/(av_raster_poststim(3,responsiveRedcells(ii))-av_raster_poststim(1,responsiveRedcells(ii)));
    for i=1:dim_index(1)
        plot(ax3, distclosestCells, influencehighSoundDist{ii}(i,sortDistIdx),...
            'color',[(dim_index(1)-i)/dim_index(1),(dim_index(1)-i)/dim_index(1),(dim_index(1)-i)/dim_index(1)]);
        hold on;
    end    
    title(join(['Redcell ',num2str(responsiveRedcells(ii))]));
    xlabel('Distance'); ylabel('Light = 0.7 influence');
    %savefig(fig3, join(['High light response ',num2str(responsiveRedcells(ii))]));
    %saveas(fig3, strcat(join(['High light response ',num2str(responsiveRedcells(ii))]),'.png'))
    close all;
end

%% decay of influence function (to light modulation in the presence of sound) with distance 

% This is collective dataset info stored in 'Across many datasets'
% subfolders. An assumption here is that the closest red cell has the greatest
% effect on any given cell which is why we cluster cells based on their
% closest red cell. We are also assuming that its effect is greater than
% that of the other red cells. 

close all;

% load('influencehighSoundDist_allsounds.mat');
% load('influencemidSoundDist_allsounds.mat');

influencehighSoundDist_allsounds = []; influencehighSoundDist_distance = [];
influencemidSoundDist_allsounds = []; influencemidSoundDist_distance = [];

for midii = 1:length(responsiveRedcellstoMidLight)
    ii = find(responsiveRedcells==responsiveRedcellstoMidLight(midii),1);
    closestCells = cellIdxs(closestRedcell(:,2)==ii);
    [distclosestCells, sortDistIdx] = sort(closestRedcell(closestCells,1));
    influencemidSoundDist{ii} = (av_raster_poststim([2,5,8,11,14,17,20],closestCells(sortDistIdx))-av_raster_poststim([1,4,7,10,13,16,19],closestCells(sortDistIdx)));
    influencemidSoundDist{ii} = influencemidSoundDist{ii}/(av_raster_poststim(2,responsiveRedcellstoMidLight(midii))-av_raster_poststim(1,responsiveRedcellstoMidLight(midii)));
    influencemidSoundDist_allsounds = [influencemidSoundDist_allsounds, influencemidSoundDist{ii}];       
    influencemidSoundDist_distance = [influencemidSoundDist_distance, distclosestCells'];
end

for highii = 1:length(responsiveRedcellstoHighLight)
    ii = find(responsiveRedcells==responsiveRedcellstoHighLight(highii),1);
    closestCells = cellIdxs(closestRedcell(:,2)==ii);
    [distclosestCells, sortDistIdx] = sort(closestRedcell(closestCells,1));
    influencehighSoundDist{ii} = (av_raster_poststim([3,6,9,12,15,18,21],closestCells(sortDistIdx))-av_raster_poststim([1,4,7,10,13,16,19],closestCells(sortDistIdx)));
    influencehighSoundDist{ii} = influencehighSoundDist{ii}/(av_raster_poststim(3,responsiveRedcellstoHighLight(highii))-av_raster_poststim(1,responsiveRedcellstoHighLight(highii)));
    influencehighSoundDist_allsounds = [influencehighSoundDist_allsounds, influencehighSoundDist{ii}];       
    influencehighSoundDist_distance = [influencehighSoundDist_distance, distclosestCells'];
end


%% decay of activity (to modulation by both light and sound) with distance

% This is collective dataset info stored in 'Across many datasets'
% subfolders. The control expt is calculating activity only in the presence 
% of sound. Note that influence function as calculated above is different 
% from just activity.

% load('activityhighLightDist_allsounds.mat');
% load('activitymidLightDist_allsounds.mat');
% load('activitynoLightDist_forMid_allsounds.mat');
% load('activitynoLightDist_forHigh_allsounds.mat');

activitynoLightDist_forMid_allsounds = []; activitynoLightDist_forMid_distance = [];
activitynoLightDist_forHigh_allsounds = []; activitynoLightDist_forHigh_distance = [];
activitymidLightDist_allsounds = []; activitymidLightDist_distance = [];
activityhighLightDist_allsounds = []; activityhighLightDist_distance = [];

for highii = 1:length(responsiveRedcellstoHighLight)
    ii = find(responsiveRedcells==responsiveRedcellstoHighLight(highii),1);
    closestCells = cellIdxs(closestRedcell(:,2)==ii);
    [distclosestCells, sortDistIdx] = sort(closestRedcell(closestCells,1));
    disp(length(closestCells))
    activityhighLightDist_allsounds = [activityhighLightDist_allsounds, av_raster_poststim([3,6,9,12,15,18,21],closestCells(sortDistIdx))]; 
    activitynoLightDist_forHigh_allsounds = [activitynoLightDist_forHigh_allsounds, av_raster_poststim([1,4,7,10,13,16,19],closestCells(sortDistIdx))];     
    activityhighLightDist_distance = [activityhighLightDist_distance, distclosestCells']; 
    activitynoLightDist_forHigh_distance = [activitynoLightDist_forHigh_distance, distclosestCells'];  
end

for midii = 1:length(responsiveRedcellstoMidLight)
    ii = find(responsiveRedcells==responsiveRedcellstoMidLight(midii),1);
    closestCells = cellIdxs(closestRedcell(:,2)==ii);
    [distclosestCells, sortDistIdx] = sort(closestRedcell(closestCells,1));
    disp(length(closestCells))
    activitymidLightDist_allsounds = [activitymidLightDist_allsounds, av_raster_poststim([2,5,8,11,14,17,20],closestCells(sortDistIdx))]; 
    activitynoLightDist_forMid_allsounds = [activitynoLightDist_forMid_allsounds, av_raster_poststim([1,4,7,10,13,16,19],closestCells(sortDistIdx))];     
    activitymidLightDist_distance = [activitymidLightDist_distance, distclosestCells']; 
    activitynoLightDist_forMid_distance = [activitynoLightDist_forMid_distance, distclosestCells']; 
end

% save('influencehighSoundDist_allsounds.mat', 'influencehighSoundDist_allsounds', 'influencehighSoundDist_distance');
% save('influencemidSoundDist_allsounds.mat','influencemidSoundDist_allsounds', 'influencemidSoundDist_distance');
% save('activityhighLightDist_allsounds.mat','activityhighLightDist_allsounds', 'activityhighLightDist_distance');
% save('activitymidLightDist_allsounds.mat','activitymidLightDist_allsounds', 'activitymidLightDist_distance');
% save('activitynoLightDist_forMid_allsounds.mat','activitynoLightDist_forMid_allsounds', 'activitynoLightDist_forMid_distance');
% save('activitynoLightDist_forHigh_allsounds.mat','activitynoLightDist_forHigh_allsounds', 'activitynoLightDist_forHigh_distance');


%% Create correlation tables for sound, light conditions for neurons within cluster regions of red cells
% This can either take the form of corr coeffs / be cross correlations of
% which I take the maximum/minimum / be cross correlation of which I take the value at time lag tau = 0. 
% The second method is useful when activities of neurons change with a delay. The third when you are looking for simulataneous activity.   

corrCells = cell(length(responsiveRedcells),1);
for ii = 1:length(responsiveRedcells)
    closestCells = cellIdxs(closestRedcell(:,2)==ii);
    corrCells{ii} = zeros(length(closestCells)*(length(closestCells)-1),22);
    for jj = 1:length(closestCells)
        othercloseCells = closestCells;
        othercloseCells(jj) = [];
        %othercloseCells(closestCells==responsiveRedcells(ii)) = [];
        corrCells{ii}(length(othercloseCells)*(jj-1)+1:length(othercloseCells)*jj,1) ...
            = vecnorm(Allcells_ROI(othercloseCells,:)-Allcells_ROI(closestCells(jj),:),2,2);
        for kk = 1:length(othercloseCells)
            for cc = 1:21
                %coeff = corrcoef(mean_Trace(cc,:,othercloseCells(kk)),mean_Trace(cc,:,closestCells(jj))); 
                %corrCells{ii}(length(othercloseCells)*(jj-1)+kk,cc+1) = coeff(1,2);
                correlationMeanTraces = xcorr(mean_Trace(cc,:,othercloseCells(kk)),...
                    mean_Trace(cc,:,closestCells(jj)),'coeff'); 
                %[coeff,corrIdx] = max(abs(correlationMeanTraces));
                %corrCells{ii}(length(othercloseCells)*(jj-1)+kk,cc+1) = coeff.*sign(correlationMeanTraces(corrIdx));
                corrCells{ii}(length(othercloseCells)*(jj-1)+kk,cc+1) = correlationMeanTraces(size(mean_Trace,2));
            end
        end
    end
end

%% variation of correlation with distance for each light condition considering cells within each cluster
% Note that this subsection is used for calculating correlations across
% sound conditions (Mean,sem comparing xcorr(tau=0) between cells in each redcell neighbourhood diff sounds), 
% light conditions (Mean,sem comparing xcorr(tau=0) between cells in each redcell neighbourhood different lights)
% and to look for the nonlinearity that may arise when combining both 
% light and sound conditions (Mean,sem comparing xcorr(tau=0) between cells in each redcell neighbourhood sound-high light nonlinearity)

close all;

fig6 = figure(); ax6=axes(fig6);
noSoundLightActivitycorr = [];
noSoundLightDistBetCells = [];
for ii = 1:length(responsiveRedcells)
    [distclosestCells, sortDistIdx] = sort(corrCells{ii}(:,1));
    for i=1
        plot(ax6, distclosestCells, corrCells{ii}(sortDistIdx,i+1),'r.');
        noSoundLightActivitycorr = [noSoundLightActivitycorr; corrCells{ii}(sortDistIdx,i+1)];
        noSoundLightDistBetCells = [noSoundLightDistBetCells; distclosestCells];
        hold on;
    end       
end
xlabel('Distance'); ylabel('Correlation coefficients');
ylim([-1,1]);
%savefig(fig6, 'No light correlation between cells in each redcell neighbourhood');
%saveas(fig6, 'No light correlation between cells in each redcell neighbourhood','png');
%close all;

fig7 = figure(); ax7=axes(fig7);
noSoundActivitycorr = [];
noSoundDistBetCells = [];
for ii = 1:length(responsiveRedcells)
    [distclosestCells, sortDistIdx] = sort(corrCells{ii}(:,1));
    for i=3
        plot(ax7, distclosestCells, corrCells{ii}(sortDistIdx,i+1),'g.');
        noSoundActivitycorr = [noSoundActivitycorr; corrCells{ii}(sortDistIdx,i+1)];
        noSoundDistBetCells = [noSoundDistBetCells; distclosestCells];
        hold on;
    end   
end
xlabel('Distance'); ylabel('Correlation coefficients');
ylim([-1,1]);
%savefig(fig7, 'Mid light correlation between cells in each redcell neighbourhood');
%saveas(fig7, 'Mid light correlation between cells in each redcell neighbourhood','png');
%close all;

fig8 = figure(); ax8=axes(fig8);
noLightActivitycorr = [];
noLightDistBetCells = [];
for ii = 1:length(responsiveRedcells)
    [distclosestCells, sortDistIdx] = sort(corrCells{ii}(:,1));
    for i=[4,7,10,13,16,19]
        plot(ax8, distclosestCells, corrCells{ii}(sortDistIdx,i+1),'k.');
        noLightActivitycorr = [noLightActivitycorr; corrCells{ii}(sortDistIdx,i+1)];
        noLightDistBetCells = [noLightDistBetCells; distclosestCells];
        hold on;
    end   
end
xlabel('Distance'); ylabel('Correlation coefficients');
ylim([-1,1]);
%savefig(fig8, 'High light correlation between cells in each redcell neighbourhood');
%saveas(fig8, 'High light correlation between cells in each redcell neighbourhood','png');
%close all;

fig9 = figure(); ax9=axes(fig9);
SoundLightActivitycorr = [];
SoundLightDistBetCells = [];
for ii = 1:length(responsiveRedcells)
    [distclosestCells, sortDistIdx] = sort(corrCells{ii}(:,1));
    for i=[6,9,12,15,18,21]
        plot(ax9, distclosestCells, corrCells{ii}(sortDistIdx,i+1),'k.');
        SoundLightActivitycorr = [SoundLightActivitycorr; corrCells{ii}(sortDistIdx,i+1)];
        SoundLightDistBetCells = [SoundLightDistBetCells; distclosestCells];
        hold on;
    end   
end
xlabel('Distance'); ylabel('Correlation coefficients');
ylim([-1,1]);
%savefig(fig8, 'High light correlation between cells in each redcell neighbourhood');
%saveas(fig8, 'High light correlation between cells in each redcell neighbourhood','png');
close all;

minNoSoundLightDistBetCells = min(noSoundLightDistBetCells);
maxNoSoundLightDistBetCells = max(noSoundLightDistBetCells);
rangeNoSoundLightDistBetCells = [floor(minNoSoundLightDistBetCells):10:maxNoSoundLightDistBetCells+10];
binnedNoSoundLightDistBetCells = discretize(noSoundLightDistBetCells,rangeNoSoundLightDistBetCells);
uniqueNoSoundLightDistBins = unique(binnedNoSoundLightDistBetCells);

minSoundLightDistBetCells = min(SoundLightDistBetCells);
maxSoundLightDistBetCells = max(SoundLightDistBetCells);
rangeSoundLightDistBetCells = [floor(minSoundLightDistBetCells):10:maxSoundLightDistBetCells+10];
binnedSoundLightDistBetCells = discretize(SoundLightDistBetCells,rangeSoundLightDistBetCells);
uniqueSoundLightDistBins = unique(binnedSoundLightDistBetCells);

mean_noSoundLightActivitycorr = []; sem_noSoundLightActivitycorr = [];
mean_noSoundActivitycorr = []; sem_noSoundActivitycorr = [];
mean_noLightActivitycorr = []; sem_noLightActivitycorr = [];
mean_SoundLightActivitycorr = []; sem_SoundLightActivitycorr = [];

figbox = figure(); axbox=axes(figbox); hold on;

for ii = 1:length(uniqueNoSoundLightDistBins)
    mean_noSoundLightActivitycorr = [mean_noSoundLightActivitycorr, mean(noSoundLightActivitycorr(binnedNoSoundLightDistBetCells==ii))];
    sem_noSoundLightActivitycorr = [sem_noSoundLightActivitycorr, ...
        std(noSoundLightActivitycorr(binnedNoSoundLightDistBetCells==ii))/sqrt(sum(binnedNoSoundLightDistBetCells==ii))];
    
    mean_noSoundActivitycorr = [mean_noSoundActivitycorr, mean(noSoundActivitycorr(binnedNoSoundLightDistBetCells==ii))];
    sem_noSoundActivitycorr = [sem_noSoundActivitycorr, std(noSoundActivitycorr(binnedNoSoundLightDistBetCells==ii))/sqrt(sum(binnedNoSoundLightDistBetCells==ii))];
end

for ii = 1:length(uniqueSoundLightDistBins)
    mean_noLightActivitycorr = [mean_noLightActivitycorr, mean(noLightActivitycorr(binnedSoundLightDistBetCells==ii))];
    sem_noLightActivitycorr = [sem_noLightActivitycorr, std(noLightActivitycorr(binnedSoundLightDistBetCells==ii))/sqrt(sum(binnedSoundLightDistBetCells==ii))];
    
    mean_SoundLightActivitycorr = [mean_SoundLightActivitycorr, mean(SoundLightActivitycorr(binnedSoundLightDistBetCells==ii))];
    sem_SoundLightActivitycorr = [sem_SoundLightActivitycorr, std(SoundLightActivitycorr(binnedSoundLightDistBetCells==ii))/sqrt(sum(binnedSoundLightDistBetCells==ii))];
end

plot(rangeNoSoundLightDistBetCells(1:23), mean_noSoundLightActivitycorr(1:23), 'k'); 
plot(rangeNoSoundLightDistBetCells(1:23), mean_noSoundLightActivitycorr(1:23) + sem_noSoundLightActivitycorr(1:23), 'k--');
plot(rangeNoSoundLightDistBetCells(1:23), mean_noSoundLightActivitycorr(1:23) - sem_noSoundLightActivitycorr(1:23), 'k--');

plot(rangeNoSoundLightDistBetCells(1:23), mean_noSoundActivitycorr(1:23), 'b'); 
plot(rangeNoSoundLightDistBetCells(1:23), mean_noSoundActivitycorr(1:23) + sem_noSoundActivitycorr(1:23), 'b--');
plot(rangeNoSoundLightDistBetCells(1:23), mean_noSoundActivitycorr(1:23) - sem_noSoundActivitycorr(1:23), 'b--');

plot(rangeSoundLightDistBetCells(1:23), mean_noLightActivitycorr(1:23), 'r'); 
plot(rangeSoundLightDistBetCells(1:23), mean_noLightActivitycorr(1:23) + sem_noLightActivitycorr(1:23), 'r--');
plot(rangeSoundLightDistBetCells(1:23), mean_noLightActivitycorr(1:23) - sem_noLightActivitycorr(1:23), 'r--');

plot(rangeSoundLightDistBetCells(1:23), mean_SoundLightActivitycorr(1:23), 'color', [1,0,1]); 
plot(rangeSoundLightDistBetCells(1:23), mean_SoundLightActivitycorr(1:23) + sem_SoundLightActivitycorr(1:23), 'color', [1,0,1], 'Linestyle', '--');
plot(rangeSoundLightDistBetCells(1:23), mean_SoundLightActivitycorr(1:23) - sem_SoundLightActivitycorr(1:23), 'color', [1,0,1], 'Linestyle', '--');

xlabel('Distance(um)'); ylabel('Correlation coefficients');
xlim([0,170]);
set(gca,'FontSize',20);
set(gcf, 'PaperUnits', 'inches');
x_width=10 ;y_width=8;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
savefig(figbox, 'Mean,sem comparing xcorr(tau=0) between cells in each redcell neighbourhood sound-high light nonlinearity');
saveas(figbox, 'Mean,sem comparing xcorr(tau=0) between cells in each redcell neighbourhood sound-high light nonlinearity','png');

%% Create correlation tables for sound, light conditions for neurons across cluster regions of red cells
% This can either take the form of corr coeffs / be cross correlations of
% which I take the maximum/minimum / be cross correlation of which I take the value at time lag tau = 0. 
% The second method is useful when activities of neurons change with a delay. The third when you are looking for simulataneous activity.   
%{
corrCellsAcrossRegs = cell(length(responsiveRedcells),1);
for ii = 1:length(responsiveRedcells)
    cellsReg1 = cellIdxs(closestRedcell(:,2)==ii);
    for jj = ii+1:length(responsiveRedcells)
        cellsReg2 = cellIdxs(closestRedcell(:,2)==jj);
        corrCellsAcrossRegs{ii,jj} = zeros(length(cellsReg1)*length(cellsReg2),22);
        for kk = 1:length(cellsReg1)
            corrCellsAcrossRegs{ii,jj}(length(cellsReg2)*(kk-1)+1:length(cellsReg2)*kk,1) ...
                = vecnorm(Allcells_ROI(cellsReg2,:)-Allcells_ROI(cellsReg1(kk),:),2,2);
            for ll = 1:length(cellsReg2)
                for cc = 1:21
                    %coeff = corrcoef(mean_Trace(cc,:,cellsReg2(ll)),mean_Trace(cc,:,cellsReg1(kk))); 
                    %corrCellsAcrossRegs{ii}(length(cellsReg2)*(kk-1)+ll,cc+1) = coeff(1,2);
                    correlationMeanTraces = xcorr(mean_Trace(cc,:,cellsReg2(ll)), mean_Trace(cc,:,cellsReg1(kk)),'coeff'); 
                    %[coeff,corrIdx] = max(abs(correlationMeanTraces));
                    %corrCellsAcrossRegs{ii}(length(cellsReg2)*(kk-1)+ll,cc+1) = coeff.*sign(correlationMeanTraces(corrIdx));
                    corrCellsAcrossRegs{ii,jj}(length(cellsReg2)*(kk-1)+ll,cc+1) = correlationMeanTraces(size(mean_Trace,2));
                end
            end
        end
    end
end

%% variation of correlation with distance for each light condition considering cells within each cluster
% note the change in for loops, I took out with sound expts

noLightActivitycorrAcrossRegs = [];
distBetCellsAcrossRegs = [];
for ii = 1:size(corrCellsAcrossRegs,1)
    for jj = ii+1:size(corrCellsAcrossRegs,2)
        [distclosestCells, sortDistIdx] = sort(corrCellsAcrossRegs{ii,jj}(:,1));
        for i=[1,4,7,10,13,16,19]
            noLightActivitycorrAcrossRegs = [noLightActivitycorrAcrossRegs;  corrCellsAcrossRegs{ii,jj}(sortDistIdx,i+1)];
            distBetCellsAcrossRegs = [distBetCellsAcrossRegs; distclosestCells];
        end   
    end
end

midLightActivitycorrAcrossRegs = [];
for ii = 1:size(corrCellsAcrossRegs,1)
    for jj = ii+1:size(corrCellsAcrossRegs,2)
        [distclosestCells, sortDistIdx] = sort(corrCellsAcrossRegs{ii,jj}(:,1));
        for i=[2,5,8,11,14,17,20]
            midLightActivitycorrAcrossRegs = [midLightActivitycorrAcrossRegs;  corrCellsAcrossRegs{ii,jj}(sortDistIdx,i+1)];
        end   
    end
end

highLightActivitycorrAcrossRegs = [];
for ii = 1:size(corrCellsAcrossRegs,1)
    for jj = ii+1:size(corrCellsAcrossRegs,2)
        [distclosestCells, sortDistIdx] = sort(corrCellsAcrossRegs{ii,jj}(:,1));
        for i=[3,6,9,12,15,18,21]
            highLightActivitycorrAcrossRegs = [highLightActivitycorrAcrossRegs;  corrCellsAcrossRegs{ii,jj}(sortDistIdx,i+1)];
        end   
    end
end

minDistBetCellsAcrossRegs = min(distBetCellsAcrossRegs);
maxDistBetCellsAcrossRegs = max(distBetCellsAcrossRegs);
rangeDistBetCellsAcrossRegs = [floor(minDistBetCellsAcrossRegs):10:maxDistBetCellsAcrossRegs+10];
binnedDistBetCellsAcrossRegs = discretize(distBetCellsAcrossRegs, rangeDistBetCellsAcrossRegs);
uniqueDistBinsAcrossRegs = unique(binnedDistBetCellsAcrossRegs);

mean_noLightActivitycorrAcrossRegs = []; sem_noLightActivitycorrAcrossRegs = [];
mean_midLightActivitycorrAcrossRegs = []; sem_midLightActivitycorrAcrossRegs = [];
mean_highLightActivitycorrAcrossRegs = []; sem_highLightActivitycorrAcrossRegs = [];
distCorrAcrossRegs = [];

for ii = 1:length(uniqueDistBinsAcrossRegs)
    mean_noLightActivitycorrAcrossRegs = [mean_noLightActivitycorrAcrossRegs, ...
        mean(noLightActivitycorrAcrossRegs(binnedDistBetCellsAcrossRegs==ii))];
    sem_noLightActivitycorrAcrossRegs = [sem_noLightActivitycorrAcrossRegs, ...
        std(noLightActivitycorrAcrossRegs(binnedDistBetCellsAcrossRegs==ii))/sqrt(sum(binnedDistBetCellsAcrossRegs==ii))];
    
    mean_midLightActivitycorrAcrossRegs = [mean_midLightActivitycorrAcrossRegs, ...
        mean(midLightActivitycorrAcrossRegs(binnedDistBetCellsAcrossRegs==ii))];
    sem_midLightActivitycorrAcrossRegs = [sem_midLightActivitycorrAcrossRegs, ...
        std(midLightActivitycorrAcrossRegs(binnedDistBetCellsAcrossRegs==ii))/sqrt(sum(binnedDistBetCellsAcrossRegs==ii))];
    
    mean_highLightActivitycorrAcrossRegs = [mean_highLightActivitycorrAcrossRegs, ...
        mean(highLightActivitycorrAcrossRegs(binnedDistBetCellsAcrossRegs==ii))];
    sem_highLightActivitycorrAcrossRegs = [sem_highLightActivitycorrAcrossRegs, ...
        std(highLightActivitycorrAcrossRegs(binnedDistBetCellsAcrossRegs==ii))/sqrt(sum(binnedDistBetCellsAcrossRegs==ii))];
    
    distCorrAcrossRegs = [distCorrAcrossRegs, uniqueDistBinsAcrossRegs(ii)];
end

% plot(rangeDistBetCellsAcrossRegs(1:23), mean_noLightActivitycorrAcrossRegs(1:23), 'color',[0.6,0.6,1]); 
% plot(rangeDistBetCellsAcrossRegs(1:23), mean_noLightActivitycorrAcrossRegs(1:23) + sem_noLightActivitycorrAcrossRegs(1:23), 'color',[0.6,0.6,1], 'LineStyle', '--');
% plot(rangeDistBetCellsAcrossRegs(1:23), mean_noLightActivitycorrAcrossRegs(1:23) - sem_noLightActivitycorrAcrossRegs(1:23), 'color',[0.6,0.6,1], 'LineStyle', '--');
% 
% plot(rangeDistBetCellsAcrossRegs(1:23), mean_midLightActivitycorrAcrossRegs(1:23), 'color',[1,0.6,0.6]); 
% plot(rangeDistBetCellsAcrossRegs(1:23), mean_midLightActivitycorrAcrossRegs(1:23) + sem_midLightActivitycorrAcrossRegs(1:23), 'color',[1,0.6,0.6], 'LineStyle', '--');
% plot(rangeDistBetCellsAcrossRegs(1:23), mean_midLightActivitycorrAcrossRegs(1:23) - sem_midLightActivitycorrAcrossRegs(1:23), 'color',[1,0.6,0.6], 'LineStyle', '--');

plot(rangeDistBetCellsAcrossRegs(1:23), mean_highLightActivitycorrAcrossRegs(1:23), 'color',[0.6,0.6,0.6]); 
plot(rangeDistBetCellsAcrossRegs(1:23), mean_highLightActivitycorrAcrossRegs(1:23) + sem_highLightActivitycorrAcrossRegs(1:23), 'color',[0.6,0.6,0.6], 'LineStyle', '--');
plot(rangeDistBetCellsAcrossRegs(1:23), mean_highLightActivitycorrAcrossRegs(1:23) - sem_highLightActivitycorrAcrossRegs(1:23), 'color',[0.6,0.6,0.6], 'LineStyle', '--');

xlabel('Distance(um)'); ylabel('Correlation coefficients');
xlim([0,170]);
set(gca,'FontSize',20);
set(gcf, 'PaperUnits', 'inches');
x_width=10 ;y_width=8;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
savefig(figbox, 'Mean,sem comparing xcorr(tau=0) between cells across redcell neighbourhoods high lights');
saveas(figbox, 'Mean,sem comparing xcorr(tau=0) between cells across redcell neighbourhoods high lights','png');

%% average activity change in the presence of highest sound and light conditions

fig12 = figure(); ax12=axes(fig12); hold on;
mtresponsive =0;
for ii = 1:length(responsiveRedcells)
    closestCells = cellIdxs(closestRedcell(:,2)==ii);
    mt=0;
    for jj = 1:length(closestCells)
        mt = mt+mean_Trace(3,:,closestCells(jj));
    end
    mtresponsive = mtresponsive+mt/length(closestCells); 
end
plot(ax12, mtresponsive/length(responsiveRedcells), 'r');

mtresponsive =0;
for ii = 1:length(responsiveRedcells)
    closestCells = cellIdxs(closestRedcell(:,2)==ii);
    mt=0;
    for jj = 1:length(closestCells)
        mt = mt+mean_Trace(19,:,closestCells(jj));
    end
    mtresponsive = mtresponsive+mt/length(closestCells); 
end
plot(ax12, mtresponsive/length(responsiveRedcells), 'b');

mtresponsive =0;
for ii = 1:length(responsiveRedcells)
    closestCells = cellIdxs(closestRedcell(:,2)==ii);
    mt=0;
    for jj = 1:length(closestCells)
        mt = mt+mean_Trace(21,:,closestCells(jj));
    end
    mtresponsive = mtresponsive+mt/length(closestCells); 
end
plot(ax12, mtresponsive/length(responsiveRedcells), 'k');

xlabel('Time'); ylabel('Average activity'); 
savefig(fig12, 'Average activity given 90db sound and highest light');
saveas(fig12, 'Average activity given 90db sound and highest light','png');

%% variation of correlation with distance of cells that lie on the perpendicular bisector
% The last time I checked this (02-05), it did not give anything
% insightful. Maybe look at response functions for these cells next.

close all;

perpBisectorCells = cell(length(responsiveRedcells),length(responsiveRedcells));
for ii=1:size(raster,3)
    [disttoRedCells,idxRedCells] = sort(DistfromRedcells(ii,:)); 
    closestTwoRedCells = idxRedCells(1:2);
    ratioDisttoRedCells = disttoRedCells(2)/disttoRedCells(1);
    if ratioDisttoRedCells < 1.2
        perpBisectorCells{closestTwoRedCells(1), closestTwoRedCells(2)} = ...
            [ii,perpBisectorCells{closestTwoRedCells(1), closestTwoRedCells(2)}];
    end 
end

openfig('Cell layout.fig'); hold on;
for ii = 1:length(perpBisectorCells)
    for jj = 1:length(perpBisectorCells)
        if ~isempty(perpBisectorCells{ii,jj})
            for kk = 1:length(perpBisectorCells{ii,jj})
                h = patch(spatialInfo.ROIs{perpBisectorCells{ii,jj}(kk)}(:,2),...
                    spatialInfo.ROIs{perpBisectorCells{ii,jj}(kk)}(:,1),'y','facealpha',0.8,'edgealpha',0,...
                    'LineWidth',1,'facecolor','k');
            end
        end
    end
end

corrCellsPerpBisector = cell(size(perpBisectorCells));
for ii = 1:length(perpBisectorCells)
    for jj = 1:length(perpBisectorCells)
        if length(perpBisectorCells{ii,jj})>1
            corrCellsPerpBisector{ii,jj} = zeros(length(perpBisectorCells{ii,jj})*(length(perpBisectorCells{ii,jj})-1),22);
            for kk = 1:length(perpBisectorCells{ii,jj})
                othercloseCells = perpBisectorCells{ii,jj};
                othercloseCells(kk) = [];
                corrCellsPerpBisector{ii,jj}(length(othercloseCells)*(kk-1)+1:length(othercloseCells)*kk,1) ...
                    = vecnorm(Allcells_ROI(othercloseCells,:)-Allcells_ROI(perpBisectorCells{ii,jj}(kk),:),2,2);
                for ll = 1:length(othercloseCells)
                    for cc = 1:21
                        coeff = corrcoef(mean_Trace(cc,:,othercloseCells(ll)),mean_Trace(cc,:,perpBisectorCells{ii,jj}(kk))); 
                        corrCellsPerpBisector{ii,jj}(length(othercloseCells)*(kk-1)+ll,cc+1) = coeff(1,2);
                    end
                end
            end
        end
    end
end

fig13 = figure(); ax13=axes(fig13);
noLightcorrPerp = []; distBetPerpCells = [];
for ii = 1:length(perpBisectorCells)
    for jj = 1:length(perpBisectorCells)
        if length(perpBisectorCells{ii,jj})>1
            [distclosestCells, sortDistIdx] = sort(corrCellsPerpBisector{ii,jj}(:,1));
            for i=[1,4,7,10,13,16,19]
                plot(ax13, distclosestCells, corrCellsPerpBisector{ii,jj}(sortDistIdx,i+1),'k.');
                noLightcorrPerp = [noLightcorrPerp; corrCellsPerpBisector{ii,jj}(sortDistIdx,i+1)];
                distBetPerpCells = [distBetPerpCells; distclosestCells];
                hold on;
            end       
        end
    end
end
xlabel('Distance'); ylabel('Correlation coefficients');
ylim([-1,1]);

fig14 = figure(); ax14=axes(fig14);
midLightcorrPerp = []; 
for ii = 1:length(perpBisectorCells)
    for jj = 1:length(perpBisectorCells)
        if length(perpBisectorCells{ii,jj})>1
            [distclosestCells, sortDistIdx] = sort(corrCellsPerpBisector{ii,jj}(:,1));
            for i=[2,5,8,11,14,17,20]
                plot(ax14, distclosestCells, corrCellsPerpBisector{ii,jj}(sortDistIdx,i+1),'k.');
                midLightcorrPerp = [midLightcorrPerp; corrCellsPerpBisector{ii,jj}(sortDistIdx,i+1)];
                hold on;
            end       
        end
    end
end
xlabel('Distance'); ylabel('Correlation coefficients');
ylim([-1,1]);

fig15 = figure(); ax15=axes(fig15);
highLightcorrPerp = []; 
for ii = 1:length(perpBisectorCells)
    for jj = 1:length(perpBisectorCells)
        if length(perpBisectorCells{ii,jj})>1
            [distclosestCells, sortDistIdx] = sort(corrCellsPerpBisector{ii,jj}(:,1));
            for i=[3,6,9,12,15,18,21]
                plot(ax15, distclosestCells, corrCellsPerpBisector{ii,jj}(sortDistIdx,i+1),'k.');
                highLightcorrPerp = [highLightcorrPerp; corrCellsPerpBisector{ii,jj}(sortDistIdx,i+1)];                
                hold on;
            end       
        end
    end
end
xlabel('Distance'); ylabel('Correlation coefficients');
ylim([-1,1]);

minDistBetPerpCells = min(distBetPerpCells);
maxDistBetPerpCells = max(distBetPerpCells);
rangeDistBetPerpCells = [minDistBetCells:10:maxDistBetCells+10];
binnedDistBetPerpCells = discretize(distBetPerpCells,rangeDistBetPerpCells);

figbox = figure(); axbox=axes(figbox);
for ii = 1:length(unique(binnedDistBetPerpCells))
    plot(rangeDistBetPerpCells(ii),mean(noLightcorrPerp(binnedDistBetPerpCells==ii)),'r.')
    errorbar(rangeDistBetPerpCells(ii), mean(noLightcorrPerp(binnedDistBetPerpCells==ii)), ...
        std(noLightcorrPerp(binnedDistBetPerpCells==ii)), std(noLightcorrPerp(binnedDistBetPerpCells==ii)),...
    'Color', 'r');
    hold on;
end

for ii = 1:length(unique(binnedDistBetPerpCells))
    plot(rangeDistBetPerpCells(ii),mean(midLightcorrPerp(binnedDistBetPerpCells==ii)),'g.')
    errorbar(rangeDistBetPerpCells(ii), mean(midLightcorrPerp(binnedDistBetPerpCells==ii)), ...
        std(midLightcorrPerp(binnedDistBetPerpCells==ii)), std(midLightcorrPerp(binnedDistBetPerpCells==ii)),...
    'Color', 'g');
    hold on;
end

for ii = 1:length(unique(binnedDistBetPerpCells))
    plot(rangeDistBetPerpCells(ii),mean(highLightcorrPerp(binnedDistBetPerpCells==ii)),'k.')
    errorbar(rangeDistBetPerpCells(ii), mean(highLightcorrPerp(binnedDistBetPerpCells==ii)), ...
        std(highLightcorrPerp(binnedDistBetPerpCells==ii)), std(highLightcorrPerp(binnedDistBetPerpCells==ii)),...
    'Color', 'k');
    hold on;
end

xlabel('Distance'); ylabel('Correlation coefficients');

%% original code for influence and correlation functions

av_activity = zeros(size(av_raster_poststim));
for jj=1:size(raster,3)
    for ii=1:length(stimInfo.index)  
        av_activity(ii,jj) = nanmean(squeeze(nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,:,jj),2)));
    end
end

ActivitysurroundingCells = cell(length(closestCells),1);
midResponse = cell(length(closestCells),1);
highResponse = cell(length(closestCells),1);
signalCorrelations = cell(length(closestCells),1);
for i = 1:length(ActivitysurroundingCells)
    ActivitysurroundingCells{i} = av_raster_poststim([1,4,7,10,13,16,19],closestCells{i});
    midResponse{i} = av_raster_poststim([2,5,8,11,14,17,20],closestCells{i});
    highResponse{i} = av_raster_poststim([3,6,9,12,15,18,21],closestCells{i});
    signalCorrelations{i} = corrcoef(av_raster_poststim([1,4,7,10,13,16,19],closestCells{i}));
end

signalCorrelationsRedcell = [];
ActivitygivenRedCell = [];
midResponseRedcell= [];
highResponseRedcell = [];

for i = 1:length(closestCells)
    if length(closestCells{i})>1
       logicalRedCell = closestCells{i}==Ind_cells_red(i);
       logicalOtherCells = ~logicalRedCell;
       ActivitygivenRedCell = [ActivitygivenRedCell, ActivitysurroundingCells{i}(:,logicalOtherCells)];
       highResponseRedcell = [highResponseRedcell, highResponse{i}(:,logicalOtherCells)];
       midResponseRedcell = [midResponseRedcell, midResponse{i}(:,logicalOtherCells)];
       signalCorrelationsRedcell = [signalCorrelationsRedcell, ...
           signalCorrelations{i}(logicalRedCell, logicalOtherCells)];
    end
end

highResponseRedcell = (highResponseRedcell-ActivitygivenRedCell)/0.7;
midResponseRedcell = (midResponseRedcell-ActivitygivenRedCell)/0.304;

fsig = figure();
histogram(signalCorrelationsRedcell);
title('MT064 20190822, signal correlations');
savefig(fsig, 'MT064_20190822_signalCorrelations_av_raster_poststim.fig')

fhigh = figure();
for i=1:7
    plot(signalCorrelationsRedcell,highResponseRedcell(i,:),'.');
    hold on;
end
[out_sig, idx_sig] = sort(signalCorrelationsRedcell);
out_highresp = highResponseRedcell(:,idx_sig);
plot(out_sig,mean(out_highresp,1),'k');
title('MT064 20190822, high response');
savefig(fhigh, 'MT064_20190822_highResponse_av_raster_poststim.fig')

fmid = figure();
for i=1:7
    plot(signalCorrelationsRedcell,midResponseRedcell(i,:),'.');
    hold on;
end
[out_sig, idx_sig] = sort(signalCorrelationsRedcell);
out_midresp = midResponseRedcell(:,idx_sig);
plot(out_sig,mean(out_midresp,1),'k');
title('MT064 20190822, mid response');
savefig(fmid,'MT064_20190822_midResponse_av_raster_poststim.fig')

%% original code for closest cells
DistbetRedcells = zeros(length(responsiveRedcells),length(responsiveRedcells));
for i = 1:size(DistbetRedcells,1)
    DistbetRedcells(i,:) = vecnorm(Redcells_ROI-Redcells_ROI(i,:),2,2);
end

sortDist = sort(DistbetRedcells(DistbetRedcells>0));
RadialDist = sortDist(1)/2;

Allcells_ROI = zeros(length(spatialInfo.ROIs),2);
for i = 1:size(Allcells_ROI,1)
    Allcells_ROI(i,:) = mean(spatialInfo.ROIs{i});
end

closestCells = cell(length(responsiveRedcells),1);
useMinDist = 1;
for i = 1:size(Redcells_ROI,1)
    if sum(vecnorm(Allcells_ROI-Redcells_ROI(i,:),2,2)<sortDist(1)*3/4)<2
        useMinDist = 0;
    end
end

for i = 1:size(Redcells_ROI,1)
    if useMinDist
        closestCells{i} = find(vecnorm(Allcells_ROI-Redcells_ROI(i,:),2,2)<sortDist(1)*3/4);
    else
        closestCells{i} = find(vecnorm(Allcells_ROI-Redcells_ROI(i,:),2,2)<sortDist(3)*3/4);
    end
end
%}