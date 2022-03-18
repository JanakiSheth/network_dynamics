%% Goal of the script

% To create activity datasets found in the nnmf folder on taco. 
% "activity data" is the original recordings, with some processing: remove laser bands, create rasters, figure out sound+ /sound- / laser+ cells etc 

%% Noise click + optogenetics Analysis : 

% OUTLINE
% V load the data
% V make raster triggered on events
% V do stats on each condition to see what cells respond, and from that
% get which cells are sound-responsive (from max sound no opto condition)
% or laser-responsive
% - plot curve for sound-responsive cells of response versus sound
% amplitude for different laser levels (look at individual cells and average)
% - from that 1) fit slope to get gain 2) measure max amplitude of response

% - use 99% criterion, and compare using [0 1] or [1 2] seconds after sound
% onset 

%%  Load data of interest - under /data/ folder

clear variables
close all

% Select the files to open
cd '/home/janaki/Dropbox/project_with_melanie/DataForJanaki/ImagingDatasetsAndCorrespondingInfo_DirectlyFromExpts/'
[Filename,Pathname] = uigetfile('*.mat','Select the file to open');
load(sprintf('%s%s%s',Pathname,'/',Filename)); 

mouse=Filename(3:7);
date=str2double(Filename(9:16));
mouseline = 'SOM';

% look if there's a RedCells.mat file
RedCellfiles = dir([mouse,'_',num2str(date),'_','RedCells*.mat']);
load(RedCellfiles.name);

% load FRA file
FRAfiles = dir([mouse,'_',num2str(date),'_','FRA*.mat']);
load(FRAfiles.name);

%% Calculate when laser bands are, then remove them with registering of concerned cells

 % I will modify calcium directly, so make a copy in calcium_raw before
 % modifications
calcium_raw=calcium;

fr=exptInfo.fr; % frame rate

% find in original recording where the laser bands are

fs=stimInfo.fs;

index_laserevents=find(stimInfo.index(:,2)~=0);

% find laser band indices in terms of fs
indfs_laser=[];
for i=1:length(stimInfo.order)
    if ismember(stimInfo.order(i),index_laserevents)
        begin_indfs=exptInfo.preStimSilence*fs+(i-1)*fs*(stimInfo.ISI/1000+stimInfo.tDur_opto);
        for j=1:stimInfo.nPulses % number of laser pulses  
            interval=stimInfo.evOn_opto(j):stimInfo.evOff_opto(j);
            indfs_laser=[indfs_laser begin_indfs+interval];
        end
        
    end    
end
   
% convert the indices in terms of frame lines
indlinefr_laser=indfs_laser*(fr*dat.ops.Ly)/fs;
indlinefr_laser=unique(round(indlinefr_laser));   

% indfr_laser=unique(ceil(indlinefr_laser/dat.ops.Ly)); % ceil because really what we are looking for is mod

% go to the frame and lines where there is laser, 

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
            calcium.npilSubTraces(j,indfr_laser(i))=NaN;
            indcell_reg(j,indfr_laser(i))=1;
        end
    end
end

%% Extract Parameters

fr=exptInfo.fr; % frame rate
relevant_info=['t_{sound} = ',num2str(stimInfo.nClicks/stimInfo.clickRate),'s ; ','t_{opto} = ', num2str(stimInfo.tDur_opto),'s ; '];
relevant_info2 = ['max Amp_{sound} = ', num2str(max(stimInfo.intensity)),'dB; ','max Amp_{opto} = ', num2str(max(stimInfo.amplitude_opto*100)),'%'];

% make the raster, eg cut up time series in each trial

preOnsetTime=1; %1 second
postOnsetTime=round((events.eventsOn(2)-events.eventsOn(1))/fr)-1; % duration between two stimuli
doDFoverF0=1;
raster = makeCaRaster_NaN_MT(calcium.npilSubTraces,events.eventsOn,round(preOnsetTime*fr),round(postOnsetTime*fr),doDFoverF0);

% average over all events and plot the raster cell versus time

av_raster = squeeze(nanmean(raster,1));

% sort depending on components of index

repeats=stimInfo.repeats;
[raster_ordered, mean_Trace, std_mean_Trace]=makeOrdRasterAndAveTrace_MT(raster,stimInfo.index,stimInfo.order,stimInfo.repeats);

% number of components along each dimension of index
for i=1:size(stimInfo.index,2)
    dim_index(i)=length(unique(stimInfo.index(:,i)));
end

dim_relevant=find(dim_index~=1);

% plot all rasters for different conditions when only one varying dimension ordered in a  matrix, on several figures if necessary    

fig(1)=figure(1);
hold on
for ii=1:length(stimInfo.index)
    subplot(dim_index(dim_relevant(1)),dim_index(dim_relevant(2)),ii)
    imagesc([1:size(raster_ordered,2)]/fr,[1:size(raster_ordered,3)],squeeze(mean_Trace(ii,:,:))')
    colormap gray
    set(gca,'TickDir','out','FontSize',14)
    xlabel('Time (s)')
    ylabel('Cells')
    title({['index 1 = ',num2str(stimInfo.index(ii,dim_relevant(1)))],['index 2 = ',num2str(stimInfo.index(ii,dim_relevant(2)))]})
    Ind_cells_red=find(RedCells.isCell_red);
    for i=1:length(Ind_cells_red)
        line([0 0.5],[Ind_cells_red(i) Ind_cells_red(i)],'Color','r','LineWidth',2)
    end
end

% figure with just the max parameter values

for i=1:length(dim_relevant)
    % find indices where dim_relevant is max
    b=find(stimInfo.index(:,dim_relevant(i))==max(stimInfo.index(:,dim_relevant(i))));
    % find indices where dim non relevant is zero
    if dim_relevant==1; % condition added for only laser changing and sound always on
        ind_max_dim_relevant(i)=b(find(stimInfo.index(b)));
    else
        ind_max_dim_relevant(i)=b(find(stimInfo.index(b,setdiff(dim_relevant,dim_relevant(i)))==0));
    end
end


if dim_relevant > 1
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
        Ind_cells_red=find(RedCells.isCell_red);
        for i=1:length(Ind_cells_red)
            line([0 0.5],[Ind_cells_red(i) Ind_cells_red(i)],'Color','r','LineWidth',2)
        end
    end
end


%% sort for every condition depending on sound response and red intensity and show rasters

Tdelay_windur=1; % window duration of 1s for averaging.
Tdelay_test=(0:0.03:postOnsetTime-Tdelay_windur);
alpha=0.01;

Tdelay=[];
tstats_time=zeros(length(stimInfo.index),size(raster,3),length(Tdelay_test));
pstats_time=zeros(length(stimInfo.index),size(raster,3),length(Tdelay_test));
zscore=zeros(length(stimInfo.index),size(raster,3),length(Tdelay_test));
av_raster_poststim_time=zeros(length(stimInfo.index),size(raster,3),length(Tdelay_test));
sem_raster_poststim_time=zeros(length(stimInfo.index),size(raster,3),length(Tdelay_test));

% calculate stats for all cells, all stims, all windows
for jj=1:size(raster,3)
    for ii=1:length(stimInfo.index)  
        % calculate stats with all Tdelay_test
        for kk=1:length(Tdelay_test)
            [tstats_time(ii,jj,kk),pstats_time(ii,jj,kk)]=ttest(nanmean(raster_ordered((ii-1)*10+1:10*(ii-1)+repeats,1:round(fr*preOnsetTime),jj),2),...
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

Nsignifcells_allstims_time = zeros(length(Tdelay_test),1);
for kk=1:length(Tdelay_test)
    Nsignifcells_allstims_time(kk)=length(find(sum(tstats_time(:,:,kk))>=1));
end
Tdelay_allstims=Tdelay_test(find(Nsignifcells_allstims_time(:)==max(Nsignifcells_allstims_time(:)),1,'first'));

cm=[0 0 0;
    1 0.7 0.7;
    1 0 0];

set(groot,'defaultAxesColorOrder',cm)

figure, 
% plot(Tdelay_test,Nsignifcells_time,'LineWidth',1)
hold on, 
plot(Tdelay_test,Nsignifcells_allstims_time,'b','LineWidth',2)
plot(Tdelay_allstims,max(Nsignifcells_allstims_time),'bo','LineWidth',2,'MarkerSize',12)
xlabel('Onset of response window (s)')
ylabel('Number of significantly responsive cells')
title({[mouse,' - ',num2str(date),' - ',num2str(size(raster,3)),' cells'],...
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
    [a, Ind_sort_poststim{i}]=sort(av_raster_poststim(ind_max_dim_relevant(i),:),'descend');
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


%% Plot curves of cell response vs sound amplitude for cell responsive-cells
% Figure with overlap individual sound-responsive cell curves, different
% laser for different panels
% figure with average for the different laser amplitudes overlapped
% modified 2019/09/26 : make cells be only Sound+ or only Sound-

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
        % ['Red cells = red : n = ',num2str(length(Ind_CellSelection{2})), '/', num2str(size(raster,3))] }

%% save workspace

close all
if mouseline == 'SOM'
    path_workspace='/home/janaki/Dropbox/project_with_melanie/DataForJanaki/networkFiles/SOM_ActivityData/';
elseif mouseline == 'VIP'
    path_workspace='/home/janaki/Dropbox/project_with_melanie/DataForJanaki/networkFiles/VIP_ActivityData/';
end
cd(path_workspace)
name_workspace=['workspace_',mouse,'_',num2str(date),'_compil',regexprep(datestr(today,'yyyy-mm-dd'),'[-]','')];
save (name_workspace)

     %% stats for change sound / sound+laser

     % if different sound levels, take the same sound level wihtout the
     % laser as the baseline
     
for ii=1:length(stimInfo.index) 
    b=find(stimInfo.index(:,1)==stimInfo.index(ii,1) & stimInfo.index(:,2)==0);
soundrepeats_baseline=[repeats*(b-1)+1:repeats*b];
% just for that SOM recording where I combined data together
% soundrepeats_baseline=[sum(stimInfo.repeats(1:b-1))+1:sum(stimInfo.repeats(1:b))];
% repeats=stimInfo.repeats(ii);
for jj=1:size(raster,3)
[tstats_SoundAsBaseline(ii,jj) pstats_SoundAsBaseline(ii,jj)]=ttest2(mean(raster_ordered(soundrepeats_baseline,round(fr*(preOnsetTime+Tdelay(jj))):round(fr*(preOnsetTime+Tdelay(jj)+Tdelay_windur)),jj),2),mean(raster_ordered(repeats*(ii-1)+1:repeats*ii,round(fr*(preOnsetTime+Tdelay(jj))):round(fr*(preOnsetTime+Tdelay(jj)+Tdelay_windur)),jj),2),'Alpha',alpha);
% [tstats_SoundAsBaseline(ii,jj) pstats_SoundAsBaseline(ii,jj)]=ttest2(mean(raster_ordered(soundrepeats_baseline,round(fr*(preOnsetTime+Tdelay(1))):round(fr*(preOnsetTime+Tdelay(2))),jj),2),mean(raster_ordered(sum(stimInfo.repeats(1:ii-1))+1:sum(stimInfo.repeats(1:ii)),round(fr*(preOnsetTime+Tdelay(1))):round(fr*(preOnsetTime+Tdelay(2))),jj),2),'Alpha',alpha);
end
end

% mark which cells show an activity change, not ordered
Cells_ActivityChange_SoundAsBaseline=cell(length(stimInfo.index),1);
for kk=1:length(stimInfo.index)
    Cells_ActivityChange_SoundAsBaseline{kk}=zeros(1,size(raster,3));
    b=find(stimInfo.index(:,1)==stimInfo.index(kk,1) & stimInfo.index(:,2)==0);
    for jj=1:size(raster,3)
        if pstats_SoundAsBaseline(kk,jj)<alpha && av_raster_poststim(kk,jj)-av_raster_poststim(b,jj) >0
            Cells_ActivityChange_SoundAsBaseline{kk}(jj)=1;
        elseif pstats_SoundAsBaseline(kk,jj)<alpha && av_raster_poststim(kk,jj)-av_raster_poststim(b,jj) <0
            Cells_ActivityChange_SoundAsBaseline{kk}(jj)=-1;
        end
    end
end
     
criterion_SoundAsBaseline=cell(4,1);
Ind_CellSelection_SoundAsBaseline=cell(4,1);

%%%% 1) All cells %%%%
    Ind_CellSelection_SoundAsBaseline{1}=[1:size(av_raster,2)];
    criterion_SoundAsBaseline{1} = 'All cells';
    criterion_SoundAsBaseline_short{1} = 'All';
% %%%% 2) All red cells %%%%
if exist('RedCells')==1
    Ind_CellSelection_SoundAsBaseline{2}=Ind_cells_red;
    criterion_SoundAsBaseline{2} = 'All red cells';
    criterion_SoundAsBaseline_short{2} = 'Red';
end
%%%% 3) All cells increasing their response from sound to sound+laser %%%%
    b=find(stimInfo.index(:,1)~=0);
    Ind_CellSelection_SoundAsBaseline{3}=[];
    for i=2:length(b)
        Ind_CellSelection_SoundAsBaseline{3}=union(Ind_CellSelection_SoundAsBaseline{3},find(Cells_ActivityChange_SoundAsBaseline{b(i)}==1));
    end
    criterion_SoundAsBaseline{3} = 'All cells increasing their response from sound to sound+laser';
    criterion_SoundAsBaseline_short{3} = 'SoundLaser+';
%%%% 4) All cells their response from sound to sound+laser %%%%
    b=find(stimInfo.index(:,1)~=0);
    Ind_CellSelection_SoundAsBaseline{4}=[];
    for i=2:length(b)
        Ind_CellSelection_SoundAsBaseline{4}=union(Ind_CellSelection_SoundAsBaseline{4},find(Cells_ActivityChange_SoundAsBaseline{b(i)}==-1));
    end
    criterion_SoundAsBaseline{4} = 'All cells decreasing their response from sound to sound+laser';
    criterion_SoundAsBaseline_short{4} = 'SoundLaser-';

    
    
         %%%%%%%%%%% stats for change laser / sound+laser

     % if different sound levels, take the same sound level wihtout the
     % laser as the baseline
     
for ii=1:length(stimInfo.index) 
    b=find(stimInfo.index(:,2)==stimInfo.index(ii,2) & stimInfo.index(:,1)==0);
laserrepeats_baseline=[repeats*(b-1)+1:repeats*b];
for jj=1:size(raster,3)
[tstats_LaserAsBaseline(ii,jj) pstats_LaserAsBaseline(ii,jj)]=ttest2(mean(raster_ordered(laserrepeats_baseline,round(fr*(preOnsetTime+Tdelay(jj))):round(fr*(preOnsetTime+Tdelay(jj))+Tdelay_windur),jj),2),mean(raster_ordered(repeats*(ii-1)+1:repeats*ii,round(fr*(preOnsetTime+Tdelay(jj))):round(fr*(preOnsetTime+Tdelay(jj)+Tdelay_windur)),jj),2),'Alpha',alpha);
end
end

% mark which cells show an activity change, not ordered
Cells_ActivityChange_LaserAsBaseline=cell(length(stimInfo.index),1);
for kk=1:length(stimInfo.index)
    Cells_ActivityChange_LaserAsBaseline{kk}=zeros(1,size(raster,3));
    b=find(stimInfo.index(:,2)==stimInfo.index(ii,2) & stimInfo.index(:,1)==0);
    for jj=1:size(raster,3)
        if pstats_LaserAsBaseline(kk,jj)<alpha && av_raster_poststim(kk,jj)-av_raster_poststim(b,jj) >0
            Cells_ActivityChange_LaserAsBaseline{kk}(jj)=1;
        elseif pstats_LaserAsBaseline(kk,jj)<alpha && av_raster_poststim(kk,jj)-av_raster_poststim(b,jj) <0
            Cells_ActivityChange_LaserAsBaseline{kk}(jj)=-1;
        end
    end
end
     
criterion_LaserAsBaseline=cell(4,1);
Ind_CellSelection_LaserAsBaseline=cell(4,1);

%%%% 1) All cells %%%%
    Ind_CellSelection_LaserAsBaseline{1}=[1:size(av_raster,2)];
    criterion_LaserAsBaseline{1} = 'All cells';
% %%%% 2) All red cells %%%%
if exist('RedCells')==1
    Ind_CellSelection_LaserAsBaseline{2}=Ind_cells_red;
    criterion_LaserAsBaseline{2} = 'All red cells';
end
%%%% 3) All cells increasing their response from laser to sound+laser %%%%
    b=find(stimInfo.index(:,2)~=0);
    Ind_CellSelection_LaserAsBaseline{3}=[];
    for i=2:length(b)
        Ind_CellSelection_LaserAsBaseline{3}=union(Ind_CellSelection_LaserAsBaseline{3},find(Cells_ActivityChange_LaserAsBaseline{b(i)}==1));
    end
    criterion_LaserAsBaseline{3} = 'All cells increasing their response from laser to sound+laser';
%%%% 4) All cells their response from laser to sound+laser %%%%
    b=find(stimInfo.index(:,2)~=0);
    Ind_CellSelection_LaserAsBaseline{4}=[];
    for i=2:length(b)
        Ind_CellSelection_LaserAsBaseline{4}=union(Ind_CellSelection_LaserAsBaseline{4},find(Cells_ActivityChange_LaserAsBaseline{b(i)}==-1));
    end
    criterion_LaserAsBaseline{4} = 'All cells decreasing their response from laser to sound+laser';

%% save cell groups for Venn diagrams

name_fig=[num2str(date),'_',mouse,'_',num2str(Filename(16:end-4))];
cd 'C:\Users\Melanie\Dropbox\Postdoc\Analysis\'

CellSelection=struct;
CellSelection.Ind_CellSelection=Ind_CellSelection;
CellSelection.Ind_CellSelection_SoundAsBaseline=Ind_CellSelection_SoundAsBaseline;
CellSelection.Ind_CellSelection_LaserAsBaseline=Ind_CellSelection_LaserAsBaseline;
CellSelection.criterion=criterion;
CellSelection.criterion_SoundAsBaseline=criterion_SoundAsBaseline;
CellSelection.criterion_LaserAsBaseline=criterion_LaserAsBaseline;

save([name_fig,'_CellSelection'],'CellSelection')

%% concatenate all cell selection files for different recordings in one

CellSelectionAllVIP=struct;
CellSelectionAllVIP.FilesOrder={'20190117_MT032_08';'20190122_MT032_03';'20190130_MT032_03';'20190121_MT037_03'};
CellSelectionAllVIP.ncell_PerFile=[];
CellSelectionAllVIP.Ind_CellSelection=cell(6,1);
CellSelectionAllVIP.Ind_CellSelection_SoundAsBaseline=cell(4,1);
CellSelectionAllVIP.Ind_CellSelection_LaserAsBaseline=cell(4,1);

%manually load all the CellSelection files in order defined above, then
%compile the lines below the load line
load('20190117_MT032_FRA_OptoStim_08_CellSelection.mat')


ncellall=sum(CellSelectionAllVIP.ncell_PerFile);
CellSelectionAllVIP.ncell_PerFile=[CellSelectionAllVIP.ncell_PerFile; length(CellSelection.Ind_CellSelection{1})];

for i=1:6
CellSelectionAllVIP.Ind_CellSelection{i}=[CellSelectionAllVIP.Ind_CellSelection{i}; reshape(CellSelection.Ind_CellSelection{i},[],1)+ncellall];
end
for i=1:4
CellSelectionAllVIP.Ind_CellSelection_SoundAsBaseline{i}=[CellSelectionAllVIP.Ind_CellSelection_SoundAsBaseline{i}; reshape(CellSelection.Ind_CellSelection_SoundAsBaseline{i},[],1)+ncellall];
CellSelectionAllVIP.Ind_CellSelection_LaserAsBaseline{i}=[CellSelectionAllVIP.Ind_CellSelection_LaserAsBaseline{i}; reshape(CellSelection.Ind_CellSelection_LaserAsBaseline{i},[],1)+ncellall];
end
CellSelectionAllVIP.criterion=CellSelection.criterion;
CellSelectionAllVIP.criterion_SoundAsBaseline=CellSelection.criterion_SoundAsBaseline;
CellSelectionAllVIP.criterion_LaserAsBaseline=CellSelection.criterion_LaserAsBaseline;


%% distance between cells

% build matrix with distance between every cell (symmetric or just
% triangluar for calculation and then fill in other half so that index doesn't matter)
 % for cell position, taking the middle of the max extension
 micronsPerPixel=exptInfo.recInfo.micronsPerPixel(1); 
 Dist_betweenCells=zeros(size(spatialInfo.ROIs,2),size(spatialInfo.ROIs,2));
 
 for ii=1:size(spatialInfo.ROIs,2)
     for jj=ii:size(spatialInfo.ROIs,2)
         Centerix=(max(spatialInfo.ROIs{ii}(:,1))+min(spatialInfo.ROIs{ii}(:,1)))/2* micronsPerPixel;
         Centeriy=(max(spatialInfo.ROIs{ii}(:,2))+min(spatialInfo.ROIs{ii}(:,2)))/2* micronsPerPixel;
         Centerjx=(max(spatialInfo.ROIs{jj}(:,1))+min(spatialInfo.ROIs{jj}(:,1)))/2* micronsPerPixel;
         Centerjy=(max(spatialInfo.ROIs{jj}(:,2))+min(spatialInfo.ROIs{jj}(:,2)))/2* micronsPerPixel;
         
        Dist_betweenCells(ii,jj)=sqrt((Centerix-Centerjx)^2+(Centeriy-Centerjy)^2);
        Dist_betweenCells(jj,ii)=Dist_betweenCells(ii,jj);
     end
 end

 %% plot histogram of distance cells of a same group
  Dist_intragroup=cell(6,1);
  
 figure
 for q=1:6
 subplot(3,2,q)
 b=[];
 for ii=1:length(Ind_CellSelection{q})-1 % red cells
     b=[b Dist_betweenCells(Ind_CellSelection{q}(ii),Ind_CellSelection{q}(ii+1:end))];
 end
 Dist_intragroup{q}=b;
     histogram(b,'BinLimits',[min(Dist_intragroup{1}) max(Dist_intragroup{1})],'BinWidth',20,'DisplayStyle','stairs','Normalization','probability')
     hold on
     title({criterion_short{q}, ['mean = ',num2str(mean(b),3),' µm; median = ',num2str(median(b),3),' µm'],['ndist=',num2str(length(b)),'; ncell=',num2str(length(Ind_CellSelection{q}))]})
     xlabel('Distance between cells (µm)')
     ylabel('Counts / number of inputs')
 end

 % next step is to concatenate all distances from the different recordings
 % and do the same figure with all the data
 
 
  %% plot histogram of distance cells between different groups
  
  % 6 groups, so distance between groups is 
  
  Dist_intergroup=cell(6,6);
  
 figure
 for q=2:6
     for q2=1:6
         if q2>q
 subplot(4,4,4*(q-2)+q2-2)
 b=[];
 for ii=1:length(Ind_CellSelection{q}) % red cells
     for jj=1:length(Ind_CellSelection{q2})
     b=[b Dist_betweenCells(Ind_CellSelection{q}(ii),Ind_CellSelection{q2}(jj))];
     end
 end
 Dist_intergroup{q,q2}=b;
     histogram(b)
histogram(b,'BinLimits',[min(Dist_intragroup{1}) max(Dist_intragroup{1})],'BinWidth',20,'DisplayStyle','bar','Normalization','probability')

     hold on
      title({[criterion_short{q},' vs ',criterion_short{q2}], ['mean = ',num2str(mean(b),3),' µm; median = ',num2str(median(b),3),' µm'],['ndist=',num2str(length(b)),';']})
     xlabel('Distance between cells (µm)')
     ylabel('Counts / number of inputs')
         end
     end
 end

 % next step is to concatenate all distances from the different recordings
 % and do the same figure with all the data
 
 % next step put the intragroup distances also on the diagonal
 % compare to what theoretical distribution would be if evenly spaced out
 % implement closest neighbour
 % see whether there is a better measure for spatial clusering
 
%% calculate first neighbour for each cell and see group belongings
 
% actually might be interesting to calculate the first ~5 neighbours and
% then plot them as a histogram of group belonging

% one array will be the identity of the first 5 or 10 neighbours fro every
% cell

% another array will be the count of nth neighbour depending on their
% identity

FirstNeighboursCount=cell(7,1);
for q=1:7
    FirstNeighboursCount{q}=zeros(length(Ind_CellSelection{1})-1,7);
end
    
 for i=1:length(Ind_CellSelection{1}) %all cells
  [SortedDistance(i,:) FirstNeighbours(i,:)]=sort(Dist_betweenCells(i,:));
  for q=1:7; %make 7 for all the cells that do not respond
      if ismember(i,Ind_CellSelection{q});
          for j=2:length(Ind_CellSelection{1}); %1st cell will be identity, jth closest neighbour
              for q2=1:7;
                  if ismember(FirstNeighbours(i,j),Ind_CellSelection{q2})
                    FirstNeighboursCount{q}(j-1,q2)=FirstNeighboursCount{q}(j-1,q2)+1;
                  end
              end
          end
      end
  end
 end
 
 for q=1:7
     figure
     for i=1:5 % first 5 neighbours
     subplot(5,1,i)
     bar(FirstNeighboursCount{q}(i,:),'FaceColor',[0.9290 0.6940 0.1250])
     hold on,
        bar(FirstNeighboursCount{1}(i,:)*FirstNeighboursCount{q}(1,1)/FirstNeighboursCount{1}(1,1),'FaceAlpha',0,'LineWidth',2)
     xticklabels(criterion_short)
     ylabel(['Neighbour #',num2str(i)])
     end
     subplot(511)
     title(criterion_short{q})
 end
 
 figure
 for q=2:7
        subplot(2,3,q-1)  
      bar(FirstNeighboursCount{q}(1,:),'FaceColor',[0.9290 0.6940 0.1250])
     hold on,
        bar(FirstNeighboursCount{1}(1,:)*FirstNeighboursCount{q}(1,1)/FirstNeighboursCount{1}(1,1),'FaceAlpha',0,'LineWidth',2)
     xticklabels(criterion_short)
     ylabel(['Neighbour #',num2str(1)])
          title(criterion_short{q})
 end
 
 % with last bit commented out it's just the first neighbour counts in a
 % matrix
 for q=1:7
FirstNeighboursProb(q,:)=FirstNeighboursCount{q}(1,:); %/FirstNeighboursCount{q}(1,1);
 end
 
 figure
 imagesc(FirstNeighboursProb)
 xticklabels(criterion_short)
 yticklabels(criterion_short)
 
for q=2:6
    for q2=2:6
    CellBehaviorOverlap(q-1,q2-1)=length(intersect(Ind_CellSelection{q},Ind_CellSelection{q2}));
    CellBehaviorOverlapProb(q-1,q2-1)=length(intersect(Ind_CellSelection{q},Ind_CellSelection{q2}))/length(Ind_CellSelection{q});
    end
end

 figure
 imagesc(CellBehaviorOverlap)
 xticks([1:5])
 xticklabels({criterion_short{2:6}})
  yticks([1:5])
 yticklabels({criterion_short{2:6}})
 title('Absolute count of overlap between groups')
 
  figure
 imagesc(CellBehaviorOverlapProb)
 xticks([1:5])
 xticklabels({criterion_short{2:6}})
  yticks([1:5])
 yticklabels({criterion_short{2:6}})
  title('Probability that the row group belongs also to a column group')
 
  %% For sound as baseline stats
  
  FirstNeighboursCount_SoundAsBaseline=cell(4,1);
for q=1:4
    FirstNeighboursCount_SoundAsBaseline{q}=zeros(length(Ind_CellSelection_SoundAsBaseline{1})-1,4);
end


 for i=1:length(Ind_CellSelection_SoundAsBaseline{1}) %all cells
  [SortedDistance_SoundAsBaseline(i,:) FirstNeighbours_SoundAsBaseline(i,:)]=sort(Dist_betweenCells(i,:));
  for q=1:4; %make 7 for all the cells that do not respond
      if ismember(i,Ind_CellSelection_SoundAsBaseline{q});
          for j=2:length(Ind_CellSelection_SoundAsBaseline{1}); %1st cell will be identity, jth closest neighbour
              for q2=1:4;
                  if ismember(FirstNeighbours_SoundAsBaseline(i,j),Ind_CellSelection_SoundAsBaseline{q2})
                    FirstNeighboursCount_SoundAsBaseline{q}(j-1,q2)=FirstNeighboursCount_SoundAsBaseline{q}(j-1,q2)+1;
                  end
              end
          end
      end
  end
 end
 
 for q=1:4
     figure
     for i=1:5 % first 5 neighbours
     subplot(5,1,i)
     bar(FirstNeighboursCount_SoundAsBaseline{q}(i,:),'FaceColor',[0.9290 0.6940 0.1250])
     hold on,
        bar(FirstNeighboursCount_SoundAsBaseline{1}(i,:)*FirstNeighboursCount_SoundAsBaseline{q}(1,1)/FirstNeighboursCount_SoundAsBaseline{1}(1,1),'FaceAlpha',0,'LineWidth',2)
     xticklabels(criterion_SoundAsBaseline_short)
     ylabel(['Neighbour #',num2str(i)])
     end
     subplot(511)
     title(criterion_SoundAsBaseline_short{q})
 end
 
 figure
 for q=2:4
        subplot(2,3,q-1)  
      bar(FirstNeighboursCount_SoundAsBaseline{q}(1,:),'FaceColor',[0.9290 0.6940 0.1250])
     hold on,
        bar(FirstNeighboursCount_SoundAsBaseline{1}(1,:)*FirstNeighboursCount_SoundAsBaseline{q}(1,1)/FirstNeighboursCount_SoundAsBaseline{1}(1,1),'FaceAlpha',0,'LineWidth',2)
     xticklabels(criterion_SoundAsBaseline_short)
     ylabel(['Neighbour #',num2str(1)])
          title(criterion_SoundAsBaseline_short{q})
 end
 
 % with last bit commented out it's just the first neighbour counts in a
 % matrix
 for q=1:4
FirstNeighboursProb_SoundAsBaseline(q,:)=FirstNeighboursCount_SoundAsBaseline{q}(1,:); %/FirstNeighboursCount{q}(1,1);
 end
 
 figure
 imagesc(FirstNeighboursProb_SoundAsBaseline)
 xticklabels(criterion_SoundAsBaseline_short)
 yticklabels(criterion_SoundAsBaseline_short)
 
for q=2:4
    for q2=2:4
    CellBehaviorOverlap_SoundAsBaseline(q-1,q2-1)=length(intersect(Ind_CellSelection_SoundAsBaseline{q},Ind_CellSelection_SoundAsBaseline{q2}));
    CellBehaviorOverlapProb_SoundAsBaseline(q-1,q2-1)=length(intersect(Ind_CellSelection_SoundAsBaseline{q},Ind_CellSelection_SoundAsBaseline{q2}))/length(Ind_CellSelection_SoundAsBaseline{q});
    end
end

 figure
 imagesc(CellBehaviorOverlap_SoundAsBaseline)
 xticks([1:3])
 xticklabels({criterion_SoundAsBaseline_short{2:4}})
  yticks([1:3])
 yticklabels({criterion_SoundAsBaseline_short{2:4}})
 title('Absolute count of overlap between groups')

  figure
 imagesc(CellBehaviorOverlapProb_SoundAsBaseline)
 xticks([1:3])
 xticklabels({criterion_SoundAsBaseline_short{2:4}})
  yticks([1:3])
 yticklabels({criterion_SoundAsBaseline_short{2:4}})
  title('Probability that the row group belongs also to a column group')
  
    %% plot tissue with cell number
 
  figure
  imagesc(brighten(brighten(dat.mimg(:,:,2),1),0.5))
colormap gray
hold on
for kk = 1:size(spatialInfo.ROIs,2)
    h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',1,'edgecolor',[0.7 0.7 0.7]);
    text(mean(spatialInfo.ROIs{kk}(:,2)),mean(spatialInfo.ROIs{kk}(:,1)),num2str(kk),'Color','m')
end
axis tight equal
set(gca,'TickDir','out','FontSize',14)
xlabel('Dorsal-ventral (pixels)')
ylabel('Rostral-caudal (pixels)')
title({[num2str(date),' - ', mouse],['Green channel'] })
  
%% for one cell plot cell response with colours
% and do also figure with individual cell traces
% see how red cells are

%% make cell tissue plot for max stim conditions
 
micronsPerPixel=exptInfo.recInfo.micronsPerPixel(1);
%  
% if dim_relevant > 1
% figure
% a=[1 fliplr(ind_max_dim_relevant) length(stimInfo.index)];
%     hold on
% for ii=1:4
% subplot(2,2,ii)
% imagesc([1:size(raster_ordered,2)]/fr,[1:size(raster_ordered,3)],squeeze(mean_Trace(a(ii),:,:))')
% colormap gray
% set(gca,'TickDir','out','FontSize',14)
% xlabel('Time (s)')
% ylabel('Cells')
% 
% end
% end


%%%%

if length(dim_relevant) > 1
figure
a=[1 fliplr(ind_max_dim_relevant) length(stimInfo.index)];
    hold on
for ii=1:4
 subplot(2,2,ii)

imagesc(brighten(brighten(dat.mimg(:,:,2),1),0.5))
colormap gray
hold on



 for kk = 1:size(spatialInfo.ROIs,2) %all cells grey outline
    if Cells_ActivityChange{a(ii)}(kk)==1
            h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0.8,'edgealpha',0,'LineWidth',1,'facecolor',[1 0.6 0.6]); 
    elseif Cells_ActivityChange{a(ii)}(kk)==-1
             h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0.8,'edgealpha',0,'LineWidth',1,'facecolor',[0 0 1]);
    else
     h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',1,'edgecolor',[0.7 0.7 0.7]);
    end
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

      title({[num2str(date),'_',mouse,'_',num2str(Filename(16:end-4))], ...
          ['index 1 = ',num2str(stimInfo.index(a(ii),dim_relevant(1))),'; index 2 = ',num2str(stimInfo.index(a(ii),dim_relevant(2)))], ...
         ['n = ',num2str(size(raster,3)),' total cells, Red cells = red : n = ',num2str(length(Ind_CellSelection{2})), '/', num2str(size(raster,3))] ...
        ['Increase = red interior : n = ',num2str(length(find(Cells_ActivityChange{a(ii)}(:)==1))), '/', num2str(size(raster,3))], ...
        ['Decrease = blue interior : n = ',num2str(length(find(Cells_ActivityChange{a(ii)}(:)==-1))), '/', num2str(size(raster,3))] })
         

end
end


%% make a movie of averaged network responses ?

mov_mean_Trace=cell(size(mean_Trace,1),1);
movmat_mean_Trace=cell(size(mean_Trace,1),1);
Nframes=size(mean_Trace,2);
% cmap=colormap(gray(256));

maxc=max(max(max(abs(mean_Trace)))); %scaling of colormap

tic
for ii=1:size(mean_Trace,1)
% movmat_mean_Trace{ii}=zeros([size(dat.mimg(:,:,2),1),size(dat.mimg(:,:,2),2),1,Nframes],class(dat.mimg));
 movmat_mean_Trace{ii}=ceil(255/2)*ones(size(dat.mimg(:,:,2),1),size(dat.mimg(:,:,2),2),1,Nframes);

for jj=1:Nframes % for every time point
    movmat_mask=ceil(255/2)*ones(size(dat.mimg(:,:,2),1),size(dat.mimg(:,:,2),2));
    for kk=1:size(spatialInfo.ROIs,2) % for every cell
        % create a mask for each cell
     BW = poly2mask(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),size(dat.mimg(:,:,2),1),size(dat.mimg(:,:,2),2));
     BW_idx = find(BW);
     movmat_mask(BW_idx)= ceil(255/2)+mean_Trace(ii,jj,kk)*255/(2*maxc);
        
%         
%         for ll=1:length(spatialInfo.ROIs{kk}(:,2)), 
%             %for real value - but right now colormap problem
% %              movmat_mean_Trace{ii}(spatialInfo.ROIs{kk}(ll,1),spatialInfo.ROIs{kk}(ll,2),1,jj)=mean_Trace(ii,jj,kk);
%             %scaled between 0 and 1
%              movmat_mean_Trace{ii}(spatialInfo.ROIs{kk}(ll,1),spatialInfo.ROIs{kk}(ll,2),1,jj)=ceil(255/2)+mean_Trace(ii,jj,kk)*255/(2*maxc);        
%         end
        
% h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0.8,'edgealpha',0,'LineWidth',1,'facecolor',[1 0.6 0.6]); 
    
    
    end
    movmat_mean_Trace{ii}(:,:,1,jj)=movmat_mask;
end  
mov_mean_Trace{ii}=immovie(movmat_mean_Trace{ii},colormap(gray(256)));
end

toc

% Save videos

for i=1:size(mean_Trace,1)
v=VideoWriter([mouse,'_',num2str(date),'_graymov',num2str(i),'.avi'],'Uncompressed AVI'); v.FrameRate=30; open(v)
writeVideo(v,mov_mean_Trace{i})
close(v)
end

%% across laser or sound paramters, number of cells that increase, decrease

laseramp=unique(stimInfo.index(:,2));
soundamp=unique(stimInfo.index(:,1));

NCellChange=zeros(dim_index(1),dim_index(2),3); % third dimension is increase(+1) / none(0) / decrease(-1)

for i=1:size(stimInfo.index,1)
    a=find(stimInfo.index(i,1)==soundamp);
    b=find(stimInfo.index(i,2)==laseramp);
    NCellChange(a,b,1)=length(find(Cells_ActivityChange{i}(:)==1));
    NCellChange(a,b,2)=length(find(Cells_ActivityChange{i}(:)==0));
    NCellChange(a,b,3)=length(find(Cells_ActivityChange{i}(:)==-1));
end

figure
subplot(121)
hold on
for i=1:length(laseramp)
    plot(soundamp,squeeze(NCellChange(:,i,1)),'Color',[1 1 1]-i/length(laseramp)*[0 1 1]) %red
    plot(soundamp,squeeze(NCellChange(:,i,2)),'Color',1-i/length(laseramp)*[0.7 0.7 0.7]) %gray
    plot(soundamp,squeeze(NCellChange(:,i,3)),'Color',[1 1 1]-i/length(laseramp)*[1 1 0]) %blue
end
title('number of cells modulated at different laser powers')
xlabel('Sound amplitude (dB)')
ylabel('Number of cells')

subplot(122)
hold on
for i=1:length(soundamp)
    plot(laseramp,squeeze(NCellChange(i,:,1)),'Color',[1 1 1]-i/length(soundamp)*[0 1 1]) %red
    plot(laseramp,squeeze(NCellChange(i,:,2)),'Color',1-i/length(soundamp)*[0.7 0.7 0.7]) %gray
    plot(laseramp,squeeze(NCellChange(i,:,3)),'Color',[1 1 1]-i/length(soundamp)*[1 1 0]) %blue
end
title('number of cells modulated at different sound amplitudes')
xlabel('Laser amplitude (Vc)')
ylabel('Number of cells')

figure
title_fig=cell(3,1);
title_fig{1}='Increasing cells';
title_fig{2}='No mod';
title_fig{3}='Decreasing cells';

for i=1:3
subplot(3,length(dim_relevant),1+length(dim_relevant)*(i-1))
bar(soundamp,squeeze(NCellChange(:,:,i)))
title({'Different laser powers',title_fig{i}})
xlabel('Sound amplitude (dB)')
ylabel('Number of cells')

subplot(3,length(dim_relevant),2+length(dim_relevant)*(i-1))
bar(squeeze(NCellChange(:,:,i)'))
title({'Different sound amplitudes',title_fig{i}})
xlabel('Laser amplitude (Vc)')
ylabel('Number of cells')
end


%% plotting trajectories of individual cells

% plot +1 / 0 / -1 for cells that actually are modulated at all
% with offset to each cell to separate them

figure
for ii=1:length(soundamp)
    a=find(stimInfo.index(:,1)==soundamp(ii));
subplot(1,length(soundamp),ii)
hold on
for kk=1:size(av_raster,2)
    b=[];
    for jj=1:length(a)
        b=[b Cells_ActivityChange{a(jj)}(kk)];
    end
    if sum(abs(b))~=0
 plot(b+1/3*(kk-size(av_raster,2)/2)/(size(av_raster,2)))
    end
end
xlabel('laser amplitude')
ylabel('response to stim')
title({['Sound amplitude = ',num2str(soundamp(ii)),'dB'],'1=increase, -1=decrease'})
end

figure
for ii=1:length(laseramp)
    a=find(stimInfo.index(:,2)==laseramp(ii));
subplot(1,length(laseramp),ii)
hold on
for kk=1:size(av_raster,2)
    b=[];
    for jj=1:length(a)
        b=[b Cells_ActivityChange{a(jj)}(kk)];
    end
    if sum(abs(b))~=0
 plot(soundamp,b+1/3*(kk-size(av_raster,2)/2)/(size(av_raster,2)))
    end
end
xlabel('sound amplitude')
ylabel('response to stim')
title({['laser amplitude = ',num2str(laseramp(ii)),'dB'],'1=increase, -1=decrease'})
end


%% same figure but with line width and text box saying how many cells go
% from one to another


%%%%%%%%%%%% across laser

% make matrix (denombrer) of all possible trajectories
values = -1:1;                               %//  data
k = length(laseramp);                                      %//  data
n = numel(values);                          %//  number of values
trajec_laserspace = values(dec2base(0:n^k-1,n)-'0'+1); 
count_trajec_laserspace=zeros(size(trajec_laserspace,1),length(soundamp));

figure
for ii=1:length(soundamp)
    a=find(stimInfo.index(:,1)==soundamp(ii));
subplot(1,length(soundamp),ii)
hold on
for kk=1:size(av_raster,2)
    b=[];
    for jj=1:length(a)
        b=[b Cells_ActivityChange{a(jj)}(kk)];
    end
    if sum(abs(b))~=0
        c=find(trajec_laserspace(:,1)==b(1)&trajec_laserspace(:,2)==b(2)&trajec_laserspace(:,3)==b(3));
        count_trajec_laserspace(c,ii)=count_trajec_laserspace(c,ii)+1;
    end
end
%  plot(b+1/3*(kk-size(av_raster,2)/2)/(size(av_raster,2)))
c=colormap(lines(size(count_trajec_laserspace,1)));
for jj=1:size(count_trajec_laserspace,1)
    if count_trajec_laserspace(jj,ii)~=0
plot(trajec_laserspace(jj,:),'LineWidth',100*count_trajec_laserspace(jj,ii)/(size(av_raster,2)),'Color',c(jj,:))
x=mod(jj,2*(length(laseramp)-1));
% text(1.25+x,(0.75-0.5*)*trajec_laserspace(jj,1+mod(jj,length(laseramp)-1))+0.2*trajec_laserspace(jj,2+mod(jj,length(laseramp)-1)),num2str(count_trajec_laserspace(jj)))
 text(1.25+0.5*x,0.05+(1-0.25-0.5*mod(x,2))*trajec_laserspace(jj,1+floor(x/2))+(0.25+0.5*mod(x,2))*trajec_laserspace(jj,2+floor(x/2)),num2str(count_trajec_laserspace(jj,ii)),'EdgeColor',c(jj,:))

    end
end
xlabel('laser amplitude')
ylabel('response to stim')
ylim([-1.25 1.25])
title({['Sound amplitude = ',num2str(soundamp(ii)),'dB'],'1=increase, -1=decrease'})
end

%%%%%%%%%%%% across sound

% make matrix (denombrer) of all possible trajectories
values = -1:1;                               %//  data
k = length(soundamp);                                      %//  data
n = numel(values);                          %//  number of values
trajec_laserspace = values(dec2base(0:n^k-1,n)-'0'+1); 
count_trajec_laserspace=zeros(size(trajec_laserspace,1),length(laseramp));

figure
for ii=1:length(laseramp)
    a=find(stimInfo.index(:,2)==laseramp(ii));
subplot(1,length(laseramp),ii)
hold on
for kk=1:size(av_raster,2)
    b=[];
    for jj=1:length(a)
        b=[b Cells_ActivityChange{a(jj)}(kk)];
    end
    if sum(abs(b))~=0
        c=find(trajec_laserspace(:,1)==b(1)&trajec_laserspace(:,2)==b(2)); % not generalized - only if there are 2 soudn amps
        count_trajec_laserspace(c,ii)=count_trajec_laserspace(c,ii)+1;
    end
end
%  plot(b+1/3*(kk-size(av_raster,2)/2)/(size(av_raster,2)))
c=colormap(lines(size(count_trajec_laserspace,1)));
for jj=1:size(count_trajec_laserspace,1)
    if count_trajec_laserspace(jj,ii)~=0
plot(trajec_laserspace(jj,:),'LineWidth',100*count_trajec_laserspace(jj,ii)/(size(av_raster,2)),'Color',c(jj,:))
x=mod(jj,2*(length(soundamp)-1));
% text(1.25+x,(0.75-0.5*)*trajec_laserspace(jj,1+mod(jj,length(laseramp)-1))+0.2*trajec_laserspace(jj,2+mod(jj,length(laseramp)-1)),num2str(count_trajec_laserspace(jj)))
 text(1.25+0.5*x,0.05+(1-0.25-0.5*mod(x,2))*trajec_laserspace(jj,1+floor(x/2))+(0.25+0.5*mod(x,2))*trajec_laserspace(jj,2+floor(x/2)),num2str(count_trajec_laserspace(jj,ii)),'EdgeColor',c(jj,:))

    end
end
xlabel('sound amplitude')
xticks([1:length(soundamp)])
xticklabels({num2str(soundamp)})
ylabel('response to stim')
ylim([-1.25 1.25])
title({['Laser amplitude = ',num2str(laseramp(ii)),'%'],'1=increase, -1=decrease'})
end


%% Number of cells changing their behavior ebtween two stim conditions and direction of change

Cells_ActivityChange_BtwConditions=zeros(size(Cells_ActivityChange,1),size(Cells_ActivityChange,1),5);
label=cell(size(Cells_ActivityChange,1),1);
% 5 as last dimension for -2/-1/0/1/2 as difference between two states

for ii=1:size(Cells_ActivityChange,1)
    for jj=1:size(Cells_ActivityChange,1)
        a=Cells_ActivityChange{ii}(:)-Cells_ActivityChange{jj}(:);
        Cells_ActivityChange_BtwConditions(ii,jj,1)=length(find(a==-2));
        Cells_ActivityChange_BtwConditions(ii,jj,2)=length(find(a==-1));
        Cells_ActivityChange_BtwConditions(ii,jj,3)=length(find(a==0));
        Cells_ActivityChange_BtwConditions(ii,jj,4)=length(find(a==1));        
        Cells_ActivityChange_BtwConditions(ii,jj,5)=length(find(a==2));        
    end
     label{ii}=[num2str(stimInfo.index(ii,1)),'/',num2str(round(stimInfo.index(ii,2)*100))];
end

figure
for ii=1:5
subplot(2,5,ii)
imagesc(Cells_ActivityChange_BtwConditions(:,:,ii))
xticks([1:size(Cells_ActivityChange,1)])
xticklabels(label)
yticks([1:size(Cells_ActivityChange,1)])
yticklabels(label)
title({'activity change matrix (row - col)',['sorted with sound amps regrouped'],['change = ',num2str(ii-3)]})
end
%ordered with laser conditions regrouped
ord=[];
for jj=1:length(laseramp)
ord=[ord jj:length(laseramp):size(stimInfo.index,1)];
end

for ii=1:5
subplot(2,5,5+ii)
imagesc(Cells_ActivityChange_BtwConditions(ord,ord,ii))
xticks([1:size(Cells_ActivityChange,1)])
xticklabels({label{ord}})
yticks([1:size(Cells_ActivityChange,1)])
yticklabels({label{ord}})
title({'activity change matrix (row - col)',['sorted with laser amps regrouped'],['change = ',num2str(ii-3)]})
end

%% untouched from general analysis underneath

%% plot increase/decrease for red, non red and all cells
%figure 1 w/0 laser, figure 2 w/ laser

if exist('stimInfo.index_names')==0;
    stimInfo.index_names={'sound','opto'};
end

CellTracesInfo.dim_relevant=dim_relevant;
CellTracesInfo.index=stimInfo.index;
CellTracesInfo.av_raster=av_raster;
CellTracesInfo.fr=fr;
CellTracesInfo.mean_Trace_ord=mean_Trace;
CellTracesInfo.Cells_ActivityChange_Ord=Cells_ActivityChange;
CellTracesInfo.Ind_cells_red_sort=Ind_cells_red;
% CellTracesInfo.mean_Trace_ord=mean_Trace_ord;
% CellTracesInfo.Cells_ActivityChange_Ord=Cells_ActivityChange_Ord;
% CellTracesInfo.Ind_cells_red_sort=Ind_cells_red_sort;

if length(dim_relevant)==1;
        ind=[min(stimInfo.index(:,dim_relevant))];      
        makeFig_CellTraces_General_MT(CellTracesInfo,ind)
        hold on
        suptitle({[exptInfo.recDate,' - ', exptInfo.mouse, ' - ', regexprep(Filename(16:end-4), '_','-')], ...
            [relevant_info,'; index 1 = ',num2str(ind(1)),' (',char(stimInfo.index_names(1)),');'],...
            [' '],[' ']})
        
        ind=[max(stimInfo.index(:,dim_relevant))];      
        makeFig_CellTraces_General_MT(CellTracesInfo,ind)
        hold on
        suptitle({[exptInfo.recDate,' - ', exptInfo.mouse, ' - ', regexprep(Filename(16:end-4), '_','-')], ...
            [relevant_info,'; index 1 = ',num2str(ind(1)),' (',char(stimInfo.index_names(1)),');'],...
            [' '],[' ']})
else length(dim_relevant)==2;
    for i=min(stimInfo.index(:,dim_relevant(1))):max(stimInfo.index(:,dim_relevant(1)))-min(stimInfo.index(:,dim_relevant(1))):max(stimInfo.index(:,dim_relevant(1)));
        ind=[i,min(stimInfo.index(:,dim_relevant(2)))];      
        makeFig_CellTraces_General_MT(CellTracesInfo,ind)
        hold on
        suptitle({[exptInfo.recDate,' - ', exptInfo.mouse, ' - ', regexprep(Filename(16:end-4), '_','-')], ...
            [relevant_info,'; index 1 = ',num2str(ind(1)),' (',char(stimInfo.index_names(1)),'); index 2 = ',num2str(ind(2)),' (',char(stimInfo.index_names(2)),')'],...
            [' '],[' ']})
        
        ind=[i,max(stimInfo.index(:,dim_relevant(2)))];      
        makeFig_CellTraces_General_MT(CellTracesInfo,ind)
        hold on
        suptitle({[exptInfo.recDate,' - ', exptInfo.mouse, ' - ', regexprep(Filename(16:end-4), '_','-')], ...
            [relevant_info,'; index 1 = ',num2str(ind(1)),' (',char(stimInfo.index_names(1)),'); index 2 = ',num2str(ind(2)),' (',char(stimInfo.index_names(2)),')'],...
            [' '],[' ']})

    end
end
%% bar graph of sound responsive cells with and without laser



figure
suptitle({[exptInfo.recDate,' - ', exptInfo.mouse, ' - ', regexprep(Filename(16:end-4), '_','-')], ...
[' '],[' ']})
hold on
subplot(121)
hold on
%a=setdiff(find(Cells_ActivityChange_Ord{1}==1),intersect(Ind_cells_red_sort,find(Cells_ActivityChange_Ord{1}==1)));

plot([av_raster_poststim(1,Ind_sort_poststim(a)); av_raster_poststim(2,Ind_sort_poststim(a))],'o-')

bar([1 2],[mean(av_raster_poststim(1,Ind_sort_poststim(a))) mean(av_raster_poststim(2,Ind_sort_poststim(a)))],'FaceColor','none','LineWidth',1.2)
errorbar([1 2],[mean(av_raster_poststim(1,Ind_sort_poststim(a))) mean(av_raster_poststim(2,Ind_sort_poststim(a)))],[std(av_raster_poststim(1,Ind_sort_poststim(a)))/sqrt(length(a)) std(av_raster_poststim(2,Ind_sort_poststim(a)))/sqrt(length(a))],'k.','LineWidth',1.2)
[h2 p2] = ttest(av_raster_poststim(1,Ind_sort_poststim(a)), av_raster_poststim(2,Ind_sort_poststim(a)));
title({[' '],'average value for sound-responsive,', 'non-red cells with and without opto',['paired t-test : p = ',num2str(p2),'; n = ',num2str(length(a))]})

subplot(122)
hold on
a=setdiff(find(Cells_ActivityChange_Ord{1}==1),intersect(Ind_cells_red_sort,find(Cells_ActivityChange_Ord{1}==1)));
b=av_raster_poststim(2,Ind_sort_poststim(a))./av_raster_poststim(1,Ind_sort_poststim(a));
boxplot(b)
title({'sound-responsive, non-red cells', 'average value for individual ratio with and without opto'})
text(0.1,0,{['mean = ',num2str(round(mean(b),3))],['std = ',num2str(round(std(b),3))],['median = ',num2str(round(median(b),3))],})
xlim([0 1.5])

%% red cells before and after opto

figure
subplot(231)
hold on
a=intersect(Ind_cells_red_sort,find(Cells_ActivityChange_Ord{1}==1));
if length(a)~=0
plot([1:size(av_raster,1)]/fr,squeeze(mean_Trace_ord(1,:,a)),'LineWidth',0.5,'Color',[1 0.7 0.7])
plot([1:size(av_raster,1)]/fr,mean(squeeze(mean_Trace_ord(1,:,a)),2),'r-','LineWidth',1)
end
xlabel('Time (s)')
ylabel('Change in fluorescence \Delta F / F_0')
title({[exptInfo.recDate,' - ', exptInfo.mouse, ' - ', regexprep(Filename(16:end-4), '_','-')], ...
'All red cell that increase - index = 0', ...
['thick line = average, n = ',num2str(length(a)), '/', num2str(size(Ind_cells_red_sort,2))]})

subplot(232)
hold on
a=intersect(Ind_cells_red_sort,find(Cells_ActivityChange_Ord{1}==-1));
if length(a)~=0
plot([1:size(av_raster,1)]/fr,squeeze(mean_Trace_ord(1,:,a)),'LineWidth',0.5,'Color',[1 0.7 0.7])
plot([1:size(av_raster,1)]/fr,mean(squeeze(mean_Trace_ord(1,:,a)),2),'r-','LineWidth',1)
end
xlabel('Time (s)')
ylabel('Change in fluorescence \Delta F / F_0')
title({[exptInfo.recDate,' - ', exptInfo.mouse, ' - ', regexprep(Filename(16:end-4), '_','-')], ...
'All red cell that decrease - index = 0', ...
['thick line = average, n = ',num2str(length(a)), '/', num2str(size(Ind_cells_red_sort,2))]})

subplot(233)
hold on
a=intersect(Ind_cells_red_sort,find(Cells_ActivityChange_Ord{1}==0));
plot([1:size(av_raster,1)]/fr,squeeze(mean_Trace_ord(1,:,a)),'LineWidth',0.5,'Color',[1 0.7 0.7])
plot([1:size(av_raster,1)]/fr,mean(squeeze(mean_Trace_ord(1,:,a)),2),'r-','LineWidth',1)
xlabel('Time (s)')
ylabel('Change in fluorescence \Delta F / F_0')
title({[exptInfo.recDate,' - ', exptInfo.mouse, ' - ', regexprep(Filename(16:end-4), '_','-')], ...
'All red cell with no activity change - index = 0', ...
['thick line = average, n = ',num2str(length(a)), '/', num2str(size(Ind_cells_red_sort,2))]})

subplot(234)
hold on
a=intersect(Ind_cells_red_sort,find(Cells_ActivityChange_Ord{2}==1));
if length(a)~=0
plot([1:size(av_raster,1)]/fr,squeeze(mean_Trace_ord(2,:,a)),'LineWidth',0.5,'Color',[1 0.7 0.7])
plot([1:size(av_raster,1)]/fr,mean(squeeze(mean_Trace_ord(2,:,a)),2),'r-','LineWidth',1)
end
xlabel('Time (s)')
ylabel('Change in fluorescence \Delta F / F_0')
title({[exptInfo.recDate,' - ', exptInfo.mouse, ' - ', regexprep(Filename(16:end-4), '_','-')], ...
'All red cell that increase - index = 1', ...
['thick line = average, n = ',num2str(length(a)), '/', num2str(size(Ind_cells_red_sort,2))]})

subplot(235)
hold on
a=intersect(Ind_cells_red_sort,find(Cells_ActivityChange_Ord{2}==-1));
if length(a)~=0
plot([1:size(av_raster,1)]/fr,squeeze(mean_Trace_ord(2,:,a)),'LineWidth',0.5,'Color',[1 0.7 0.7])
plot([1:size(av_raster,1)]/fr,mean(squeeze(mean_Trace_ord(2,:,a)),2),'r-','LineWidth',1)
end
xlabel('Time (s)')
ylabel('Change in fluorescence \Delta F / F_0')
title({[exptInfo.recDate,' - ', exptInfo.mouse, ' - ', regexprep(Filename(16:end-4), '_','-')], ...
'All red cell that decrease - index = 1', ...
['thick line = average, n = ',num2str(length(a)), '/', num2str(size(Ind_cells_red_sort,2))]})

subplot(236)
hold on
a=intersect(Ind_cells_red_sort,find(Cells_ActivityChange_Ord{2}==0));
plot([1:size(av_raster,1)]/fr,squeeze(mean_Trace_ord(2,:,a)),'LineWidth',0.5,'Color',[1 0.7 0.7])
plot([1:size(av_raster,1)]/fr,mean(squeeze(mean_Trace_ord(2,:,a)),2),'r-','LineWidth',1)
xlabel('Time (s)')
ylabel('Change in fluorescence \Delta F / F_0')
title({[exptInfo.recDate,' - ', exptInfo.mouse, ' - ', regexprep(Filename(16:end-4), '_','-')], ...
'All red cell with no activity change - index = 1', ...
['thick line = average, n = ',num2str(length(a)), '/', num2str(size(Ind_cells_red_sort,2))]})



%% for test cells - not ordered

test_cell=17;
% test_cell=find(Ind_sort_poststim==1);
figure
hold on
for ii=1:length(stimInfo.index)
    if length(dim_relevant)==1
subplot(1,length(stimInfo.index),ii)
    else
subplot(dim_index(dim_relevant(1)),dim_index(dim_relevant(2)),ii)
    end

plot([1:size(raster_ordered,2)]/fr,mean_Trace(ii,:,test_cell))
end
%endfigure,

%repeats
% figure
% subplot(121)
% imagesc([1:size(raster_ordered,2)]/fr,[1:repeats],raster_ordered(1+repeats*stimcond:repeats*(stimcond+1),:,test_cell))
% xlabel('Time (s)')
% ylabel('Repeats')
% title({['Repeats for cell ',num2str(test_cell)],['index 1 = ',num2str(stimInfo.index(1,dim_relevant(1)))]})
% 
% subplot(122)
% imagesc([1:size(raster_ordered,2)]/fr,[1:repeats],raster_ordered(repeats+1:2*repeats,:,test_cell))
% xlabel('Time (s)')
% ylabel('Repeats')
% title({['Repeats for cell ',num2str(test_cell)],['index 1 = ',num2str(stimInfo.index(2,dim_relevant(1)))]})

figure
hold on
for ii=1:length(stimInfo.index)
    if length(dim_relevant)==1
subplot(1,length(stimInfo.index),ii)
    else
subplot(dim_index(dim_relevant(1)),dim_index(dim_relevant(2)),ii)
imagesc([1:size(raster_ordered,2)]/fr,[1:repeats],raster_ordered(1+repeats*(ii-1):repeats*ii,:,test_cell))
xlabel('Time (s)')
ylabel('Repeats')
% title({['Repeats for cell ',num2str(test_cell)],['index 1 = ',num2str(stimInfo.index(ii,dim_relevant(1)))],['index 2 = ',num2str(stimInfo.index(ii,dim_relevant(2)))]})
    end
end

%%%%%%% plot statistics over different window intervals
% and mean + sem of fluorescence trace

figure

stimcond=3;
smoothn=31;


    subplot(311)

Time=[0:1:size(raster,2)-1]/fr-preOnsetTime;
    hold on
    x = repmat([Time(1:end) fliplr(Time)],1,1)';
    y = [smooth(mean_Trace(stimcond,:,test_cell),smoothn,'sgolay')'+smooth(std_mean_Trace(stimcond,:,test_cell),smoothn,'sgolay')'./sqrt(repeats) fliplr(smooth(mean_Trace(stimcond,:,test_cell),smoothn,'sgolay')'-smooth(std_mean_Trace(stimcond,:,test_cell),smoothn,'sgolay')'./sqrt(repeats))]';
    ph = patch(x,y,1);
    set(ph,'FaceVertexCData',1, ...
    'EdgeColor','none', ...
    'FaceAlpha',0.4)
plot(Time,smooth(mean_Trace(stimcond,:,test_cell),smoothn,'sgolay'),'k')
    hold off
    xlim([Time(1) Time(end)])
    xlabel('Time (s)');
    ylabel('Fluorescence variation');
    title({[exptInfo.recDate,' - ', exptInfo.mouse], ['Cell ',num2str(test_cell),' - Stimulus condition ',num2str(stimcond),' - smoothed sgolay, span = ',num2str(smoothn)]})
    hold off
colormap(gray)

subplot(312)

semilogy(Tdelay_test,squeeze(pstats_time(stimcond,test_cell,:)),'k'),
xlim([Time(1) Time(end)])
    xlabel('Time (s)');
    ylabel('p-value over 1s window');

subplot(313)

plot(Tdelay_test,squeeze(abs(zscore(stimcond,test_cell,:))),'k'),
xlim([Time(1) Time(end)])
    xlabel('Time (s)');
    ylabel('dprime over 1s window');


%% plot cells with their stats in outline



figure
 micronsPerPixel=exptInfo.recInfo.micronsPerPixel(1);
      suptitle({[exptInfo.recDate,' - ', exptInfo.mouse, ' - ', regexprep(Filename(16:end-4), '_','-')], ...
         relevant_info}) 
 for ii=1:dim_index
     subplot(1,dim_index,ii)
 imagesc(brighten(brighten(dat.mimg(:,:,2),1),0.5))
colormap gray
hold on
% c = colormap('jet');
 %   c = c(1:floor(size(c,1)/length(uModules{jj})):floor(size(c,1)/length(uModules{jj}))*length(uModules{jj}),:);
 %   colormap gray
 if exist('RedCells')==1
     for kk=1:size(Ind_cells_red_sort,2)
         h = patch(spatialInfo.ROIs{Ind_sort_poststim(Ind_cells_red_sort(kk))}(:,2),spatialInfo.ROIs{Ind_sort_poststim(Ind_cells_red_sort(kk))}(:,1),'r','facealpha',0.8,'edgealpha',0,'LineWidth',0.5,'facecolor',[1 0.6 0.6]);
     end
 end
 for kk = 1:size(spatialInfo.ROIs,2)
     if Cells_ActivityChange_Ord{ii}(kk)==1;
         h = patch(spatialInfo.ROIs{Ind_sort_poststim(kk)}(:,2),spatialInfo.ROIs{Ind_sort_poststim(kk)}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',2,'edgecolor',[1 0 0]);
     elseif Cells_ActivityChange_Ord{ii}(kk)==-1;
         h = patch(spatialInfo.ROIs{Ind_sort_poststim(kk)}(:,2),spatialInfo.ROIs{Ind_sort_poststim(kk)}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',2,'edgecolor',[0 0 1]);
     else
         h = patch(spatialInfo.ROIs{Ind_sort_poststim(kk)}(:,2),spatialInfo.ROIs{Ind_sort_poststim(kk)}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',1,'edgecolor',[0.7 0.7 0.7]);
     end
 end

axis tight equal
set(gca,'TickDir','out','FontSize',14)
xlabel('Dorsal-ventral (pixels)')
ylabel('Rostral-caudal (pixels)')
title({['index = ',num2str(stimInfo.index(ii)),' (opto amplitude)' ], ...
        ['Increased activity = red cells n = ',num2str(length(find(Cells_ActivityChange_Ord{ii}==1))), '/', num2str(size(Cells_ActivityChange_Ord{ii},2))], ...
        ['Decreased activity = blue cells n = ',num2str(length(find(Cells_ActivityChange_Ord{ii}==-1))), '/', num2str(size(Cells_ActivityChange_Ord{ii},2))] })
      if exist('RedCells')==1
      title({['index = ',num2str(stimInfo.index(ii)),' (opto amplitude)' ], ...
        ['Increased activity = red cells n = ',num2str(length(find(Cells_ActivityChange_Ord{ii}==1))), '/', num2str(size(Cells_ActivityChange_Ord{ii},2))], ...
        ['Decreased activity = blue cells n = ',num2str(length(find(Cells_ActivityChange_Ord{ii}==-1))), '/', num2str(size(Cells_ActivityChange_Ord{ii},2))], ...
        ['Red cells = pink interior n = ',num2str(length(Ind_cells_red_sort)), '/', num2str(size(Cells_ActivityChange_Ord{ii},2))] })
      end

 end
 
 
 %% fluorescence value of red cell versus response to laser stim
 
figure,
hold on,
for i=1:dim_index(2)
plot(RedCells.AvGrayValue_red(Ind_cells_red),av_raster_poststim(i,Ind_cells_red),'+','Color',(i-1)/dim_index(2)*[1 1 1])
end
 
 