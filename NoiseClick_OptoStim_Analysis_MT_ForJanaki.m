%% Noise click + optogenetics Analysis : 

% OUTLINE
% - load the data
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
%{
clear all
close all

% Select the file to open: e.g. MT068_20190830_FRA_OptoStim_03.mat
 
cd 'E:\2Pdata\data\'
[Filename,Pathname] = uigetfile('*.mat','Select the file to open');
load(sprintf('%s%s%s',Pathname,'\',Filename)); 
% Janaki: don't worry about tifStacks2
if exist([Pathname(1:end-5),Filename(1:5),'\',Filename(7:14),Filename(1:6),'tifStacks2'])~=0
    folder_proc='tifStacks2';
else
    folder_proc='tifStacks';
end
% load corresponding _proc file, e.g.:
% F_MT068_20190830MT068_tifStacks_plane1_proc.mat
procFile=['F_',Filename(1:14),Filename(1:5),'_',folder_proc,'_plane1_proc.mat'];
folder=num2str(str2num(Filename(end-5:end-4))); % to avoid changing code when the folder has below or above 10
load([Pathname(1:end-5),Filename(1:5),'\',Filename(7:14),Filename(1:6),folder_proc,'\',folder,'\',procFile]); 

mouse=Filename(1:5);
date=str2num(Filename(7:14));

% look if there's a RedCells.mat file, e.g:
% MT068_20190830_RedCells_SlopeRemoval_Movie.mat
% Janaki: as I changed my identification for red cells, the name of the
% redcell file changed. Now I'm using the RedCells_SlopeRemoval_Movie
RedCellfiles = dir(['E:\2Pdata\',mouse,'\',num2str(date),mouse,'_tifStacks','\**\RedCells*.mat']);
if length(RedCellfiles)~=0
    if any(strcmp({RedCellfiles.name}, 'RedCells_SlopeRemoval_Movie.mat'))
    i=find(strcmp({RedCellfiles.name}, 'RedCells_SlopeRemoval_Movie.mat')==1);
        if length(i)==1
        load([RedCellfiles(i).folder,'\RedCells_SlopeRemoval_Movie.mat'])
        RedCellLoadfile='RedCells_SlopeRemoval_Movie.mat';
        elseif length(i)>1
           [indx,tf] = listdlg('ListString',{RedCellfiles.folder},'PromptString',{'Select the folder',...
             'containing the Red Cell files'},'SelectionMode','single');
            i=indx;
        load([RedCellfiles(i).folder,'\RedCells_SlopeRemoval_Movie.mat'])
        RedCellLoadfile='RedCells_SlopeRemoval_Movie.mat';    
        end
    elseif any(strcmp({RedCellfiles.name}, 'RedCells_Movie.mat'))
    i=find(strcmp({RedCellfiles.name}, 'RedCells_Movie.mat')==1);
    load([RedCellfiles(i).folder,'\RedCells_Movie.mat'])
    RedCellLoadfile='RedCells_Movie.mat';
    else
    i=find(strcmp({RedCellfiles.name}, 'RedCells.mat')==1);
    load([RedCellfiles(i).folder,'\RedCells.mat'])  
    RedCellLoadfile='RedCells.mat';
    end
end
%}
%% calculate when laser bands are 
 % Janaki: when the laser shines, it makes bands of higher intensity in the
 % image, so this section is to find out where these bands are, which cells
 % are within the bands at that time frame (changes because of registration/movement), 
 % and replace the fluorescence at that time frame by a NaN. 

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
% Qs for melanie : does this conversion imply that the laser is basically
% scanning the frame? so we have a problem if the laser passes through a
% ROI?
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
            calcium.npilSubTraces(j,indfr_laser(i))=NaN;
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
postOnsetTime=round((events.eventsOn(2)-events.eventsOn(1))/fr)-1; % end raster x seconds after stim onset; here = duration between two stimuli = 5 seconds; 
doDFoverF0=1; % =1 means it does the z-scoring -> (F-mean(Fbaseline))/std(Fbaseline)
raster = makeCaRaster_NaN_MT(calcium.npilSubTraces,events.eventsOn,round(preOnsetTime*fr),round(postOnsetTime*fr),doDFoverF0);

% average over all events

av_raster = squeeze(nanmean(raster,1));



%% sort depending on components of index

% Janaki: uses function makeOrdRasterAndAveTrace_MT.m
% this rearranges the raster by stim combination and calculates average and
% std for each stim combination

repeats=stimInfo.repeats;

[raster_ordered mean_Trace std_mean_Trace]=makeOrdRasterAndAveTrace_MT(raster,stimInfo.index,stimInfo.order,stimInfo.repeats);
% the function makeOrdRasterAndAveTrace_MT.m replaces the following few
% lines; In the function, TrialOrder is not an output as it's not used
% afterwards but it could be modified to save it
% 
% for ii=1:length(stimInfo.index)
%     TrialOrder(ii,:) = find(stimInfo.order==ii);
%     raster_ordered(repeats*(ii-1)+1:repeats*ii,:,:)=raster(TrialOrder(ii,:),:,:);
% end
% 
% % average the trials for each test frequency
% 
% mean_Trace=zeros(length(stimInfo.index),size(raster,2),size(raster,3));
% 
% for ii=1:length(stimInfo.index)  
%     mean_Trace(ii,:,:) = nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,:,:),1);
%     std_mean_Trace(ii,:,:) = nanstd(raster_ordered(repeats*(ii-1)+1:repeats*ii,:,:),1);
% end



% number of components along each dimension of index
for i=1:size(stimInfo.index,2)
dim_index(i)=length(unique(stimInfo.index(:,i)));
end

dim_relevant=find(dim_index~=1);


% ordered in a  matrix, on several figures if necessary    

if length(dim_relevant)==2
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
    if dim_relevant==1; % condition added for only laser changing and sound always on
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


%%%% Option 1 : varying interval then take peak of average; copied from
%%%% GeneralStimulus_Analysis_MT

% Tdelay_test=[0.1:0.03:2];
% Tdelay_windur=1; % window duration of 1s for averaging.

% for jj=1:size(raster,3)
%     a=[];
%     b1=[];
%     for ii=1:length(Tdelay_test)
%         a(:,ii)=nanmean(mean_Trace(:,round(fr*(preOnsetTime+Tdelay_test(ii))):round(fr*(preOnsetTime+Tdelay_test(ii)+Tdelay_windur)),jj),2);
%     end
%     b1=max(abs(a(1,:)));
%     Tdelay(jj)=min(Tdelay_test(find(abs(a(1,:))==b1)));
%     av_raster_poststim(:,jj)=a(:,find(abs(a(1,:))==b1,1));
%     for ii=1:length(stimInfo.index)  
%         b2=squeeze(nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,round(fr*(preOnsetTime+Tdelay(jj))):round(fr*(preOnsetTime+Tdelay(jj)+Tdelay_windur)),jj),2));
%         sem_raster_poststim(ii,jj)=nanstd(b2)/sqrt(repeats);
%     end
% end
% 
% for i=1:length(dim_relevant)
%   [a Ind_sort_poststim{i}]=sort(av_raster_poststim(ind_max_dim_relevant(i),:),'descend');
%   mean_Trace_ord{i}=mean_Trace(:,:,Ind_sort_poststim{i});
% end
% 
% stats along 1st/2nd dimension
% alpha=0.01; 
% for ii=1:length(stimInfo.index) 
%     for jj=1:size(raster,3)
%       [tstats(ii,jj) pstats(ii,jj)]=ttest(nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,1:round(fr*preOnsetTime),jj),2),nanmean(raster_ordered(repeats*(ii-1)+1:repeats*ii,round(fr*(preOnsetTime+Tdelay(jj))):round(fr*(preOnsetTime+Tdelay(jj)+Tdelay_windur)),jj),2),'Alpha',alpha);
%     end
% end

%%%% Option 2 : varying interval then take min of pstat t-test; adapated from
%%%% above; no correction to t-test because I then take only one window.

Tdelay_test=[0.1:0.03:1.5];
Tdelay_windur=1; % window duration of 1s for averaging.
alpha=0.01;
 
Tdelay=[];
tstats_time=[];
pstats_time=[];
zscore_time=[];

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
        av_raster_poststim(ii,jj)=nanmean(b2);
        sem_raster_poststim(ii,jj)=nanstd(b2)/sqrt(repeats);
    end
end

for i=1:length(dim_relevant)
    [a Ind_sort_poststim{i}]=sort(av_raster_poststim(ind_max_dim_relevant(i),:),'descend');
    mean_Trace_ord{i}=mean_Trace(:,:,Ind_sort_poststim{i});
end



%%%% Option 3 : fixed interval for averaging stats analysis

% here different from general_analysis on 20190128 : not taking a moving
% window to find strongest response, rather take a 1-s long window starting
% either at sound onset or just after (ie after opto too)

% Tdelay=[1.1 2.1]; %interval after stim ONSET over which stats analysis is done
% Tdelay=[0.5 1.5];
% Tdelay=[0.75 1.25];

% for jj=1:size(raster,3)
%        av_raster_poststim(:,jj)=mean(mean_Trace(:,round(fr*(preOnsetTime+Tdelay(1))):round(fr*(preOnsetTime+Tdelay(2))),jj),2);
% for ii=1:length(stimInfo.index)  
%     a=squeeze(mean(raster_ordered(repeats*(ii-1)+1:repeats*ii,round(fr*(preOnsetTime+Tdelay(1))):round(fr*(preOnsetTime+Tdelay(2)))+round(fr),jj),2));
%     sem_raster_poststim(ii,jj)=std(a)/sqrt(repeats);
% end
% end
% 
% for i=1:length(dim_relevant)
% [a Ind_sort_poststim{i}]=sort(av_raster_poststim(ind_max_dim_relevant(i),:),'descend');
% mean_Trace_ord{i}=mean_Trace(:,:,Ind_sort_poststim{i});
% end
% 
% % stats along 1st/2nd dimension 
% 
% alpha=0.01; 
% 
% for ii=1:length(stimInfo.index) 
% for jj=1:size(raster,3)
% [tstats(ii,jj) pstats(ii,jj)]=ttest(mean(raster_ordered(repeats*(ii-1)+1:repeats*ii,1:round(fr*preOnsetTime),jj),2),mean(raster_ordered(repeats*(ii-1)+1:repeats*ii,round(fr*(preOnsetTime+Tdelay(1))):round(fr*(preOnsetTime+Tdelay(2))),jj),2),'Alpha',alpha);
% end
% end

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
%     for i=2:length(b)
%         Ind_CellSelection{3}=union(Ind_CellSelection{3},find(Cells_ActivityChange{b(i)}==1));
%     end
    
    for jj=1:size(av_raster,2)
        if sum(abs(Cells_ActivityChange2(b,jj)))>0
            b1=find(Cells_ActivityChange2(b,jj)~=0);
            c=find(abs(av_raster_poststim(:,jj))==max(abs(av_raster_poststim(b(b1),jj))));
            if Cells_ActivityChange2(c,jj)==1;
                Ind_CellSelection{3}=union(Ind_CellSelection{3},jj);
            elseif Cells_ActivityChange2(c,jj)==-1;
                Ind_CellSelection{4}=union(Ind_CellSelection{4},jj);
            end
        end
    end
    Ind_CellSelection{3}=reshape(Ind_CellSelection{3},[],1);
    Ind_CellSelection{4}=reshape(Ind_CellSelection{4},[],1);
    
    criterion{3} = 'All sound-responsive cells, increasing';
    criterion_short{3} = 'Sound+';
%%%% 4) All sound-responsive cells, decreasing %%%%
%     b=find(stimInfo.index(:,2)==0);
%     Ind_CellSelection{4}=[];
%     for i=2:length(b)
%         Ind_CellSelection{4}=union(Ind_CellSelection{4},find(Cells_ActivityChange{b(i)}==-1));
%     end
    criterion{4} = 'All sound-responsive cells, decreasing'; 
        criterion_short{4} = 'Sound-';
%%%% 5) All laser-responsive cells, increasing %%%%
    b=find(stimInfo.index(:,1)==0 & stimInfo.index(:,2)~=0);
    Ind_CellSelection{5}=[];
%     for i=2:length(b)
%         Ind_CellSelection{5}=union(Ind_CellSelection{5},find(Cells_ActivityChange{b(i)}==1));
%     end
    
        
    for jj=1:size(av_raster,2)
        if sum(abs(Cells_ActivityChange2(b,jj)))>0
            b1=find(Cells_ActivityChange2(b,jj)~=0);
            c=find(abs(av_raster_poststim(:,jj))==max(abs(av_raster_poststim(b(b1),jj))));
            if Cells_ActivityChange2(c,jj)==1;
            Ind_CellSelection{5}=union(Ind_CellSelection{5},jj);
            elseif Cells_ActivityChange2(c,jj)==-1;
            Ind_CellSelection{6}=union(Ind_CellSelection{6},jj);
            end
        end
    end
    Ind_CellSelection{5}=reshape(Ind_CellSelection{5},[],1);
    Ind_CellSelection{6}=reshape(Ind_CellSelection{6},[],1);    
    
    criterion{5} = 'All laser-responsive cells, increasing';
        criterion_short{5} = 'Laser+';
%%%% 6) All laser-responsive cells, decreasing %%%%
%     b=find(stimInfo.index(:,1)==0);
%     Ind_CellSelection{6}=[];
%     for i=2:length(b)
%         Ind_CellSelection{6}=union(Ind_CellSelection{6},find(Cells_ActivityChange{b(i)}==-1));
%     end
    criterion{6} = 'All laser-responsive cells, decreasing';
        criterion_short{6} = 'Laser-';
%%%% 7) Cells that don't respond and are not red %%%%
    Ind_CellSelection{7}=[1:size(av_raster,2)];
    for q=2:6;
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

% %%%% 7) Cells responsive to particular condition %%%%
%     Ind_CellSelection{7}=find(Cells_ActivityChange{10}==1);%
%     criterion{7} = 'Cells responsive to condition 10';

%% Plot curve cell response vs sound amplitude for stim- responsive-cells
% figure with overlap individual sound-responsive cell curves, different
% laser for different panels
% figure with average for the different laser amplitudes overlapped



for q=[1:6];
    fig(q+1)=figure(q+1);
    if length(Ind_CellSelection{q})~=0; 
        % plot the response curves for all individual cells
        for ii=1:dim_index(2) % dim index 2 = laser
            a=[ii:dim_index(2):length(stimInfo.index)];
            subplot(1,dim_index(2)+1,ii)
            hold  on
            if length(Ind_CellSelection{q})==1;
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

