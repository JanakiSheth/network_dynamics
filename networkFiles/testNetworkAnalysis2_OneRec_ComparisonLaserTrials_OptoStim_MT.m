%% test code for one recording - NMF analysis + incorporating other info
% to be extended to more recordings

% loading coeff and subsets for the different laser powers, then finding
% the signif ones!

% load the test data

clear variables
close all

% have default axes color black and not grey
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

% Select the file to open - only present subset files

mouseline = 'SOM';

if mouseline == 'VIP'
    cd '/home/janaki/Dropbox/project_with_melanie/networkFiles/VIP_NMFDatasets/'
elseif mouseline == 'SOM'
    cd '/home/janaki/Dropbox/project_with_melanie/networkFiles/SOM_NMFDatasets/'
end
    [subset_filename,Pathname] = uigetfile('*_subset.mat','Select the file to open');
 
mouse=subset_filename(6:10);
date=subset_filename(12:19);
recnumb=subset_filename(21:22);

if ~contains(subset_filename,'PosNeg')==0
    is_PosNeg=1;
else
    is_PosNeg=0;
end

% find which laser conditions are available
ind=strfind(subset_filename,'laser');
files = dir([subset_filename(1:ind+4),'*_subset.mat']);
laserconditions=zeros(1,3);
for i=1:length(files)
    ind=strfind(files(i).name,'laser');
    laserconditions(i) = str2double(files(i).name(ind+5));
end

% load the subsets and coeffs for the laser conditions available
% subsets - matrix of edge weights per subset (size N edges x N subsets)
% coeffs - activity of each subset over time (size N subsets x N time points).
coeff=cell(3,1);
subset=cell(3,1);
for k=1:length(laserconditions)
    subset{laserconditions(k)}=load(files(k).name);
    ind=strfind(files(k).name,'_subset.mat');
    coeff_filename=[files(k).name(1:ind),'coeff.mat'];
    coeff{laserconditions(k)}=load(coeff_filename);   
end

% load functional connectivity data
if mouseline == 'VIP'
    folderpath=['/home/janaki/Dropbox/project_with_melanie/networkFiles/VIP_connectivityData'];
elseif mouseline == 'SOM'
    folderpath=['/home/janaki/Dropbox/project_with_melanie/networkFiles/SOM_connectivityData'];
end
cd(folderpath)

connectivity_data=cell(3,1);
for k=1:length(laserconditions)
    ind=strfind(files(k).name,'_subset.mat');
    connectivity_filename=dir(['data_preproc_NMF_',files(k).name(6:ind-14),'*.mat']);       
    connectivity_data{k}=load(connectivity_filename.name);
end

% load activity data

if mouseline == 'VIP'
    folderpath=['/home/janaki/Dropbox/project_with_melanie/networkFiles/VIP_ActivityData'];
elseif mouseline == 'SOM'
    folderpath=['/home/janaki/Dropbox/project_with_melanie/networkFiles/SOM_ActivityData'];
end
cd(folderpath)
activityfile = dir(['workspace_',mouse,'_',date,'*.mat']);
activityworkspace=load(activityfile.name);

cd(Pathname)

% if PosNeg dataset, separate pos and neg component of coeff and calculate
% relative expression
if is_PosNeg ==1
    for k=laserconditions
        coeff{k,1}.coeff_all=coeff{k,1}.coeff;
        coeff{k,1}.coeff_pos=coeff{k,1}.coeff_all(:,1:end/2);
        coeff{k,1}.coeff_neg=coeff{k,1}.coeff_all(:,end/2+1:end);
        coeff{k,1}.coeff=coeff{k}.coeff_pos-coeff{k}.coeff_neg;    

        coeff{k,1}.mean_relativecoeff=mean(coeff{k,1}.coeff,2);
        for i=1:size(coeff{k,1}.coeff,1)
            coeff{k}.entropy(i)=wentropy(coeff{k,1}.coeff(i,:),'shannon');
        end
    end
end

%% plot weights of all subsets

% jet colormap - could replace 15 by max number of subsets in one laser
% conditions
cs=jet(15); % other option: hsv

set(groot,'defaultAxesColorOrder','remove');
set(groot,'defaultAxesColorOrder',cs);

figure(1)
for i=laserconditions
    subplot(2,3,i)
    plot(subset{i}.subset','LineWidth',2)
    title({[mouseline,' - ',mouse,' - ',date,'-',recnumb], ...
        ['all subsets, laser power',num2str(i)]})
    xlabel('edge number')
    ylabel('weight')
    subplot(2,3,3+i)
    plot(coeff{i}.coeff')
    xlabel('frame number')
    ylabel('coeff')
end

if is_PosNeg ==1    
    figure(2)
    for k=laserconditions
        subplot(3,3,k),
        bar(coeff{k}.mean_relativecoeff)
        xlabel('subset number')
        ylabel('mean relative coeff')
        title({[mouseline,' - ',mouse,' - ',date,'-',recnumb], ...
            ['laser power ',num2str(k)],'Mean relative coefficient = positive - negative coefficients'})
        subplot(3,3,3+k),
        [B, indm]= sort(coeff{k}.mean_relativecoeff,'descend');
        bar(B)
        set(gca,'xticklabel',{num2str(indm)})
        xlabel('subset number')
        ylabel('mean relative coeff')
        title({['Mean relative coefficient, ranked'],...
            ['mean over all subsets: ',num2str(mean(B),3)]})
        subplot(3,3,6+k),
        bar(coeff{k}.entropy)
        xlabel('subset number')
        ylabel('mean relative coeff')
        title('Entropy over whole traces')  
    end
end

%% turn subset weights into matrices, plot them

% calculate number of cells from number of edges

ncells = 1/2*(1+sqrt(1+8*size(subset{laserconditions(1)}.subset,2))); % solution of quadratic equation / Nedges = N(N-1)/2;

% fill in matrix for every subset - normalise with L2 norm
% this subsection converts a row of weights for edges into a matrix
% of weights for individual cells.

for k=laserconditions % laser conditions
    for i=1:size(subset{k}.subset,1) % number of subsets
        B=tril(ones(ncells),-1);
        B(B==1)=subset{k}.subset(i,:)/norm(subset{k}.subset(i,:));
        subset{k}.subset_mat(i,:,:)=B+B';
        subset{k}.subset_matnorm(i)=norm(subset{k}.subset(i,:));
    end
end


% plot the matrices

for k=laserconditions % one figure per laser condition
    figure(100+k)
    for i=1:size(subset{k}.subset,1)
        subplot(3,5,i),
        hold on,
        imagesc(squeeze(subset{k}.subset_mat(i,:,:)))
        title({[mouseline,' - ',mouse,' - ',date,'-',recnumb], ...
        ['matrix of subset ',num2str(i),', laser power',num2str(k)]})
        xlabel('nodes')
        ylabel('nodes')
        axis tight
        set(gca,'YDir','reverse') 
    end
    colorbar;
end


%% Order the matrix cell order by activity

Ind_CellSelection=activityworkspace.Ind_CellSelection;
 %red cells
CellID.red=Ind_CellSelection{2};
% sound+ cells
CellID.soundplaserp=setdiff(intersect(Ind_CellSelection{3},Ind_CellSelection{5}),CellID.red); % sound+ and Laser+, not red
CellID.soundplasern=setdiff(intersect(Ind_CellSelection{3},Ind_CellSelection{6}),CellID.red); % sound+ and Laser-, not red
CellID.soundplaser0=setdiff(setdiff(setdiff(Ind_CellSelection{3},Ind_CellSelection{5}),Ind_CellSelection{6}),CellID.red); % sound+ and Laser0, not red
% sound- cells
CellID.soundnlaserp=setdiff(intersect(Ind_CellSelection{4},Ind_CellSelection{5}),CellID.red); % sound- and Laser+, not red
CellID.soundnlasern=setdiff(intersect(Ind_CellSelection{4},Ind_CellSelection{6}),CellID.red); % sound- and Laser-, not red
CellID.soundnlaser0=setdiff(setdiff(setdiff(Ind_CellSelection{4},Ind_CellSelection{5}),Ind_CellSelection{6}),CellID.red); % sound- and Laser0, not red
% sound0 cells
CellID.sound0laserp=setdiff(setdiff(setdiff(Ind_CellSelection{5},Ind_CellSelection{3}),Ind_CellSelection{4}),CellID.red); % sound0 and Laser+, not red
CellID.sound0lasern=setdiff(setdiff(setdiff(Ind_CellSelection{6},Ind_CellSelection{3}),Ind_CellSelection{4}),CellID.red); % sound0 and Laser-, not red
CellID.sound0laser0=setdiff(setdiff(setdiff(setdiff(setdiff(Ind_CellSelection{1}',Ind_CellSelection{3}),Ind_CellSelection{4}),Ind_CellSelection{5}),Ind_CellSelection{6}),CellID.red); % sound0 and Laser0, not red
% put all this info in numbered cell for looping over, later on:
CellID.group=cell(10,1);
CellID.group{1}=CellID.soundplasern;
CellID.group{2}=CellID.soundplaser0;
CellID.group{3}=CellID.soundplaserp;
CellID.group{4}=CellID.sound0laserp;
CellID.group{5}=CellID.soundnlaserp;
CellID.group{6}=CellID.soundnlaser0;
CellID.group{7}=CellID.soundnlasern;
CellID.group{8}=CellID.sound0lasern;
CellID.group{9}=CellID.sound0laser0;
CellID.group{10}=CellID.red;

Order=[];
Order=[Order; CellID.soundplasern; CellID.soundplaser0; CellID.soundplaserp]; % Sound+ cells: Laser-, laser0 then laser+
Order=[Order; CellID.sound0laserp; CellID.soundnlaserp]; % laser+ cells
Order=[Order; CellID.soundnlaser0; CellID.soundnlasern]; % sound- cells
Order=[Order; CellID.sound0lasern]; % laser- cells
Order=[Order; CellID.sound0laser0; CellID.red]; %non resp cells and red cells


% plot the matrices

for k=laserconditions % one figure per laser condition
    figure(90+k)
    for i=1:size(subset{k}.subset,1)
        subplot(3,5,i),
        hold on,
        imagesc(squeeze(subset{k}.subset_mat(i,Order,Order)))
        % y axis
        line([-1 -1],[1 length([CellID.soundplasern; CellID.soundplaser0; CellID.soundplaserp])],'Color',[1 0.9 0],'LineWidth',2) %sound+ cells
        line([-2 -2],[1+length([CellID.soundplasern; CellID.soundplaser0]) 1+length([CellID.soundplasern; CellID.soundplaser0])+length(setdiff(Ind_CellSelection{5},CellID.red))],'Color','r','LineWidth',2 ) % laser+ cells
        line([-1 -1],[1+length([CellID.soundplasern; CellID.soundplaser0; CellID.soundplaserp; CellID.sound0laserp]) 1+length([CellID.soundplasern; CellID.soundplaser0; CellID.soundplaserp; CellID.sound0laserp; CellID.soundnlaserp; CellID.soundnlaser0; CellID.soundnlasern])],'Color','g','LineWidth',2) % sound- cells
        line([-2 -2],[1+length([CellID.soundplasern; CellID.soundplaser0; CellID.soundplaserp; CellID.sound0laserp; CellID.soundnlaser0;]) 1+length([CellID.soundplasern; CellID.soundplaser0; CellID.soundplaserp; CellID.sound0laserp; CellID.soundnlaserp; CellID.soundnlaser0; CellID.soundnlasern; CellID.sound0lasern])],'Color','b','LineWidth',2) % laser- cells
        line([-2 -2],[1 length(CellID.soundplasern)],'Color','b','LineWidth',2) % laser- and sound+ cells  first part
        line([-3 -3],[length(Order)-length(CellID.red) length(Order)],'Color','r','LineWidth',3) % laser- and sound+ cells  first part
        % x axis
        line([1 length([CellID.soundplasern; CellID.soundplaser0; CellID.soundplaserp])],[-1 -1],'Color',[1 0.9 0],'LineWidth',2) %sound+ cells
        line([1+length([CellID.soundplasern; CellID.soundplaser0]) 1+length([CellID.soundplasern; CellID.soundplaser0])+length(setdiff(Ind_CellSelection{5},CellID.red))],[-2 -2],'Color','r','LineWidth',2 ) % laser+ cells
        line([1+length([CellID.soundplasern; CellID.soundplaser0; CellID.soundplaserp; CellID.sound0laserp]) 1+length([CellID.soundplasern; CellID.soundplaser0; CellID.soundplaserp; CellID.sound0laserp; CellID.soundnlaserp; CellID.soundnlaser0; CellID.soundnlasern])],[-1 -1],'Color','g','LineWidth',2) % sound- cells
        line([1+length([CellID.soundplasern; CellID.soundplaser0; CellID.soundplaserp; CellID.sound0laserp; CellID.soundnlaser0;]) 1+length([CellID.soundplasern; CellID.soundplaser0; CellID.soundplaserp; CellID.sound0laserp; CellID.soundnlaserp; CellID.soundnlaser0; CellID.soundnlasern; CellID.sound0lasern])],[-2 -2],'Color','b','LineWidth',2) % laser- cells
        line([1 length(CellID.soundplasern)],[-2 -2],'Color','b','LineWidth',2) % laser- and sound+ cells  first part
        line([length(Order)-length(CellID.red) length(Order)],[-3 -3],'Color','r','LineWidth',3) % laser- and sound+ cells  first part        
        
        
        title({[mouseline,' - ',mouse,' - ',date,'-',recnumb], ...
        ['matrix of subset ',num2str(i),', laser power',num2str(k)]})
        xlabel('nodes')
        ylabel('nodes')
        axis tight
        set(gca,'YDir','reverse') 
    end
end


%% average subset matrix by cell activity

for k=laserconditions
    subset{k}.subset_grouped=zeros(size(subset{k}.subset,1),10,10); % 10 is the number of groups
    for l=1:size(subset{k}.subset,1)
        for i=1:10
            if ~isempty(CellID.group{i})
                for j=1:10
                    if ~isempty(CellID.group{i})
                        subset{k}.subset_grouped(l,i,j)=mean2(subset{k}.subset_mat(l,CellID.group{i},CellID.group{j}));                            
                    end
                end
            end
        end
    end
end

% plot the matrices

for k=laserconditions % one figure per laser condition
    figure(95+k)
    for i=1:size(subset{k}.subset,1)
        subplot(3,5,i),
        hold on,
        imagesc(squeeze(subset{k}.subset_grouped(i,:,:)),[0 0.025])
        % y axis
        line([-1 -1],[0.5 3.5],'Color',[1 0.9 0],'LineWidth',2) %sound+ cells
        line([-2 -2],[2.5 5.5],'Color','r','LineWidth',2 ) % laser+ cells
        line([-1 -1],[4.5  7.5],'Color','g','LineWidth',2) % sound- cells
        line([-2 -2],[6.5 8.5],'Color','b','LineWidth',2) % laser- cells
        line([-2 -2],[0.5 1.5],'Color','b','LineWidth',2) % laser- and sound+ cells  first part
        line([-3 -3],[9.5 10.5],'Color','r','LineWidth',3) % laser- and sound+ cells  first part
        % x axis
        line([0.5 3.5],[-1 -1],'Color',[1 0.9 0],'LineWidth',2) %sound+ cells
        line([2.5 5.5],[-2 -2],'Color','r','LineWidth',2 ) % laser+ cells
        line([4.5 7.5],[-1 -1],'Color','g','LineWidth',2) % sound- cells
        line([6.5 8.5],[-2 -2],'Color','b','LineWidth',2) % laser- cells
        line([0.5 1.5],[-2 -2],'Color','b','LineWidth',2) % laser- and sound+ cells  first part
        line([9.5 10.5],[-3 -3],'Color','r','LineWidth',3) % laser- and sound+ cells  first part        
        
        
        title({[mouseline,' - ',mouse,' - ',date,'-',recnumb], ...
    ['matrix of subset ',num2str(i),', laser power',num2str(k)]})
xlabel('nodes')
ylabel('nodes')
axis tight
set(gca,'YDir','reverse') 
colorbar
    end
end



%% to erase later? or modify anyways? but good to check on several datasets whether it's worth pursuing
% and find a way of labelling the interesting

nedges=size(subset{laserconditions(1)}.subset,2);

ang=@(a1,a2,b1,b2) 180/pi*acos(dot(subset{a1}.subset(a2,:),subset{b1}.subset(b2,:))/(norm(subset{a1}.subset(a2,:))*norm(subset{b1}.subset(b2,:))));

ind_angs=[];
angles=[];
dots=[];
normdots=[];
for k1=laserconditions
for i1=1:size(subset{k1}.subset,1)
    for k2=laserconditions
    for i2=1:size(subset{k2}.subset,1)
        if k2>k1
        ind_angs=[ind_angs; k1 i1 k2 i2];
        angles=[angles; real(ang(k1,i1,k2,i2))];
        dots=[dots; dot(subset{k1}.subset(i1,:),subset{k2}.subset(i2,:))];
        normdots=[normdots; dot(subset{k1}.subset(i1,:),subset{k2}.subset(i2,:))/(norm(subset{k1}.subset(i1,:))*norm(subset{k2}.subset(i2,:)))];
        end
    end
    end
end
end

a=find(dots>1.1/nedges); % need some work on that thresold: I had put it at 1.5 before, it excluded subset shared between laser1 and 2
b=find(angles(a)>0);
c=a(b);
selected_ind_angs=ind_angs(c,:);
selected_angles=angles(c);
selected_dots=dots(c);

figure, plot(dots*nedges,angles,'+')


%% for each laser power, create raster for sound stim, average them

% create raster
coeff_rast=cell(3,1);
for i=laserconditions
    for j=1:size(coeff{i}.coeff,2)/181
    coeff_rast{i}(j,1:181,:)=coeff{i}.coeff(:,181*(j-1)+1:181*j)';
    end
end

stimInfo=activityworkspace.stimInfo;

% average
for i=laserconditions
     index=stimInfo.index(find(stimInfo.index(:,2)==stimInfo.amplitude_opto(i)),:);
     order=1+floor((stimInfo.order(ismember(stimInfo.order,find(stimInfo.index(:,2)==stimInfo.amplitude_opto(i))))-1)/3);
[coeff_rastord{i} mean_coeff{i} std_mean_coeff{i}]=makeOrdRasterAndAveTrace_MT(coeff_rast{i},index,order,stimInfo.repeats);
end

%% plot the mean traces for every subset

co = [0.9 0.9 0.9;
    0.8 0.8 0.8;
    0.7 0.7 0.7;
    0.6 0.6 0.6;
    0.4 0.4 0.4;
    0.2 0.2 0.2;
    0 0 0];
set(groot,'defaultAxesColorOrder',co)

fr=activityworkspace.fr;

for i=laserconditions % one figure per laser condition
    figure(10+i)
    for j=1:size(mean_coeff{i},3)
        subplot(3,5,j),
        hold on,
        plot([1:181]/fr,mean_coeff{i}(:,:,j)')
        title({[mouseline,' - ',mouse,' - ',date,'-',recnumb], ...
    ['subset ',num2str(j),', laser power',num2str(i)]})
xlabel('Time (s)')
ylabel('Coeff')
xlim([0 181/fr])
    end
end

%% for each subset, find subsets that show significant stim-locked response

% oh yeah, here is the thing: with smoothing windows, there are like 8-9
% points (ie 250-300ms) before windowed functional connectivity starts rising - it would have been
% better to create raster in original dataset with 2 seconds prior to stim,
% not just one ...
% maybe this isn't so much the case with MTD?

% perform t test on each sound stim for each subset, for each laser
% condition

%%%% Let's do same analysis as for single cell: Option 2 from
%%%% NoiseClick_OptoStim_Analysis_MT

repeats=stimInfo.repeats;
preOnsetTime=2;

% %%%% Option 2 : varying interval then take min of pstat t-test; adapated from
% %%%% above; no correction to t-test because I then tak eonly one window.
% 
Tdelay_test=[0.1:0.03:1.5];
Tdelay_windur=1; % window duration of 1s for averaging.
alpha=0.01;
%  

coeff_tstats=cell(3,1);
for k=laserconditions % laser conditions
    coeff_tstats{k}.Tdelay=[];
    coeff_tstats{k}.tstats_time=zeros(length(stimInfo.intensity),size(coeff_rastord{k},3),length(Tdelay_test));
    coeff_tstats{k}.pstats_time=zeros(length(stimInfo.intensity),size(coeff_rastord{k},3),length(Tdelay_test));
    coeff_tstats{k}.zscore_time=zeros(length(stimInfo.intensity),size(coeff_rastord{k},3),length(Tdelay_test));
    for j=1:size(coeff_rastord{k},3) % number of subsets per laser condition
        for i=1:length(stimInfo.intensity) % sound conditions           
            % calculate stats with all Tdelay_test
            for l=1:length(Tdelay_test)
                [coeff_tstats{k}.tstats_time(i,j,l) coeff_tstats{k}.pstats_time(i,j,l)]=ttest(nanmean(coeff_rastord{k}(repeats*(i-1)+1:repeats*i,1:round(fr*Tdelay_windur),j),2),nanmean(coeff_rastord{k}(repeats*(i-1)+1:repeats*i,round(fr*(preOnsetTime+Tdelay_test(l))):round(fr*(preOnsetTime+Tdelay_test(l)+Tdelay_windur)),j),2),'Alpha',alpha);
                coeff_tstats{k}.zscore_time(i,j,l)=(mean(nanmean(coeff_rastord{k}(repeats*(i-1)+1:repeats*i,1:round(fr*Tdelay_windur),j),2))-mean(nanmean(coeff_rastord{k}(repeats*(i-1)+1:repeats*i,round(fr*(preOnsetTime+Tdelay_test(l))):round(fr*(preOnsetTime+Tdelay_test(l)+Tdelay_windur)),j),2)))/sqrt(0.5*(std(nanmean(coeff_rastord{k}(repeats*(i-1)+1:repeats*i,1:round(fr*Tdelay_windur),j),2))^2+std(nanmean(coeff_rastord{k}(repeats*(i-1)+1:repeats*i,round(fr*(preOnsetTime+Tdelay_test(l))):round(fr*(preOnsetTime+Tdelay_test(l)+Tdelay_windur)),j),2))^2)^2);
            end
            coeff_tstats{k}.Tdelay(i,j)=min(Tdelay_test(find(coeff_tstats{k}.pstats_time(i,j,:)==min(coeff_tstats{k}.pstats_time(i,j,:)))));
           
      % compute stats, average and std with chosen Tdelay
    [coeff_tstats{k}.tstats(i,j) coeff_tstats{k}.pstats(i,j)]=ttest(nanmean(coeff_rastord{k}(repeats*(i-1)+1:repeats*i,1:round(fr*Tdelay_windur),j),2),nanmean(coeff_rastord{k}(repeats*(i-1)+1:repeats*i,round(fr*(preOnsetTime+coeff_tstats{k}.Tdelay(i,j))):round(fr*(preOnsetTime+coeff_tstats{k}.Tdelay(i,j)+Tdelay_windur)),j),2),'Alpha',alpha);
    b2=squeeze(nanmean(coeff_rastord{k}(repeats*(i-1)+1:repeats*i,round(fr*(preOnsetTime+coeff_tstats{k}.Tdelay(i,j))):round(fr*(preOnsetTime+coeff_tstats{k}.Tdelay(i,j)+Tdelay_windur)),j),2));
    % probably put it in another structure
    coeff{k}.ave_poststim(i,j)=nanmean(b2); 
    coeff{k}.sem_poststim(i,j)=nanstd(b2)/sqrt(repeats);           
            
        end
    end
end


% plot the time time varying pstats

for i=laserconditions % one figure per laser condition
    figure(20+i)
    for j=1:size(mean_coeff{i},3) % number of subsets
        subplot(3,5,j), 
        semilogy(Tdelay_test,squeeze(coeff_tstats{i}.pstats_time(:,j,:))')
        hold on,
        line([Tdelay_test(1) Tdelay_test(end)],[alpha alpha],'Color','r')
        title({[mouseline,' - ',num2str(Filename(1:5)),'-',num2str(Filename(7:14)),'-',num2str(Filename(16:end-4))], ...
    ['subset ',num2str(j),', laser power',num2str(i)]})
xlabel('frame')
ylabel('p-stat')
xlim([0 Tdelay_test(end)])
    end
end


% find subsets that have significant stim-evoked responses

for k=1:length(stimInfo.amplitude_opto) % laser conditions
coeff_tstats{k}.isSignif=find(sum(coeff_tstats{k}.tstats)>0);
end


%% for each subset, find subsets that show AVERAGE significant stim-locked response
% copy paste the same analysis but instead of separating sound stim by
% sound stim - find subsets that are stim-evoked on average - might help
% discard the noisy ones and fish out some other ones?


% perform t test for average response to all sound stim for each subset, for each laser
% condition

%%%% Let's do same analysis as for single cell: Option 2 from
%%%% NoiseClick_OptoStim_Analysis_MT

repeats=stimInfo.repeats;
preOnsetTime=2;

% %%%% Option 2 : varying interval then take min of pstat t-test; adapated from
% %%%% above; no correction to t-test because I then tak eonly one window.
% 
Tdelay_test=[0.1:0.03:1.5];
Tdelay_windur=1; % window duration of 1s for averaging.
alpha=0.01;
%  


for k=laserconditions % laser conditions
    coeff_tstats{k}.Tdelay_ave=[];
    coeff_tstats{k}.tstats_ave_time=zeros(size(coeff_rastord{k},3),length(Tdelay_test));
    coeff_tstats{k}.pstats_ave_time=zeros(size(coeff_rastord{k},3),length(Tdelay_test));
    coeff_tstats{k}.zscore_ave_time=zeros(size(coeff_rastord{k},3),length(Tdelay_test));
    for j=1:size(coeff_rastord{k},3) % number of subsets per laser condition
%         for i=1:length(stimInfo.intensity) % sound conditions           
            % calculate stats with all Tdelay_test
            for l=1:length(Tdelay_test)
                [coeff_tstats{k}.tstats_ave_time(j,l) coeff_tstats{k}.pstats_ave_time(j,l)]=ttest(nanmean(coeff_rastord{k}(:,1:round(fr*Tdelay_windur),j),2),nanmean(coeff_rastord{k}(:,round(fr*(preOnsetTime+Tdelay_test(l))):round(fr*(preOnsetTime+Tdelay_test(l)+Tdelay_windur)),j),2),'Alpha',alpha);
                coeff_tstats{k}.zscore_ave_time(j,l)=(mean(nanmean(coeff_rastord{k}(:,1:round(fr*Tdelay_windur),j),2))-mean(nanmean(coeff_rastord{k}(:,round(fr*(preOnsetTime+Tdelay_test(l))):round(fr*(preOnsetTime+Tdelay_test(l)+Tdelay_windur)),j),2)))/sqrt(0.5*(std(nanmean(coeff_rastord{k}(:,1:round(fr*Tdelay_windur),j),2))^2+std(nanmean(coeff_rastord{k}(:,round(fr*(preOnsetTime+Tdelay_test(l))):round(fr*(preOnsetTime+Tdelay_test(l)+Tdelay_windur)),j),2))^2)^2);
            end
            coeff_tstats{k}.Tdelay_ave(j)=min(Tdelay_test(find(coeff_tstats{k}.pstats_ave_time(j,:)==min(coeff_tstats{k}.pstats_ave_time(j,:)))));
           
      % compute stats, average and std with chosen Tdelay
    [coeff_tstats{k}.tstats_ave(j) coeff_tstats{k}.pstats_ave(j)]=ttest(nanmean(coeff_rastord{k}(:,1:round(fr*Tdelay_windur),j),2),nanmean(coeff_rastord{k}(:,round(fr*(preOnsetTime+coeff_tstats{k}.Tdelay_ave(j))):round(fr*(preOnsetTime+coeff_tstats{k}.Tdelay_ave(j)+Tdelay_windur)),j),2),'Alpha',alpha);
    b2=squeeze(nanmean(coeff_rastord{k}(:,round(fr*(preOnsetTime+coeff_tstats{k}.Tdelay_ave(j))):round(fr*(preOnsetTime+coeff_tstats{k}.Tdelay_ave(j)+Tdelay_windur)),j),2));
    % probably put it in another structure
    coeff{k}.ave_avepoststim(j)=nanmean(b2); 
    coeff{k}.sem_avepoststim(j)=nanstd(b2)/sqrt(repeats);           
            
%         end
    end
end


% plot figure with average stim-evoked response + sem

for i=laserconditions % one figure per laser condition
    figure(30+i)
    for j=1:size(mean_coeff{i},3)
        subplot(3,5,j),
        hold on,
        plot([1:181]/fr,squeeze(mean(coeff_rastord{i}(:,:,j),1))','k')
        plot([1:181]/fr,squeeze(mean(coeff_rastord{i}(:,:,j),1))'+squeeze(std(coeff_rastord{i}(:,:,j),0,1))'./sqrt(size(coeff_rastord{i},1)),'Color',[0.8 0.8 0.8])
        plot([1:181]/fr,squeeze(mean(coeff_rastord{i}(:,:,j),1))'-squeeze(std(coeff_rastord{i}(:,:,j),0,1))'./sqrt(size(coeff_rastord{i},1)),'Color',[0.8 0.8 0.8])
        title({[mouseline,' - ',num2str(Filename(1:5)),'-',num2str(Filename(7:14)),'-',num2str(Filename(16:end-4))], ...
    ['subset ',num2str(j),', laser power',num2str(i)]})
xlabel('Time (s)')
ylabel('MEAN Coeff + SEM')
xlim([0 181/fr])
    end
end

% plot the time time varying pstats

for i=laserconditions % one figure per laser condition
    figure(40+i)
    for j=1:size(mean_coeff{i},3) % number of subsets
        subplot(3,5,j), 
        semilogy(Tdelay_test,squeeze(coeff_tstats{i}.pstats_ave_time(j,:))','k')
        hold on,
        line([Tdelay_test(1) Tdelay_test(end)],[alpha alpha],'Color','r')
        title({[mouseline,' - ',num2str(Filename(1:5)),'-',num2str(Filename(7:14)),'-',num2str(Filename(16:end-4))], ...
    ['subset ',num2str(j),', laser power',num2str(i)]})
xlabel('Time (s)')
ylabel('p-stat')
xlim([0 Tdelay_test(end)])
    end
end


% find subsets that have significant stim-evoked responses

for k=laserconditions % laser conditions
coeff_tstats{k}.isSignif_ave=find(coeff_tstats{k}.tstats_ave>0);
end


%% Average relative coeff expression of each subset per laser condition

set(groot,'defaultAxesColorOrder',cs);

figure(4)
for i=laserconditions % one subplot per laser condition
    for j=1:size(mean_coeff{i},3)
        subplot(1,3,i),
        hold on,
        plot([1:181]/fr,squeeze(mean(coeff_rastord{i}(:,:,j),1))','LineWidth',2)
        title({[mouseline,' - ',num2str(Filename(1:5)),'-',num2str(Filename(7:14)),'-',num2str(Filename(16:end-4))], ...
    ['laser power',num2str(i)]})
xlabel('Time (s)')
ylabel('MEAN Coeff + SEM')
xlim([0 181/fr])
    end
end


%% entropy on mean coeff trace

if is_PosNeg ==1;
    for k=laserconditions
    for i=1:size(coeff_rastord{k},3)
    coeff{k}.entropy_overmeantrace(i)=wentropy(mean(coeff_rastord{k}(:,:,i),1),'shannon');
    coeff{k}.entropy_overmeantrace_norm(i)=wentropy(mean(coeff_rastord{k}(:,:,i),1)/squeeze(max(max(abs(mean(coeff_rastord{k},1))))),'shannon');    
    coeff{k}.entropy_overmeantrace_norm2(i)=entropy((mean(coeff_rastord{k}(:,:,i),1)-min(mean(coeff_rastord{k}(:,:,i),1)))/squeeze(max(mean(coeff_rastord{k}(:,:,i),1))-min(mean(coeff_rastord{k}(:,:,i),1))));    

    end
    end
end

figure(6)
for k=laserconditions
subplot(3,3,k),
    bar(coeff{k}.entropy_overmeantrace)
    xlabel('subset number')
    ylabel('entropy over mean trace')
       title({[mouseline,' - ',mouse,' - ',date,'-',recnumb], ...
        ['laser power ',num2str(k)],...
    ['Entropy over mean trace']})
    subplot(3,3,3+k),
    bar(coeff{k}.entropy_overmeantrace_norm)
    xlabel('subset number')
    ylabel('entropy norm')
    title({['Entropy over normalised mean trace']})
        subplot(3,3,6+k),
    bar(coeff{k}.entropy_overmeantrace_norm2)
    xlabel('subset number')
    ylabel('entropy norm 2')
    title({['Entropy over normalised mean trace']})
end

%% on subset weights: measure for several thresholds number of edges and number of nodes
% over all subsets per laser condition

thres_weight=[0:0.1:100]; % threshold in folds of mean subset weight (=1/Nedges)
nedgesperthres_allsubsperlaser=zeros(length(thres_weight),3); % total number of edges above threshold for each laser condition
nnodesperthres_allsubsperlaser=zeros(length(thres_weight),3); % total number of nodes above threshold for each laser condition

nedgesperthres_allsubsperlaser_peredge=zeros(length(thres_weight),3,nedges);
nnodesperthres_allsubsperlaser_peredge=zeros(length(thres_weight),3,ncells);

nedges=size(subset{laserconditions(1)}.subset,2);

for k=laserconditions
    nsubsets(k)=size(subset{k}.subset,1);
end

% go back to tracking index for edges
ind=1;
for i=1:ncells
    for j=i+1:ncells
    ind_track(ind,:)=[i j];
    ind=ind+1;
    end
end
% tricky thing: not the same number of subsets - so number of edges at zero
% threshold won't be the same

for k=laserconditions
    for i=1:length(thres_weight)
    a=subset{k}.subset;
    a(a<thres_weight(i)/nedges)=0;
    a2=sum(a,1);
    a2(a2>0)=1;
    nedgesperthres_allsubsperlaser_peredge(i,k,:)=a2';
    nedgesperthres_allsubsperlaser(i,k)=length(find(a));
    % find index of non zero edges
    [b c]=find(a);  
    nnodesperthres_allsubsperlaser_peredge(i,k,unique(ind_track(c,:)))=1;
    nnodesperthres_allsubsperlaser(i,k)=length(unique(ind_track(c,:)));   
    end
end

cm=[0 0 0;
    1 0.7 0.7;
    1 0 0];

set(groot,'defaultAxesColorOrder',cm)

figure(5), 
subplot(221),
loglog(thres_weight,nedgesperthres_allsubsperlaser,'LineWidth',2)
        title({[mouseline,' - ',mouse,' - ',date,'-',recnumb], ...
    ['Total number of edges above threshold'],...
    ['(across subsets) per laser condition']})
xlabel('threshold / nedges')
ylabel('Number of edges above threshold')
legend('laser 1', 'laser 2', 'laser 3')
subplot(222),
semilogx(thres_weight,nnodesperthres_allsubsperlaser,'LineWidth',2)
        title({[mouseline,' - ',mouse,' - ',date,'-',recnumb], ...
    ['Total number of nodes above threshold'],...
    ['(across subsets) per laser condition']})
xlabel('threshold / nedges')
ylabel('Number of nodes above threshold')
legend('laser 1', 'laser 2', 'laser 3')
% percentage
subplot(223),
loglog(thres_weight,nedgesperthres_allsubsperlaser./(nedges*nsubsets),'LineWidth',2)
        title({[mouseline,' - ',mouse,' - ',date,'-',recnumb], ...
    ['Ratio of edges above threshold'],...
    ['(across subsets) per laser condition']})
xlabel('threshold / nedges')
ylabel('Ratio of edges above threshold')
legend('laser 1', 'laser 2', 'laser 3')
subplot(224),
semilogx(thres_weight,nnodesperthres_allsubsperlaser/ncells,'LineWidth',2)
        title({[mouseline,' - ',mouse,' - ',date,'-',recnumb], ...
    ['Ratio of nodes above threshold'],...
    ['(across subsets) per laser condition']})
xlabel('threshold / nedges')
ylabel('Ratio of nodes above threshold')
legend('laser 1', 'laser 2', 'laser 3')


figure,
for k=laserconditions
    subplot(3,2,2*(k-1)+1),
    imagesc(squeeze(nedgesperthres_allsubsperlaser_peredge(:,k,:))')
    
     subplot(3,2,2*k),
    imagesc(squeeze(nnodesperthres_allsubsperlaser_peredge(:,k,:))')   
    
end

%% on subset weights: measure for several thresholds number of edges and number of nodes
% for each subset in each laser condition

thres_weight=[0:0.1:100]; % threshold in folds of mean subset weight (=1/Nedges)


% tricky thing: not the same number of subsets - so number of edges at zero
% threshold won't be the same

for k=laserconditions
    subset{k}.nedgesperthres=zeros(length(thres_weight),size(subset{k}.subset,1)); % total number of edges above threshold for each subset
    subset{k}.nnodesperthres=zeros(length(thres_weight),size(subset{k}.subset,1)); % total number of nodes above threshold for each subset  
    for j=1:size(subset{k}.subset,1)
        for i=1:length(thres_weight)
        a=subset{k}.subset(j,:);
        a(a<thres_weight(i)/nedges)=0;
        subset{k}.nedgesperthres(i,j)=length(find(a));
        % find index of non zero edges
        [b c]=find(a);  
        subset{k}.nnodesperthres(i,j)=length(unique(ind_track(c,:)));
        end
    end
end


% clustering


for j=laserconditions
    subset{j}.clusteringperthres_bin=zeros(length(thres_weight),size(subset{j}.subset,1)); % binary clustering for each threshold per subset, mean over all nodes
    subset{j}.clusteringperthres_wei=zeros(length(thres_weight),size(subset{j}.subset,1)); % weighted clustering for each threshold per subset, mean over all nodes

for k=1:size(subset{j}.subset,1)
    for i=1:length(thres_weight)
    a=squeeze(subset{j}.subset_mat(k,:,:))*subset{j}.subset_matnorm(k);
    a(a<thres_weight(i)/nedges)=0;
    
    % clustering_coeff function gives the clustering per node - average
    % over all non-zero nodes? do both?
    b=clustering_coef_wu(a);
%     clusteringperthres_wei_nonzeromean(i,k)=mean(b(b~=0));
    subset{j}.clusteringperthres_wei(i,k)=mean(b);
    
    a2=a;
    a2(a>=thres_weight(i)/nedges)=1;
    b2=clustering_coef_bu(a2);
%     clusteringperthres_bin_nonzeromean(i,k)=mean(b2(b2~=0));
    subset{j}.clusteringperthres_bin(i,k)=mean(b2);
    end
end
end

%plot

set(groot,'defaultAxesColorOrder',cs);

figure(6), 
for k=laserconditions
subplot(4,3,k),
loglog(thres_weight,subset{k}.nedgesperthres,'LineWidth',2)
        title({[mouseline,' - ',mouse,' - ',date,'-',recnumb], ...
    ['laser ',num2str(k)],...
    ['Total number of edges above threshold']})
xlabel('threshold / nedges')
ylabel('Number of edges above threshold')
subplot(4,3,3+k),
semilogx(thres_weight,subset{k}.nnodesperthres,'LineWidth',2)
        title({['Total number of nodes above threshold']})
xlabel('threshold / nedges')
ylabel('Number of nodes above threshold')
subplot(4,3,6+k),
semilogx(thres_weight,subset{k}.clusteringperthres_bin,'LineWidth',2)
        title({['Binary clustering coefficient'],...
    ['per subset, mean over all nodes']})
xlabel('threshold / nedges')
ylabel('Binary clustering coeff')
xlabel('threshold / nedges')
ylabel('Number of nodes above threshold')
subplot(4,3,9+k),
semilogx(thres_weight,subset{k}.clusteringperthres_wei,'LineWidth',2)
        title({['Weighted clustering coefficient'],...
    ['per subset, mean over all nodes']})
xlabel('threshold / nedges')
ylabel('Weighted clustering coeff')
end

%% save a few figures: go to folder

% go to location

namefig_begin=[num2str(date),mouse,'_',recnumb,'_PerLaser_PosNeg'];



savefig(figure(96),[namefig_begin,'_SubsetMatrix_Grouped_laser1.fig'])
savefig(figure(97),[namefig_begin,'_SubsetMatrix_Grouped_laser2.fig'])
savefig(figure(98),[namefig_begin,'_SubsetMatrix_Grouped_laser3.fig'])


savefig(figure(1),[namefig_begin,'_SubsetsAndCoeffs.fig'])
savefig(figure(2),[namefig_begin,'_MeanRelativeCoeffs_BarGraph.fig'])
savefig(figure(4),[namefig_begin,'_MeanCoeffTrace.fig'])
savefig(figure(5),[namefig_begin,'_EdgesNodesAboveThreshold_PerLaser.fig'])
savefig(figure(6),[namefig_begin,'_NetworkMeasures_Clustering.fig'])

savefig(figure(11),[namefig_begin,'_CoeffTraces_PerSoundAmp_laser1.fig'])
savefig(figure(12),[namefig_begin,'_CoeffTraces_PerSoundAmp_laser2.fig'])
savefig(figure(13),[namefig_begin,'_CoeffTraces_PerSoundAmp_laser3.fig'])

savefig(figure(31),[namefig_begin,'_CoeffTraces_MeanOverSoundAmp_laser1.fig'])
savefig(figure(32),[namefig_begin,'_CoeffTraces_MeanOverSoundAmp_laser2.fig'])
savefig(figure(33),[namefig_begin,'_CoeffTraces_MeanOverSoundAmp_laser3.fig'])

%% select subsets on criterion on subset weights - basically discard subsets 
% where all weights are between n standard deviations from the mean.

% calculate mean and std for every subset / for all subsets at once?
% decision: for all subsets at once, all laser conditions
% sum (weights) = 1 so mean should always be the same = 1/n_edges
% pick out subsets with edges with a weight about a given threshold

thres_subset=4; %5

% calculate a grand std over all subsets, all laser conditions
a=[];
b=[];
for i=laserconditions
    a=[a std(subset{i}.subset')];
   b=[b reshape(subset{i}.subset,1,[])];
end
subset_avestd=std(b);
subset_avestd2=mean(a);
subset_mean=1/size(subset{1}.subset,2);


for i=laserconditions % laser conditions
    subset{i}.subset_thres=cell(size(subset{i}.subset,1),1);
%     subset{i}.subset_threspersub=cell(size(subset{i}.subset,1),1);
%     subset{i}.subset_thresperlaser=cell(size(subset{i}.subset,1),1);
%     subset{i}.subset_thresperiqr=cell(size(subset{i}.subset,1),1);
%     subset{i}.subset_avestdperlaser=mean(std(subset{i}.subset'));
    for j=1:size(subset{i}.subset,1) % number of subsets        
    % threshold with std over all subsets all alser ocnditions (std(all edges))
        subset{i}.subset_thres{j}=find(subset{i}.subset(j,:)>subset_mean+thres_subset*subset_avestd);
%     % threshold per laser condition
%     subset{i}.thresperlaser=subset_mean+thres_subset*subset{i}.subset_avestdperlaser;
%     subset{i}.subset_thresperlaser{j}=find(subset{i}.subset(j,:)>subset_mean+thres_subset*subset{i}.subset_avestdperlaser);
%         % threshold per subset
%         subset{i}.threspersub(j)=subset_mean+thres_subset*std(subset{i}.subset(j,:));
%         subset{i}.subset_threspersub{j}=find(subset{i}.subset(j,:)>subset_mean+thres_subset*std(subset{i}.subset(j,:)));
%     % using the interquartile distance
%         subset{i}.thresperiqr(j)=prctile(subset{i}.subset(j,:),75)+1.5*iqr(subset{i}.subset(j,:));
%         subset{i}.subset_thresperiqr{j}=find(subset{i}.subset(j,:)>subset{i}.thresperiqr(j));
%     
    end
    figure(1)
    subplot(2,3,i)
    hold on,
    line([1 size(subset{1}.subset,2)],[subset_mean+thres_subset*subset_avestd subset_mean+thres_subset*subset_avestd],'Color','r')

end


%% calculate a number of measures on the subsets: 1) clustering

% data to calculate this on: dataset used for nmf algo
% - abs(functional connectivity) calculated from moving average = time
% dependent
% - it's already divided by laser condition, all trials concatenated
% -> do I calculate measure then average? Or the opposite? Do I sort by
% sound amplitude or not?


% measures : 
% - size of subset - already dealt with in preceding section with subset
% thresholding
% - connection density (weighted?)
% - weight distribution, average edge weight within the subset
% - weight distribution/average outside subset, across subset (look at "parcellation")
% - clustering (weighted undirected, signed or not, starting with unsigned)


% Action plan:
% 1) make raster per sound amplitude for all edges, laser conditions
% 2) make average
% 3) take thresholded subsets, figure out which ones are non zero
% 3) calculate clustering

% also included the raster and averaging in there
% Now that the clustering has been calculated for thresholded subsets over
% all laser conditions, let's now make the rasters and average per sound
% and overall average.
% Question: should I do the average clustering per subset already or still
% keep it per node? -----> kept per node




% go back to tracking index for edges
ind=1;
for i=1:ncells
    for j=i+1:ncells
    ind_track(ind,:)=[i j];
    ind=ind+1;
    end
end


% do I do one matrix with plenty of zeros? Will it be more computationally
% demanding to do that to get rid of all nodes that are not involved in
% subset?
clustering=cell(length(stimInfo.amplitude_opto),1);

for i=laserconditions % laser conditions
    %     for kk=1:size(data{i}.data,1) % number of edges without repetition
    subset{i}.subset_thresnodes=cell(size(subset{i}.subset_thres));
    subset{i}.FC_matrix=cell(size(subset{i}.subset_thres,1),length(stimInfo.amplitude_opto));
    subset{i}.subset_sizenodes=zeros(size(subset{i}.subset_thres,1),1);
    subset{i}.subset_connectiondensity=zeros(size(subset{i}.subset_thres,1),1);
    clustering{i}.clustering=cell(size(subset{i}.subset_thres,1),length(stimInfo.amplitude_opto));
    clustering{i}.clust_rast=cell(size(subset{i}.subset_thres,1),length(stimInfo.amplitude_opto));
    clustering{i}.clust_rast=cell(size(subset{i}.subset_thres,1),length(stimInfo.amplitude_opto));
    clustering{i}.clust_rastord=cell(size(subset{i}.subset_thres,1),length(stimInfo.amplitude_opto));
    clustering{i}.mean_clust=cell(size(subset{i}.subset_thres,1),length(stimInfo.amplitude_opto));
    clustering{i}.std_mean_coeff=cell(size(subset{i}.subset_thres,1),length(stimInfo.amplitude_opto));
    for j=1:size(subset{i}.subset_thres,1) % subset j in laser condtiion i
            subset{i}.subset_thresnodes{j}=unique(ind_track(subset{i}.subset_thres{j}(:),:));
            subset{i}.subset_sizenodes(j)=length(subset{i}.subset_thresnodes{j});
        if length(subset{i}.subset_thres{j}) >1; % if thresholded subset is not empty or there is more than one edge (one edge = 2 nodes = no triangles for clustering)
            subset{i}.subset_connectiondensity(j)=length(subset{i}.subset_thres{j})*2/(subset{i}.subset_sizenodes(j)*(subset{i}.subset_sizenodes(j)-1));
            % right now : funct connectivity matrix with just the relevant nodes, but coudl maybe
            % be replaced with matrix with plenty of zeros for further
            % readability if I want to save these matrices later
            tic
            for l=laserconditions % laser conditions
            subset{i}.FC_matrix{j,l}=zeros(length(subset{i}.subset_thresnodes{j}),length(subset{i}.subset_thresnodes{j}),size(data{l}.data,2));
            for k=1:length(subset{i}.subset_thres{j})
                ind1 = find(subset{i}.subset_thresnodes{j}(:)==ind_track(subset{i}.subset_thres{j}(k),1));
                ind2 = find(subset{i}.subset_thresnodes{j}(:)==ind_track(subset{i}.subset_thres{j}(k),2));
                subset{i}.FC_matrix{j,l}(ind1,ind2,:)=data{l}.data(subset{i}.subset_thres{j}(k),:);
                subset{i}.FC_matrix{j,l}(ind2,ind1,:)=data{l}.data(subset{i}.subset_thres{j}(k),:);
            end
            % now that FC matrix is filled up, calculate clustering time point
            % per time point
            
            for m=1:size(subset{i}.FC_matrix{j},3) % total number of frames in that laser condition
                clustering{i}.clustering{j,l}(:,m)=clustering_coef_wu(subset{i}.FC_matrix{j,l}(:,:,m));
            end
            
            % create raster
            for jj=1:size(coeff{i}.coeff,2)/181 %number of trials per laser condtiion
            clustering{i}.clust_rast{j,l}(jj,1:181,:)=clustering{i}.clustering{j,l}(:,181*(jj-1)+1:181*jj)';
            end
            
            % calculate average
            index=stimInfo.index(find(stimInfo.index(:,2)==stimInfo.amplitude_opto(l)),:);
            order=1+floor((stimInfo.order(ismember(stimInfo.order,find(stimInfo.index(:,2)==stimInfo.amplitude_opto(l))))-1)/3);
            [clustering{i}.clust_rastord{j,l} clustering{i}.mean_clust{j,l} clustering{i}.std_mean_coeff{j,l}]=makeOrdRasterAndAveTrace_MT(clustering{i}.clust_rast{j,l},index,order,stimInfo.repeats);
            
            end
            toc
        end
    end
end


%% plot clustering per subset

% one figure per laser condition can be laser condtiions x thresholded subsets of the given laser condition and have the
% separate sound amplitudes

% another figure can be average clustering of subset for the different
% laser codntions

% first figure: each subplot has all sound stim traces for a given subset,
% given laser condition
for i=laserconditions % laser conditions
    figure(50+i)
    for j=1:size(clustering{i}.mean_clust,1) % number of subsets
        for l=laserconditions % number of laser conditions
            if size(clustering{i}.mean_clust{j,l})~[0 0];
                subplot(length(stimInfo.amplitude_opto),size(clustering{i}.mean_clust,1),(l-1)*size(clustering{i}.mean_clust,1)+j)               
                
                hold on,
                plot([1:181]/fr,squeeze(mean(clustering{i}.mean_clust{j,l},3))')
                title({[mouseline,' - ',num2str(Filename(1:5)),'-',num2str(Filename(7:14)),'-',num2str(Filename(16:end-4))], ...
                ['Average clustering over nodes'],...
                ['laser',num2str(i),' subset ',num2str(j),'; laser power',num2str(l)]})
                xlabel('frame')
                ylabel('Weighted clustering per sound amp')
                xlim([0 181/fr])
                
            end
            
            
        end
        
    end
end


% second figure: each subplot average stim response for a given subset,
% for the different laser condition

cm=[0 0 0;
    1 0.7 0.7;
    1 0 0];

for i=laserconditions % laser conditions
    figure(54+i)
    for j=1:size(clustering{i}.mean_clust,1) % number of subsets
        for l=laserconditions % number of laser conditions
            if size(clustering{i}.mean_clust{j,l})~[0 0];
                subplot(1,size(clustering{i}.mean_clust,1),j)    
                
                hold on,
                plot([1:181]/fr,squeeze(mean(mean(clustering{i}.clust_rast{j,l},3),1)),'Color',cm(l,:),'LineWidth',1)
                plot([1:181]/fr,squeeze(mean(mean(clustering{i}.clust_rast{j,l},3),1))+squeeze(std(mean(clustering{i}.clust_rast{j,l},3),0,1))./sqrt(size(clustering{i}.clust_rast{j,l},1)),'Color',cm(l,:),'LineStyle','--')
                plot([1:181]/fr,squeeze(mean(mean(clustering{i}.clust_rast{j,l},3),1))-squeeze(std(mean(clustering{i}.clust_rast{j,l},3),0,1))./sqrt(size(clustering{i}.clust_rast{j,l},1)),'Color',cm(l,:),'LineStyle','--')

%                 plot([1:181]/fr,squeeze(mean(coeff_rastord{i}(:,:,j),1))','k')
%         plot([1:181]/fr,squeeze(mean(coeff_rastord{i}(:,:,j),1))'+squeeze(std(coeff_rastord{i}(:,:,j),0,1))'./sqrt(size(coeff_rastord{i},1)),'Color',[0.8 0.8 0.8])
%         plot([1:181]/fr,squeeze(mean(coeff_rastord{i}(:,:,j),1))'-squeeze(std(coeff_rastord{i}(:,:,j),0,1))'./sqrt(size(coeff_rastord{i},1)),'Color',[0.8 0.8 0.8])

                title({[mouseline,' - ',num2str(Filename(1:5)),'-',num2str(Filename(7:14)),'-',num2str(Filename(16:end-4))], ...
                    ['Average clustering over nodes and sound amplitudes'],...
                    ['laser',num2str(i),' subset ',num2str(j),'; laser power',num2str(l)]})
                xlabel('frame')
                ylabel('Average weighted clustering')
                xlim([0 181/fr])
            end
        end
    end
end


%% code from previous analysis file with similar title that I'm modifying for this new purpose




%% plot figure with cells and subgraphs

% plot them first separately because the subgraphs overlap
% by placing loop elsewhere (see commented off lines), you can have the
% subgraphs overlap


 for ii=1:size(subset,1) % here for separate figures
    figure
 
%  micronsPerPixel=exptInfo.recInfo.micronsPerPixel(1);
%  imagesc(brighten(brighten(dat.mimg(:,:,2),1),0.5))
colormap gray
hold on

 for kk = 1:size(spatialInfo.ROIs,2) %all cells grey outline
     h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',1,'edgecolor',[0.7 0.7 0.7]);
 end
 
% for ii=[2 5 9];%1:size(subset,1) % here for overlap of subsets on one figure 
 
% define colour of filling
%%%% blue if stim traces - both laser and stim are in it
if subset(ii,1)>0 & subset(ii,2)>0;
    c=[0 0 1]; %blue
    info='blue: both stim in subgraph';
%%%% yellow if sound stim is in it but not laser
elseif subset(ii,1)>0 & subset(ii,2)==0;
    c=[1 1 0]; %yellow
    info='yellow: sound stim in subgraph';
%%%% red if laser is in it but not sound
elseif subset(ii,1)==0 & subset(ii,2)>0;
    c=[1 0 0]; %red
    info='red: laser stim in subgraph';
%%%% grey for the rest - noise subgraphs
elseif subset(ii,1)==0 & subset(ii,2)==0;
    c=[0.7 0.7 0.7]; %grey-black
    info='grey: no stim in subgraph';
end
    
a=find(subset(ii,3:end)); % traces 1 and 2 are for laser    
 for kk=1:length(a) 
     h = patch(spatialInfo.ROIs{a(kk)}(:,2),spatialInfo.ROIs{a(kk)}(:,1),'y','facealpha',0.6,'edgealpha',0,'LineWidth',1,'facecolor',c);
 end

  if exist('RedCells')==1 % red cells, red outline
 a=find(RedCells.isCell_red);     
     for kk=1:length(a)
         h = patch(spatialInfo.ROIs{a(kk)}(:,2),spatialInfo.ROIs{a(kk)}(:,1),'r','facealpha',0,'edgealpha',1,'LineWidth',2,'edgecolor',[1 0 0]);
     end
  end
 
% end % for overlap subsets on one figure 

axis tight equal
set(gca,'TickDir','out','FontSize',14)
xlabel('Dorsal-ventral (pixels)')
ylabel('Rostral-caudal (pixels)')

      title({[mouseline,' - ',num2str(Filename(1:5)),'-',num2str(Filename(7:14)),'-',num2str(Filename(16:end-4))], ...
         ['n = ',num2str(size(raster_ordered,3)),' total cells'], ...
         ['ng = ',num2str(size(subset,1)),' total subgraphs'], ...
         ['Subgraph ',num2str(ii),'; ',info,'; n = ',num2str(length(find(subset(ii,3:end)))), '/', num2str(size(raster_ordered,3))], ...
         ['Red cells = red : n = ',num2str(length(find(RedCells.isCell_red))), '/', num2str(size(raster_ordered,3))] }) 
 end % for separate figures

%% figure with distribution of how many subgraphs each cell belongs to

% stim traces excluded
for i=1:size(subset,2)-2 
    nsubbelong(i)=length(find(subset(:,i+2))); 
end

figure
subplot(121),
histogram(nsubbelong,[0:1:10],'FaceColor','none','LineWidth',2)
set(gca,'TickDir','out','FontSize',14)
xlabel('Number of subgraphs')
ylabel('Counts')
title({[mouseline,' - ',num2str(Filename(1:5)),'-',num2str(Filename(7:14)),'-',num2str(Filename(16:end-4))], ...
         ['Histogram of the number of subgraphs a cell belongs to'] }) 

subplot(122),
barwitherr(std(nsubbelong),mean(nsubbelong),'FaceColor','none','LineWidth',2)
set(gca,'TickDir','out','FontSize',14)
xlabel('(mean + std)')
ylabel('Number of subgraphs')
title({[mouseline,' - ',num2str(Filename(1:5)),'-',num2str(Filename(7:14)),'-',num2str(Filename(16:end-4))], ...
         ['Mean+std = ',num2str(round(mean(nsubbelong),1)),'+',num2str(round(std(nsubbelong),1))] }) 

    
%% Calculate correlations of each subgraph with stim traces

corr_stim=zeros(2,size(coeff,1));

for i=1:size(coeff,1)
    corr_stim(1,i)=corr(coeff(i,:)',data(1,:)');
    corr_stim(2,i)=corr(coeff(i,:)',data(2,:)');
end

figure,
bar_handle = bar(corr_stim','grouped');
bar_handle(1).FaceColor = [0.9290 0.6940 0.1250];
bar_handle(2).FaceColor = [0.8500 0.3250 0.0980];
xlabel('subgraph number')
ylabel('correlation')
title({['correlation with stimulus trace'],...
    ['sound = yellow / laser = red']})

%% Binomial test : are red cells more/less in subgraph?


% using the function myBinomTest.m found on wikipedia

for i=1:size(subset,1)
    s=length(intersect(find(RedCells.isCell_red),find(subset(i,3:end))));
    n=length(find(subset(i,3:end)));
    p=length(find(RedCells.isCell_red))*n/size(raster_ordered,3)^2;
    pout(i)=myBinomTest(s,n,p,'one');
end

figure
plot(pout,'ko-')
xlabel('subgraph number')
ylabel('p-value')
xlim([0 11])
ylim([0 1])
title({['Results from binomial test'], ...
    [num2str(length(find(pout<0.05))),' subgraphs with signif high/low number of red cells'], ...
    [num2str(find(pout<0.05))]})

%% Calculate sound+/- laser+/- cells and such to compare with subgraphs


% copied directly from Analysis_AllRec_NoiseClick_OptoStim_MT.m

Tdelay_test=[0.1:0.03:1.5];
Tdelay_windur=1; % window duration of 1s for averaging.
alpha=0.01;

[Tdelay pstats av_raster_poststim sem_raster_poststim]=findResponseWindowWithDelay_MT(raster_ordered,stimInfo.repeats,rasterInfo.preOnsetTime,fr,Tdelay_test,Tdelay_windur,alpha);
Cells_ActivityChange = findSignificantChanges_MT(av_raster_poststim,pstats,alpha);
[Ind_CellSelection criterion criterion_short] = findIndCellSelection_MT(Cells_ActivityChange,RedCells.isCell_red,stimInfo.index,av_raster_poststim);



%% figure with each activity graph separately

 for ii=2:6 % here for separate figures
    figure
 
%  micronsPerPixel=exptInfo.recInfo.micronsPerPixel(1);
%  imagesc(brighten(brighten(dat.mimg(:,:,2),1),0.5))
colormap gray
hold on

 for kk = 1:size(spatialInfo.ROIs,2) %all cells grey outline
     h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',1,'edgecolor',[0.7 0.7 0.7]);
 end
  
% define colour of filling
%%%% red cells
if ii==2;
    c=[1 1 1]; %white filling
%%%% yellow if sound+
elseif ii==3;
    c=[1 1 0]; %yellow
%%%% green if sound-
elseif ii==4;
    c=[0 1 0]; %green
%%%% red if laser+
elseif ii==5;
    c=[1 0 0]; %red   
%%%% blue if laser-
elseif ii==6;
    c=[0 0 1]; %green
end
    
a=Ind_CellSelection{ii}; % traces 1 and 2 are for laser    
 for kk=1:length(a) 
     h = patch(spatialInfo.ROIs{a(kk)}(:,2),spatialInfo.ROIs{a(kk)}(:,1),'y','facealpha',0.6,'edgealpha',0,'LineWidth',1,'facecolor',c);
 end

  if exist('RedCells')==1 % red cells, red outline
 a=find(RedCells.isCell_red);     
     for kk=1:length(a)
         h = patch(spatialInfo.ROIs{a(kk)}(:,2),spatialInfo.ROIs{a(kk)}(:,1),'r','facealpha',0,'edgealpha',1,'LineWidth',2,'edgecolor',[1 0 0]);
     end
  end

axis tight equal
set(gca,'TickDir','out','FontSize',14)
xlabel('Dorsal-ventral (pixels)')
ylabel('Rostral-caudal (pixels)')

      title({[mouseline,' - ',num2str(Filename(1:5)),'-',num2str(Filename(7:14)),'-',num2str(Filename(16:end-4))], ...
         ['n = ',num2str(size(raster_ordered,3)),' total cells'], ...
         [criterion_short{ii},'; n = ',num2str(length(Ind_CellSelection{ii})), '/', num2str(size(raster_ordered,3))], ...
         ['Red cells = red : n = ',num2str(length(find(RedCells.isCell_red))), '/', num2str(size(raster_ordered,3))] }) 
 end % for separate figures


%% Calculate overlap between the activity graphs and subgraphs

% I'm going to make a bar version of a venn diagram

figure(60)
for i=2:6
    subplot(5,1,i-1),
    for j=1:size(subset,1)
    sl(j)=length(setdiff(find(subset(j,3:end)),Ind_CellSelection{i}))/length(find(subset(j,3:end)));
    sa(j)=length(intersect(find(subset(j,3:end)),Ind_CellSelection{i}))/length(Ind_CellSelection{i});
    al(j)=length(setdiff(Ind_CellSelection{i},find(subset(j,3:end))))/length(Ind_CellSelection{i});
    end
    bar([sl' sa' al'])
    title(criterion_short{i})
    xlabel('subgraph number')
end

figure(61)
for j=1:size(subset,1)
    subplot(5,2,j),
    clear sl sa al
    for i=2:6
    sl(i-1)=length(setdiff(find(subset(j,3:end)),Ind_CellSelection{i}))/length(find(subset(j,3:end)));
    sa(i-1)=length(intersect(find(subset(j,3:end)),Ind_CellSelection{i}))/length(Ind_CellSelection{i});
    al(i-1)=length(setdiff(Ind_CellSelection{i},find(subset(j,3:end))))/length(Ind_CellSelection{i});
    end
    bar([sl' sa' al'])
    title(['Subgraph ',num2str(j)])
    xlabel('activity')
    xticklabels({criterion_short{2:6}})
end

%% Calculate dynamic connectivity using MTD

tic
% window of 3 for MTD / dynamic connectivity calculation

%smooth the data : 
wind_sm=5;
for i=1:size(data,1);
data_sm(i,:)=smooth(data(i,:),wind_sm);
end
toc 
tic
window = 3;
mtd=coupling_withNaNs_MT(data_sm(3:end,:)',window);
toc

% tic
% make dynamic connectivity a 2d matrix
% ind=1;
% for i=1:size(mtd,1)
%     for j=i:size(mtd,1)   
%     mtd_vec(ind,:)=mtd(i,j,:);
%     ind=ind+1;
%     end
% end
% toc
% a lot faster than code just above !!
tic
m  = (1:size(mtd,1)).' >= (1:size(mtd,2));
for i=1:size(mtd,3)
    a=mtd(:,:,i);
    mtd_vec(:,i)=a(m);
end
toc


% % make another dynamic connectivity a 2d matrix with all elements
% tic
% ind=1;
% for i=1:size(mtd,1)
%     for j=1:size(mtd,1)   
%     mtd_vec2(ind,:)=mtd(i,j,:);
%     ind=ind+1;
%     end
% end
% toc
% something faster for mtd_vec2???
tic
m  = true(size(mtd,1),size(mtd,2));
for i=1:size(mtd,3)
    a=mtd(:,:,i);
    mtd_vec2(:,i)=a(m);
end
toc


%% make rasters with coeff, with dynamic connectivity

% these parameter values are already loaded

% fr=exptInfo.fr; % frame rate
% 
% %  make the raster, eg cut up time series in each trial
% 
% rasterinfo.preOnsetTime = 1; %1 second
% rasterinfo.postOnsetTime = round((events.eventsOn(2)-events.eventsOn(1))/fr)-1; % duration between two stimuli
% rasterinfo.doDFoverF0 = 1; % =1 --> subtract the mean and divide by std of baseline

% floor in postOnsetTime*fr because otherwise !st recording is one timeframe
% longer than the rest
% coeff_raster = makeCaRaster_NaN_MT(coeff,events.eventsOn,round(rasterInfo.preOnsetTime*fr),floor(rasterInfo.postOnsetTime*fr),rasterInfo.doDFoverF0);
% mtd_vec_raster = makeCaRaster_NaN_MT(mtd_vec,events.eventsOn,round(rasterInfo.preOnsetTime*fr),floor(rasterInfo.postOnsetTime*fr),rasterInfo.doDFoverF0);
% mtd_vec2_raster = makeCaRaster_NaN_MT(mtd_vec2,events.eventsOn,round(rasterInfo.preOnsetTime*fr),floor(rasterInfo.postOnsetTime*fr),rasterInfo.doDFoverF0);

% DF0OverDo = 0
coeff_raster = makeCaRaster_NaN_MT(coeff,events.eventsOn,round(rasterInfo.preOnsetTime*fr),floor(rasterInfo.postOnsetTime*fr),0);
mtd_vec_raster = makeCaRaster_NaN_MT(mtd_vec,events.eventsOn,round(rasterInfo.preOnsetTime*fr),floor(rasterInfo.postOnsetTime*fr),0);
mtd_vec2_raster = makeCaRaster_NaN_MT(mtd_vec2,events.eventsOn,round(rasterInfo.preOnsetTime*fr),floor(rasterInfo.postOnsetTime*fr),0);


%% order the trials in raster by stimulus condition and compute average and std

[coeff_rastord mean_coeff std_mean_coeff]=makeOrdRasterAndAveTrace_MT(coeff_raster,stimInfo.index,stimInfo.order,stimInfo.repeats);
[mtd_vec_rastord mean_mtd_vec std_mean_mtd_vec]=makeOrdRasterAndAveTrace_MT(mtd_vec_raster,stimInfo.index,stimInfo.order,stimInfo.repeats);
[mtd_vec2_rastord mean_mtd_vec2 std_mean_mtd_vec2]=makeOrdRasterAndAveTrace_MT(mtd_vec2_raster,stimInfo.index,stimInfo.order,stimInfo.repeats);



%% plot mean coeff for each subgraph, with all sound amplitudes, averaged across sound amplitudes

for i=1:size(subset,1)
figure(69+i)
for j=1:size(mean_coeff,1)
    subplot(length(stimInfo.intensity),length(stimInfo.amplitude_opto),j)
    plot([1:size(mean_coeff,2)]/fr,mean_coeff(j,:,i),'k')
    xlabel('Time (s)')
    ylabel('coeff')
title({['subgraph ',num2str(i)],...
    ['Sound ',num2str(stimInfo.index(j,1)),'dB; Laser ',num2str(stimInfo.index(j,2)*100),'%']})
end
end

% average across sound amplitudes

for i=1:size(subset,1)
figure(79+i)
for j=1:length(stimInfo.amplitude_opto)
    subplot(1,length(stimInfo.amplitude_opto),j)
    plot([1:size(mean_coeff,2)]/fr,mean(mean_coeff(j:length(stimInfo.amplitude_opto):end,:,i)),'k')
    xlabel('Time (s)')
    ylabel('coeff')
title({['subgraph ',num2str(i)],...
    ['All sounds; Laser ',num2str(stimInfo.index(j,2)*100),'%']})
end
end


%% average dynamic connectivity over edges within/outside/across each subgraph

a=1:size(subset,1);
for i=1:size(subset,1);
figure(29+i)
for j=1:length(stimInfo.amplitude_opto)
    
    % connectivity within subgraph
    subplot(3,length(stimInfo.amplitude_opto),j)
    subind=find(subset(i,3:end));
    b=[];
    for k1=1:length(subind)-1
        for k2=k1+1:length(subind)
            % formula I calculated is for aij -> (i-1)*N + j for
            % mtd_vec2
        b=[b (subind(k1)-1)*size(mtd,1)+subind(k2)]; 
        end
    end
    plot([1:size(mean_coeff,2)]/fr,squeeze(mean(mean(mean_mtd_vec2(j:length(stimInfo.amplitude_opto):end,:,b),1),3)),'k')
    xlabel('Time (s)')
    ylabel('dynamic connectivity')
title({['subgraph ',num2str(i)],...
    ['Average connectivity inside subgraph'], ...
    ['All sounds; Laser ',num2str(stimInfo.index(j,2)*100),'%']})

    % connectivity outside of subgraph
    subplot(3,length(stimInfo.amplitude_opto),length(stimInfo.amplitude_opto)+j)
    subind=setdiff([1:size(mtd,1)],find(subset(i,3:end)));
    b=[];
    for k1=1:length(subind)-1
        for k2=k1+1:length(subind)
            % formula I calculated is for aij -> (i-1)*N + j for
            % mtd_vec2
        b=[b (subind(k1)-1)*size(mtd,1)+subind(k2)]; 
        end
    end
    plot([1:size(mean_coeff,2)]/fr,squeeze(mean(mean(mean_mtd_vec2(j:length(stimInfo.amplitude_opto):end,:,b),1),3)),'k')
    xlabel('Time (s)')
    ylabel('dynamic connectivity')
title({['subgraph ',num2str(i)],...
    ['Average connectivity outside subgraph'], ...
    ['All sounds; Laser ',num2str(stimInfo.index(j,2)*100),'%']})

    % connectivity between inside and outside of subgraph
    subplot(3,length(stimInfo.amplitude_opto),2*length(stimInfo.amplitude_opto)+j)
    subind1=find(subset(i,3:end));
    subind2=setdiff([1:size(mtd,1)],find(subset(i,3:end)));
    b=[];
    for k1=1:length(subind1)
        for k2=1:length(subind2)
            % formula I calculated is for aij -> (i-1)*N + j for
            % mtd_vec2
        b=[b (subind1(k1)-1)*size(mtd,1)+subind2(k2)]; 
        end
    end
    plot([1:size(mean_coeff,2)]/fr,squeeze(mean(mean(mean_mtd_vec2(j:length(stimInfo.amplitude_opto):end,:,b),1),3)),'k')
    xlabel('Time (s)')
    ylabel('dynamic connectivity')
title({['subgraph ',num2str(i)],...
    ['Average connectivity across boundary of subgraph'], ...
    ['All sounds; Laser ',num2str(stimInfo.index(j,2)*100),'%']})
end
end