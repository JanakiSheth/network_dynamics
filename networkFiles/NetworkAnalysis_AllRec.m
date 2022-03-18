%% test code for ALL PosNeg recordings / per trial - NMF analysis + incorporating other info
% to be extended to more recordings

% loading coeff and subsets for the different laser powers, then finding
% the signif ones!

clear variables
close all

% have default axes color black and not grey
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

% Select the file to open - only present subset files

mouseline='VIP';
Pathname=['/home/janaki/Dropbox/project_with_melanie/DataForJanaki/networkFiles/',mouseline,'_NMFDatasets'];

cd(Pathname)
allfiles=dir('*_subset.mat');

%% find out how many recordings and which laser conditions are available

% find unique recordings
recname=cell(size(allfiles,1),1);
for i=1:size(allfiles,1)
    recname{i}=allfiles(i).name(6:22);
end
recname=unique(recname);

nrec=size(recname,1); % number of recordings

%% find laser conditions and put both name and laser conditions in a structure

recstruct = struct([]);
for j=1:nrec
    recstruct(j,1).recname=recname{j};
    files=dir(['*',recstruct(j).recname,'*_subset.mat']);
    laserconditions=[];
    for i=1:length(files)
        ind=strfind(files(i).name,'laser');
        laserconditions=[laserconditions, str2num(files(i).name(ind+5))];
    end
    recstruct(j,1).laserconditions=laserconditions;
end

%% figure out class of functional connectivity data (double or uint8)
% not modified for loading all data from all recordings

folderpath=['/home/janaki/Dropbox/project_with_melanie/DataForJanaki/networkFiles/',mouseline,'_connectivityData/'];
cd(folderpath)

for j=1:nrec
    a=strfind({allfiles.name},recstruct(j).recname);
    b=find(~cellfun(@isempty,a));
    
    disp(allfiles(b(1)).name(6:end-19))
    connectivityfile = dir(['data_preproc_NMF_',allfiles(b(1)).name(6:end-19),'*.mat']);
    connectivity_data=load(connectivityfile.name);
    recstruct(j).class=class(connectivity_data.data);
    clear connectivity_data
end

%% load the subsets and coeffs for the laser conditions available

% moment of pause, because do I really want coeff and subset inside a
% structure or do I prefer coeff=cell(#rec,3) ?

cd(Pathname)

for j=1:nrec
    recstruct(j,1).coeff=cell(3,1);
    recstruct(j,1).subset=cell(3,1);
    files=dir(['*',recstruct(j).recname,'*_subset.mat']);
    for k=1:length(recstruct(j,1).laserconditions)
        recstruct(j,1).subset{laserconditions(k)}=load(files(k).name);
        ind=strfind(files(k).name,'_subset.mat');
        filename=[files(k).name(1:ind),'coeff.mat'];
        recstruct(j,1).coeff{laserconditions(k)}=load(filename);  
        if strcmp(recstruct(j,1).class,'uint8')
           recstruct(j,1).coeff{laserconditions(k)}.coeff=recstruct(j,1).coeff{laserconditions(k)}.coeff./255; 
        end
    end
end

%% load activity data

if mouseline == 'SOM'
    responseFunctionMat = load('/home/janaki/Dropbox/project_with_melanie/DataForJanaki/Across_many_datasets/som/TrialAverageActivityAllCellsPopulationTdelay_Fitting.mat');
else
    responseFunctionMat = load('/home/janaki/Dropbox/project_with_melanie/DataForJanaki/Across_many_datasets/vip/TrialAverageActivityAllCellsPopulationTdelay_Fitting.mat');
end
num_datasets = size(responseFunctionMat.Ind_CellSelectionIncludingSEM,1);

for j=1:nrec
    folderpath=['/home/janaki/Dropbox/project_with_melanie/DataForJanaki/networkFiles/',mouseline,'_ActivityData/'];
    cd(folderpath)
    activityfile = dir(['workspace_',recname{j},'*.mat']);
    disp(['workspace_',recname{j},'*.mat'])
    recstruct(j,1).activityworkspace=load(activityfile.name);
    recstruct(j,1).response = responseFunctionMat.Ind_CellSelectionIncludingSEM{j,:};
end

cd(Pathname) 

%% for PosNeg dataset, separate pos and neg component of coeff and calculate
% relative expression
for j=1:nrec
    for k=recstruct(j).laserconditions
    recstruct(j).coeff{k}.coeff_all=recstruct(j).coeff{k}.coeff;
    recstruct(j).coeff{k}.coeff_pos=recstruct(j).coeff{k}.coeff_all(:,1:end/2);
    recstruct(j).coeff{k}.coeff_neg=recstruct(j).coeff{k}.coeff_all(:,end/2+1:end);
    recstruct(j).coeff{k}.coeff=recstruct(j).coeff{k}.coeff_pos-recstruct(j).coeff{k}.coeff_neg;
    recstruct(j).coeff{k}.mean_relativecoeff=mean(recstruct(j).coeff{k}.coeff,2);
    end
end

%% plot, per laser trial, per recording, the mean relative coeff over all subsets

Grand_meanrelcoeff = zeros(nrec,3);
Grand_meannormrelcoeff = zeros(nrec,3);

for j=1:nrec
    for k=recstruct(j).laserconditions
        % mean relative coff over all subsets, per laser condition k, per recording j
        Grand_meanrelcoeff(j,k)=mean(recstruct(j).coeff{k}.mean_relativecoeff);
        % normalised by 1/nedges
        Grand_meannormrelcoeff(j,k)=mean(recstruct(j).coeff{k}.mean_relativecoeff*size(recstruct(j).subset{k}.subset,1)/size(recstruct(j).subset{k}.subset,2));
    end
end

cs=jet(15); % other option: hsv
set(groot,'defaultAxesColorOrder','remove');
set(groot,'defaultAxesColorOrder',cs);

figure;
subplot(121); hold on
bar(mean(Grand_meanrelcoeff),'LineWidth',1,'FaceAlpha',0.2)
plot(Grand_meanrelcoeff','LineWidth',1)
title('No normalisation')
legend({'mean', recstruct.recname})
legend('hide')
subplot(122); hold on
bar(mean(Grand_meannormrelcoeff),'LineWidth',1,'FaceAlpha',0.2)
plot(Grand_meannormrelcoeff','LineWidth',1)
legend({'mean', recstruct.recname})
legend('hide')
title({'Normalisation by', 'number of edges', 'and number of subsets'})

%% plot distribution of relative coeff weights of all subsets per laser trial

Grand_allsubs_meanrelcoeff=cell(1,3);

for j=1:nrec
    for k=recstruct(j).laserconditions
        % mean relative coff for all subsets of every recording, per laser condition k
        Grand_allsubs_meanrelcoeff{k}=[Grand_allsubs_meanrelcoeff{k}; ...
            recstruct(j).coeff{k}.mean_relativecoeff*size(recstruct(j).subset{k}.subset,1)/size(recstruct(j).subset{k}.subset,2)];  
    end
end

max_meanrelcoeff = zeros(3,1);
for i=1:3 % laser conditions
    max_meanrelcoeff(i)=max(abs(Grand_allsubs_meanrelcoeff{i}));   
end

binedges=-max(max_meanrelcoeff):max(max_meanrelcoeff)/20:max(max_meanrelcoeff);

figure
for i=1:3 % laser conditions
    subplot(2,3,i)
    histogram(Grand_allsubs_meanrelcoeff{i},binedges,'Normalization','probability')
    title({[mouseline,' - ',num2str(nrec),' recordings '], ['Laser ',num2str(i),' - mean coeff distribution'], ... 
    ['mean = ',num2str(mean(Grand_allsubs_meanrelcoeff{i}))]})
    subplot(2,3,3+i)
    boxplot(Grand_allsubs_meanrelcoeff{i},{['Laser ',num2str(i)]})
end

%% plot how relative coeff weights over all subsets change with sound level and light conditions

for j=1:nrec
    % create raster
    stimInfo=recstruct(j).activityworkspace.stimInfo;
    for k=recstruct(j).laserconditions
        for i=1:size(recstruct(j).coeff{k}.coeff,2)/181
            coeff_rast(i,1:181,:)=recstruct(j).coeff{k}.coeff(:,181*(i-1)+1:181*i)'; 
        end
        index=stimInfo.index(find(stimInfo.index(:,2)==stimInfo.amplitude_opto(k)),:);
        order=1+floor((stimInfo.order(ismember(stimInfo.order,find(stimInfo.index(:,2)==stimInfo.amplitude_opto(k))))-1)/3);
        % save ordered raster, mean and std
        [coeff_rastord, mean_coeff, std_mean_coeff]=makeOrdRasterAndAveTrace_MT(coeff_rast,index,order,stimInfo.repeats);
        
        recstruct(j).coeff{k}.coeff_rastord=coeff_rastord;
        recstruct(j).coeff{k}.mean_coeff=mean_coeff;
        recstruct(j).coeff{k}.std_mean_coeff=std_mean_coeff;
        clear coeff_rast coeff_rastord mean_coeff std_mean_coeff;
    end       
end
repeats = 10;
meanCoeff_SoundAmplitudeVariation = cell(3,1);
coeffRastord_SoundAmplitudeVariation = cell(3,1);
soundResponsivenessOfCoeff = cell(3,4);

for k=1:3
    for j = 1:nrec
        number_subsets = size(recstruct(j).subset{k}.subset,1);
        for sub = 1:number_subsets
            meanCoeff_SoundAmplitudeVariation{k} = [meanCoeff_SoundAmplitudeVariation{k}, ...
                mean(recstruct(j).coeff{k}.mean_coeff(:,:,sub),2)*size(recstruct(j).subset{k}.subset,1)/size(recstruct(j).subset{k}.subset,2)];
            raster_ordered = mean(recstruct(j).coeff{k}.coeff_rastord(:,:,sub),2)*size(recstruct(j).subset{k}.subset,1)/size(recstruct(j).subset{k}.subset,2);
            coeffRastord_SoundAmplitudeVariation{k} = [coeffRastord_SoundAmplitudeVariation{k}, raster_ordered];
            isCoeffSoundResponsive = zeros(6,1);
            for amp = 2:7
                isCoeffSoundResponsive(amp) = ttest(raster_ordered(1:repeats),raster_ordered((amp-1)*10+1:10*(amp-1)+repeats),'Alpha',0.05);
            end
            if sum(isCoeffSoundResponsive)>0
                soundResponsivenessOfCoeff{k,1} = [soundResponsivenessOfCoeff{k,1},j];
                soundResponsivenessOfCoeff{k,2} = [soundResponsivenessOfCoeff{k,2},sub];
                %figure;
                %histogram(recstruct(j).subset{k}.subset(sub,:),20);
                %plot(1:7,mean(recstruct(j).coeff{k}.mean_coeff(:,:,sub),2)*size(recstruct(j).subset{k}.subset,1)/size(recstruct(j).subset{k}.subset,2));
            end
        end
    end
end

figure; hold on;
errorbar(mean(meanCoeff_SoundAmplitudeVariation{1},2),std(meanCoeff_SoundAmplitudeVariation{1},0,2)/sqrt(length(meanCoeff_SoundAmplitudeVariation{1})),'k')
errorbar(mean(meanCoeff_SoundAmplitudeVariation{2},2),std(meanCoeff_SoundAmplitudeVariation{2},0,2)/sqrt(length(meanCoeff_SoundAmplitudeVariation{2})),'b')
errorbar(mean(meanCoeff_SoundAmplitudeVariation{3},2),std(meanCoeff_SoundAmplitudeVariation{3},0,2)/sqrt(length(meanCoeff_SoundAmplitudeVariation{3})),'r')
xlabel('Sound Amplitude');
xticklabels([0,30,50,60,70,80,90])
ylabel('Average of mean coeff across subsets and recordings')

%% plot distribution of number of subsets per laser and the clustering weights of the cells in a subset vs the monotonicity of the cells
% the latter part asks the question of - Is monotonicity responsible for the
% subset's expressiveness in the presence of sound and light?

for j=1:nrec    
    for k=recstruct(j).laserconditions
        recstruct(j).nsubsets(k)=size(recstruct(j).subset{k}.subset,1);
    end
end
Grand_nsubsets=reshape([recstruct.nsubsets],3,nrec);

figure
for i=1:3 % laser conditions
    subplot(1,3,i)
    histogram(Grand_nsubsets(i,:),1.5:1:15.5)
    title({[mouseline,' - ',num2str(nrec),' recordings '], ...
        ['Laser ',num2str(i),' - number of subsets'], ...
        ['mean + SEM = ',num2str(round(mean(Grand_nsubsets(i,:)),1)),' + ',num2str(round(std(Grand_nsubsets(i,:))/sqrt(nrec),1))]})
    xlabel('number of subsets')
    ylabel('count')
end

% First step is to rearrange subset in matrix form
%copied from other mat file
% fill in matrix for every subset - normalise with L2 norm

for j=1:nrec
    recstruct(j).ncells= 1/2*(1+sqrt(1+8*size(recstruct(j).subset{recstruct(j).laserconditions(1)}.subset,2)));
    recstruct(j).nedges=size(recstruct(j).subset{laserconditions(1)}.subset,2);
    for k=recstruct(j).laserconditions
        for i=1:size(recstruct(j).subset{k}.subset,1) % number of subsets
            recstruct(j).subset{k}.subset_normL2(i)=norm(recstruct(j).subset{k}.subset(i,:)); % named subset_matnorm in other analysis file
            B=tril(ones(recstruct(j).ncells),-1);
            B(B==1)=recstruct(j).subset{k}.subset(i,:);
            recstruct(j).subset{k}.subset_matnonorm(i,:,:)=B+B';
            recstruct(j).subset{k}.subset_matnormL2(i,:,:)=recstruct(j).subset{k}.subset_matnonorm(i,:,:)/recstruct(j).subset{k}.subset_normL2(i);
        end
    end
end

for j=1:nrec  
    [miNoLight, miMidLight, miHighLight] = load_ResponseData(recstruct(j).response);
    recstruct(j).subset{1}.monotonicity_index = miNoLight;
    recstruct(j).subset{2}.monotonicity_index = miMidLight;
    recstruct(j).subset{3}.monotonicity_index = miHighLight;
    for k=recstruct(j).laserconditions    
        recstruct(j).subset{k}.clustering_wei=zeros(recstruct(j).ncells,size(recstruct(j).subset{k}.subset,1)); % weighted clustering for all subsets of each laser condition  
        for sub = 1:size(recstruct(j).subset{k}.subset,1)
            a=squeeze(recstruct(j).subset{k}.subset_matnonorm(sub,:,:));            
            % clustering_coeff function gives the clustering per node
            b=clustering_coef_wu(a);
            recstruct(j).subset{k}.clustering_wei(:,sub)=b;
        end
    end
end

for k = 1:3
    for combinations = 1:length(soundResponsivenessOfCoeff{k,1})
        j = soundResponsivenessOfCoeff{k,1}(combinations);
        sub = soundResponsivenessOfCoeff{k,2}(combinations);
        %figure; hold on;
        %plot(recstruct(j).subset{k}.monotonicity_index,recstruct(j).subset{k}.clustering_wei(:,sub),'.')
        %plot(recstruct(j).subset{k}.monotonicity_index(recstruct(j).activityworkspace.Ind_cells_red), ...
        %recstruct(j).subset{k}.clustering_wei(recstruct(j).activityworkspace.Ind_cells_red,sub),'r*')
        [rho,pval] = corr(recstruct(j).subset{k}.monotonicity_index,recstruct(j).subset{k}.clustering_wei(:,sub),'Type','Spearman');
        soundResponsivenessOfCoeff{k,3} = [soundResponsivenessOfCoeff{k,3},rho];
        soundResponsivenessOfCoeff{k,4} = [soundResponsivenessOfCoeff{k,4},pval];
        if pval<0.05
            figure;
            plot(1:7,mean(recstruct(j).coeff{k}.mean_coeff(:,:,sub),2)*size(recstruct(j).subset{k}.subset,1)/size(recstruct(j).subset{k}.subset,2));
            title(['rho ',num2str(rho),',pval ',num2str(pval)])
        end
    end
end


%% plot how subset weights (mean  +/- sem) change with laser conditions based on if the subset coeff means are positive or negative.

for k=1:3
    figure(); hold on;
    xlim([0,3]);
    for j = 1:nrec
        for sub = 1:length(recstruct(j).coeff{k}.mean_relativecoeff)
            mean_subset_weight = mean(recstruct(j).subset{k}.subset(sub,:));
            sem_subset_weight = std(recstruct(j).subset{k}.subset(sub,:))/sqrt(length(recstruct(j).subset{k}.subset(sub,:)));
            if recstruct(j).coeff{k}.mean_relativecoeff(sub) < 0
                errorbar(1,mean_subset_weight,sem_subset_weight);
            else
                errorbar(2,mean_subset_weight,sem_subset_weight);
            end
        end
    end
end

%% average nodes and edges per threshold over recordings. For each recording the average is computed over all subsets for each laser condition. 

%%% questions to solve: is thres/nedges the correct scale or should I
%%% normalise by norm2 the subset weights?
%%% questions to solve: add error bars

thres_weight=[0, logspace(-4,2,1000)];
for j=1:nrec
    recstruct(j).ncells= 1/2*(1+sqrt(1+8*size(recstruct(j).subset{recstruct(j).laserconditions(1)}.subset,2)));
    recstruct(j).nedges=size(recstruct(j).subset{laserconditions(1)}.subset,2);

    % for k=recstruct(j).laserconditions
    %     recstruct(j).nsubsets(k)=size(recstruct(j).subset{k}.subset,1);
    % end

    recstruct(j).ratioedgesperthres_allsubsperlaser=zeros(length(thres_weight),3); % total number of edges above threshold for each laser condition
    recstruct(j).rationodesperthres_allsubsperlaser=zeros(length(thres_weight),3); % total number of nodes above threshold for each laser condition

    % go back to tracking index for edges
    ind=1;
    ind_track=[];
    for i=1:recstruct(j).ncells
        for j2=i+1:recstruct(j).ncells
            ind_track(ind,:)=[i, j2];
            ind=ind+1;
        end
    end
    % tricky thing: not the same number of subsets - so number of edges at zero
    % threshold won't be the same

    for k=recstruct(j).laserconditions
        for i=1:length(thres_weight)
            a=recstruct(j).subset{k}.subset;
            a(a<thres_weight(i)/recstruct(j).nedges)=0;
            recstruct(j).ratioedgesperthres_allsubsperlaser(i,k)=length(find(recstruct(j).subset{k}.subset>=thres_weight(i)/recstruct(j).nedges))/...
                (recstruct(j).nsubsets(k)*recstruct(j).nedges); % this computes the proportion of significant edges in each subset. significance depends upon the threshold.
            % find index of non zero edges
            [~,c]=find(a);  
            recstruct(j).rationodesperthres_allsubsperlaser(i,k)=length(unique(ind_track(c,:)))/recstruct(j).ncells;   
        end
    end
end

A = zeros(nrec,length(thres_weight),3);
B = zeros(nrec,length(thres_weight),3);
for j=1:nrec
    A(j,:,:)=recstruct(j).ratioedgesperthres_allsubsperlaser;
    B(j,:,:)=recstruct(j).rationodesperthres_allsubsperlaser;
end

Grand_ratioedgesperthres=squeeze(mean(A,1));
Grand_ratioedgesperthres_sem=squeeze(std(A,0,1)./sqrt(nrec));
Grand_rationodesperthres=squeeze(mean(B,1));
Grand_rationodesperthres_sem=squeeze(std(B,0,1)./sqrt(nrec));

 %% plot
cm=[0 0 0;
    1 0.7 0.7;
    1 0 0];

set(groot,'defaultAxesColorOrder',cm)

figure;
subplot(121)
loglog(thres_weight,Grand_ratioedgesperthres,'LineWidth',2)
loglog(thres_weight,Grand_ratioedgesperthres+Grand_ratioedgesperthres_sem,'LineWidth',0.5,'LineStyle','--')
loglog(thres_weight,Grand_ratioedgesperthres-Grand_ratioedgesperthres_sem,'LineWidth',0.5,'LineStyle','--')
title({[mouseline,' - ',num2str(nrec),' recordings '], 'Ratio of edges above threshold',...
'(across subsets) per laser condition'})
xlabel('cutoff = threshold / nedges')
ylabel('Ratio of edges above threshold + SEM')
ylim([1e-6 1])
legend('laser 1', 'laser 2', 'laser 3')
subplot(122), hold on;
semilogx(thres_weight,Grand_rationodesperthres,'LineWidth',2)
semilogx(thres_weight,Grand_rationodesperthres+Grand_rationodesperthres_sem,'LineWidth',0.5,'LineStyle','--')
semilogx(thres_weight,Grand_rationodesperthres-Grand_rationodesperthres_sem,'LineWidth',0.5,'LineStyle','--')
title({[mouseline,' - ',num2str(nrec),' recordings '], 'Ratio of nodes above threshold',...
'(across subsets) per laser condition'})
xlabel('cutoff = threshold / nedges')
ylabel('Ratio of nodes above threshold + SEM')
legend('laser 1', 'laser 2', 'laser 3')

%% same figure as before but with only edges between the Sound+ cells (Average nodes and edges per threshold over recordings. For each recording the average is computed over all subsets for each laser condition.)
% do it - for all subsets and  - for most prominent subset at laser3

for j=1:nrec

    ind_SoundP=recstruct(j).activityworkspace.Ind_CellSelection{3};

    recstruct(j).ratioedgesperthres_SoundP_allsubsperlaser=zeros(length(thres_weight),3); % total number of edges above threshold for each laser condition
    recstruct(j).rationodesperthres_SoundP_allsubsperlaser=zeros(length(thres_weight),3); % total number of nodes above threshold for each laser condition

    % go back to tracking index for edges
    ind=1;
    ind_track=[];
    ind_track_SoundP=[];
    for i=1:recstruct(j).ncells
        for j2=i+1:recstruct(j).ncells
            ind_track(ind,:)=[i, j2];
            if ismember(i,ind_SoundP) && ismember(j2,ind_SoundP)
                ind_track_SoundP=[ind_track_SoundP, ind]; %this variable stores the indices corresponding to cell pairs where both cells are sound positive. 
            end   
            ind=ind+1;
        end
    end

    % tricky thing: not the same number of subsets - so number of edges at zero
    % threshold won't be the same

    for k=recstruct(j).laserconditions
        for i=1:length(thres_weight)
            a=recstruct(j).subset{k}.subset(:,ind_track_SoundP);
            a(a<thres_weight(i)/recstruct(j).nedges)=0;
            recstruct(j).ratioedgesperthres_SoundP_allsubsperlaser(i,k)=length(find(recstruct(j).subset{k}.subset(:,ind_track_SoundP)>=thres_weight(i)/recstruct(j).nedges))/...
                (recstruct(j).nsubsets(k)*length(ind_track_SoundP));
            % find index of non zero edges
            [~,c]=find(a);  
            recstruct(j).rationodesperthres_SoundP_allsubsperlaser(i,k)=length(unique(ind_track(ind_track_SoundP(c),:)))/length(ind_SoundP);   
        end
    end
end

clear A B
A = zeros(nrec,length(thres_weight),3);
B = zeros(nrec,length(thres_weight),3);
for j=1:nrec
    A(j,:,:)=recstruct(j).ratioedgesperthres_SoundP_allsubsperlaser;
    B(j,:,:)=recstruct(j).rationodesperthres_SoundP_allsubsperlaser;
end

Grand_SoundP_ratioedgesperthres=squeeze(mean(A,1));
Grand_SoundP_ratioedgesperthres_sem=squeeze(std(A,0,1)./sqrt(nrec));
Grand_SoundP_rationodesperthres=squeeze(mean(B,1));
Grand_SoundP_rationodesperthres_sem=squeeze(std(B,0,1)./sqrt(nrec));

%% plot

% cm is laser color scheme
cm=[0 0 0;
    1 0.7 0.7;
    1 0 0];

set(groot,'defaultAxesColorOrder',cm)

figure;
subplot(121), hold on;
loglog(thres_weight,Grand_SoundP_ratioedgesperthres,'LineWidth',2)
loglog(thres_weight,Grand_SoundP_ratioedgesperthres+Grand_SoundP_ratioedgesperthres_sem,'LineWidth',0.5,'LineStyle','--')
loglog(thres_weight,max(Grand_SoundP_ratioedgesperthres-Grand_SoundP_ratioedgesperthres_sem,0),'LineWidth',0.5,'LineStyle','--')
title({[mouseline,' - ',num2str(nrec),' recordings '], ...
    'Ratio of Sound+ edges above threshold','(across subsets) per laser condition'})
xlabel('cutoff = threshold / nedges')
ylabel('Ratio of Sound+ edges above threshold')
ylim([1e-6 1])
legend('laser 1', 'laser 2', 'laser 3')

subplot(122), hold on;
semilogx(thres_weight,Grand_SoundP_rationodesperthres,'LineWidth',2)
semilogx(thres_weight,Grand_SoundP_rationodesperthres+Grand_SoundP_rationodesperthres_sem,'LineWidth',0.5,'LineStyle','--')
semilogx(thres_weight,max(Grand_SoundP_rationodesperthres-Grand_SoundP_rationodesperthres_sem,0),'LineWidth',0.5,'LineStyle','--')
title({[mouseline,' - ',num2str(nrec),' recordings '], ...
     'Ratio of Sound+ nodes above threshold', '(across subsets) per laser condition'})
xlabel('cutoff = threshold / nedges')
ylabel('Ratio of Sound+ nodes above threshold')
legend('laser 1', 'laser 2', 'laser 3')

%% Edge distribution

% use thres_weights/nedges as bin edges
% normalise by total number of subsets per laser condition

edgecount=zeros(length(laserconditions),length(thres_weight)-1);
edgecount_SoundP=zeros(length(laserconditions),length(thres_weight)-1);

Grand_AllEdges_norm=cell(3,1);
Grand_AllEdges_SoundP_norm=cell(3,1);

for j=1:nrec
    
    ind_SoundP=recstruct(j).activityworkspace.Ind_CellSelection{3};

    % go back to tracking index for edges
    ind=1;
    ind_track=[];
    ind_track_SoundP=[];
    for i=1:recstruct(j).ncells
        for j2=i+1:recstruct(j).ncells
            ind_track(ind,:)=[i, j2];
            if ismember(i,ind_SoundP) && ismember(j2,ind_SoundP)
                ind_track_SoundP=[ind_track_SoundP, ind];
            end   
            ind=ind+1;
        end
    end
    
    for k=recstruct(j).laserconditions
        for i=1:size(recstruct(j).subset{k}.subset,1)
            [a, ~]= histcounts(recstruct(j).subset{k}.subset(i,:),thres_weight/recstruct(j).nedges);
            edgecount(k,:) = edgecount(k,:)+a;
            
            [a, ~]= histcounts(recstruct(j).subset{k}.subset(i,ind_track_SoundP),thres_weight/recstruct(j).nedges);
            edgecount_SoundP(k,:) = edgecount_SoundP(k,:)+a;           
        end    
    end
    
    % reshape total number of edges in one laser condition
    for k=recstruct(j).laserconditions
        Grand_AllEdges_norm{k}=[Grand_AllEdges_norm{k}; reshape(recstruct(j).subset{k}.subset*recstruct(j).nedges,[],1)];
        Grand_AllEdges_SoundP_norm{k}=[Grand_AllEdges_SoundP_norm{k}; reshape(recstruct(j).subset{k}.subset(:,ind_track_SoundP)*recstruct(j).nedges,[],1)];
    end   
end
% normalised by total number of edges per laser condiiton
edgecount_norm=edgecount./sum(edgecount,2);
edgecount_SoundP_norm=edgecount_SoundP./sum(edgecount_SoundP,2);

figure;
subplot(211), hold on;
for k=recstruct(j).laserconditions
    histogram('BinEdges',thres_weight(2:end),'BinCounts',edgecount_norm(k,2:end),'DisplayStyle','stairs')
end
title({[mouseline,' - ',num2str(nrec),' recordings '], 'Edge weight distribution per laser'})
xlabel('cutoff = threshold / nedges')
ylabel('count (normalised)')
legend('laser 1', 'laser 2', 'laser 3')
xlim([0, 10])

% boxplot of edges distributions: will give 25% interquartile: can see
% whether top end of weights is modulated
subplot(212)
boxplot([Grand_AllEdges_norm{1};Grand_AllEdges_norm{2};Grand_AllEdges_norm{3}],[repmat({'Laser 1'},length(Grand_AllEdges_norm{1}),1); repmat({'Laser 2'},length(Grand_AllEdges_norm{2}),1);  repmat({'Laser 3'},length(Grand_AllEdges_norm{3}),1);  ])

% edges between Sound+ cells
figure, 
subplot(311), hold on;
for k=recstruct(j).laserconditions
    histogram('BinEdges',thres_weight(2:end),'BinCounts',edgecount_SoundP_norm(k,2:end),'DisplayStyle','stairs')
end
title({[mouseline,' - ',num2str(nrec),' recordings '], ...
    'Edge weight distribution of Sound+ edges per laser'})
xlabel('cutoff = threshold / nedges')
ylabel('count (normalised)')
legend('laser 1', 'laser 2', 'laser 3')
xlim([0, 10])
% boxplot
subplot(312),
boxplot([Grand_AllEdges_SoundP_norm{1};Grand_AllEdges_SoundP_norm{2};Grand_AllEdges_SoundP_norm{3}],[repmat({'Laser 1'},length(Grand_AllEdges_SoundP_norm{1}),1); repmat({'Laser 2'},length(Grand_AllEdges_SoundP_norm{2}),1);  repmat({'Laser 3'},length(Grand_AllEdges_SoundP_norm{3}),1);  ])
title('box plot')
ylabel('edge weight * nedges')
% mean + STD
subplot(313),
a=[std(Grand_AllEdges_SoundP_norm{1})/sqrt(sum(edgecount_SoundP(1,:)));std(Grand_AllEdges_SoundP_norm{2})/sqrt(sum(edgecount_SoundP(2,:)));std(Grand_AllEdges_SoundP_norm{3})/sqrt(sum(edgecount_SoundP(3,:)))];
b=[mean(Grand_AllEdges_SoundP_norm{1});mean(Grand_AllEdges_SoundP_norm{2});mean(Grand_AllEdges_SoundP_norm{3})];
errorbar(b,a,'LineWidth',2)
title('edge weight * nedges: mean + SEM')
xlabel('Laser condition')
ylabel('mean edge weight * nedges + SEM')
xlim([0.5, 3.5])

%% Select one subset per recording, laser condition, with highest calcium expression in a given time window

% first, make rasters out of coeff 

for j=1:nrec
    % create raster
    stimInfo=recstruct(j).activityworkspace.stimInfo;
    
    for k=recstruct(j).laserconditions
        for i=1:size(recstruct(j).coeff{k}.coeff,2)/181
            coeff_rast(i,1:181,:)=recstruct(j).coeff{k}.coeff(:,181*(i-1)+1:181*i)'; 
        end
        index=stimInfo.index(find(stimInfo.index(:,2)==stimInfo.amplitude_opto(k)),:);
        order=1+floor((stimInfo.order(ismember(stimInfo.order,find(stimInfo.index(:,2)==stimInfo.amplitude_opto(k))))-1)/3);
        % save ordered raster, mean and std
        [coeff_rastord, mean_coeff, std_mean_coeff]=makeOrdRasterAndAveTrace_MT(coeff_rast,index,order,stimInfo.repeats);

        recstruct(j).coeff{k}.coeff_rastord=coeff_rastord;
        recstruct(j).coeff{k}.mean_coeff=mean_coeff;
        recstruct(j).coeff{k}.std_mean_coeff=std_mean_coeff;
        clear coeff_rast coeff_rastord mean_coeff std_mean_coeff;
    end       
end

% average over given time window and select subset with highest expression 
% during that time window

T_selectwindow=[2, 3]; % calculate mean from second 2 to 3 after beginning of trial i.e. during stimulus
for j=1:nrec
    for k=recstruct(j).laserconditions
        fr=recstruct(j).activityworkspace.fr;
        a=squeeze(mean(mean(recstruct(j).coeff{k}.coeff_rastord(:,round(fr*T_selectwindow(1)):round(fr*T_selectwindow(2)),:),2),1));
        [~,I]=sort(a,'descend'); 
        recstruct(j).SelectSubset(k)=I(1);
    end
end

%% average across recordings a bunch of things for each laser condition ....

% brainstorm :
% - all coeff traces of selected subsets of all recordings, and average
% trace: DONE
% - edges per threshold: DONE
% - nodes per threshold: DONE
% - weighted / binary clustering coefficient: DONE
% - weighted / binary clustering coefficient for the SoundP nodes: DONE
% - connectivity matrices averaged across activity groups
% - overlap between subsets at different laser powers?


% questions : I'm wondering whether subset weight/nedges is really the best
% way to do it - it centers all the curves around one which is the mean ...
% 

%% 1) plot all coeff traces of selected subsets of all recordings, and average
% coeff trace

Grand_meancoeff=zeros(nrec,7,181,length(laserconditions));
for j=1:nrec
    for k=recstruct(j).laserconditions
        Grand_meancoeff(j,:,:,k)=recstruct(j).coeff{k}.mean_coeff(:,:,recstruct(j).SelectSubset(k));
    end
end

% co is grayscale for different sound amps
co = [0.9 0.9 0.9;
    0.8 0.8 0.8;
    0.7 0.7 0.7;
    0.6 0.6 0.6;
    0.4 0.4 0.4;
    0.2 0.2 0.2;
    0 0 0];

% is for different recordings
cs=jet(15); 

figure
for k=laserconditions
    subplot(3,3,k),
    set(groot,'defaultAxesColorOrder',cs);
    plot((1:181)/30,squeeze(mean(Grand_meancoeff(:,:,:,k),2)))
    title({[mouseline,' - ',num2str(nrec),' recordings - laser ',num2str(k)'], ...
            'average coeff for the','different recordings'})
    xlabel('Time (s)')
    ylabel('mean coeff')

    subplot(3,3,3+k),
    set(groot,'defaultAxesColorOrder',co)
    plot((1:181)/30,squeeze(mean(Grand_meancoeff(:,:,:,k),1)))
    title({'average coeff for the','different sound amplitudes'})
    xlabel('Time (s)')
    ylabel('mean coeff')

    subplot(3,3,6+k),
    plot((1:181)/30,squeeze(mean(mean(Grand_meancoeff(:,:,:,k),2),1)),'k','LineWidth',2)
    title({'average coeff over', 'all recs, all sound amps'})
    xlabel('Time (s)')
    ylabel('mean coeff')
end

figure, hold on;
set(groot,'defaultAxesColorOrder',cm)
for k=laserconditions
    plot((1:181)/30,squeeze(mean(mean(Grand_meancoeff(:,:,:,k),2),1)),'LineWidth',2)
end
title({[mouseline,' - ',num2str(nrec),' recordings'], ...
    'average coeff over','all recs, all sound amps'})
xlabel('Time (s)')
ylabel('mean coeff')
legend('laser 1', 'laser 2', 'laser 3')

%% 2) edges per threshold for selected subset
% already calculated before for all subsets: recstruct(j).ratioedgesperthres_allsubsperlaser(i,k)

for j=1:nrec
    recstruct(j).ratioedgesperthres_SelectSubperlaser=zeros(length(thres_weight),3); % total number of edges above threshold for each laser condition
    recstruct(j).rationodesperthres_SelectSubperlaser=zeros(length(thres_weight),3); % total number of nodes above threshold for each laser condition

    % go back to tracking index for edges
    ind=1;
    ind_track=[];
    for i=1:recstruct(j).ncells
        for j2=i+1:recstruct(j).ncells
            ind_track(ind,:)=[i, j2];
            ind=ind+1;
        end
    end

    for k=recstruct(j).laserconditions
        for i=1:length(thres_weight)
            a=recstruct(j).subset{k}.subset(recstruct(j).SelectSubset(k),:);
            a(a<thres_weight(i)/recstruct(j).nedges)=0;
            recstruct(j).ratioedgesperthres_SelectSubperlaser(i,k)=length(find(recstruct(j).subset{k}.subset(recstruct(j).SelectSubset(k),:)>=thres_weight(i)/recstruct(j).nedges))...
                /(recstruct(j).nedges);
            % find index of non zero edges
            [~,c]=find(a);  
            recstruct(j).rationodesperthres_SelectSubperlaser(i,k)=length(unique(ind_track(c,:)))/recstruct(j).ncells;   
        end
    end
end

clear A B
A = zeros(nrec,length(thres_weight),3);
B = zeros(nrec,length(thres_weight),3);
for j=1:nrec
    A(j,:,:)=recstruct(j).ratioedgesperthres_SelectSubperlaser;
    B(j,:,:)=recstruct(j).rationodesperthres_SelectSubperlaser;
end

Grand_SelectSub_ratioedgesperthres=A;
Grand_SelectSub_rationodesperthres=B;
Grand_SelectSub_ratioedgesperthres_mean=squeeze(mean(A,1));
Grand_SelectSub_ratioedgesperthres_sem=squeeze(std(A,0,1)./sqrt(nrec));
Grand_SelectSub_rationodesperthres_mean=squeeze(mean(B,1));
Grand_SelectSub_rationodesperthres_sem=squeeze(std(B,0,1)./sqrt(nrec));

% plot it, for each recording, then mean

figure;
for k=laserconditions
    subplot(2,4,k), hold on;
    loglog(thres_weight,Grand_SelectSub_ratioedgesperthres(:,:,k),'LineWidth',2,'Color',cm(k,:))
    title({[mouseline,' - ',num2str(nrec),' recordings - laser ',num2str(k)], ...
            'Ratio of edges above threshold','for each recording'})
    xlabel('threshold / nedges')
    ylabel('Ratio of SelectSub edges above threshold')
    ylim([1e-5, 1])

    subplot(2,4,4), hold on;
    loglog(thres_weight,Grand_SelectSub_ratioedgesperthres_mean,'LineWidth',2)
    loglog(thres_weight,Grand_SelectSub_ratioedgesperthres_mean+Grand_SelectSub_ratioedgesperthres_sem,'LineWidth',0.5,'LineStyle','--')
    loglog(thres_weight,max(Grand_SelectSub_ratioedgesperthres_mean-Grand_SelectSub_ratioedgesperthres_sem,0),'LineWidth',0.5,'LineStyle','--')
    title({[mouseline,' - ',num2str(nrec),' recordings '], ...
        'Ratio of edges above threshold','from all highest subsets per laser condition'})
    xlabel('threshold / nedges')
    ylabel('Ratio of SelectSub edges above threshold')
    ylim([1e-5, 1])
    legend('laser 1', 'laser 2', 'laser 3')

    subplot(2,4,4+k), hold on;
    loglog(thres_weight,Grand_SelectSub_rationodesperthres(:,:,k),'LineWidth',2,'Color',cm(k,:))
    title({[mouseline,' - ',num2str(nrec),' recordings - laser ',num2str(k)], ...
            'Ratio of nodes above threshold','for each recording'})
    xlabel('threshold / nedges')
    ylabel('Ratio of SelectSub nodes above threshold')
    ylim([1e-2, 1])

    subplot(2,4,8), hold on;
    loglog(thres_weight,Grand_SelectSub_rationodesperthres_mean,'LineWidth',2)
    loglog(thres_weight,Grand_SelectSub_rationodesperthres_mean+Grand_SelectSub_rationodesperthres_sem,'LineWidth',0.5,'LineStyle','--')
    loglog(thres_weight,max(Grand_SelectSub_rationodesperthres_mean-Grand_SelectSub_rationodesperthres_sem,0),'LineWidth',0.5,'LineStyle','--')
    title({[mouseline,' - ',num2str(nrec),' recordings '], ...
         'Ratio of nodes above threshold','from all highest subsets per laser condition'})
    xlabel('threshold / nedges')
    ylabel('Ratio of SelectSub nodes above threshold')
    ylim([1e-3, 1])
    legend('laser 1', 'laser 2', 'laser 3')
end

%% 3) Clustering

% First step is to rearrange subset in matrix form
%copied from other mat file
% fill in matrix for every subset - normalise with L2 norm

for j=1:nrec
    for k=recstruct(j).laserconditions
        for i=1:size(recstruct(j).subset{k}.subset,1) % number of subsets
            recstruct(j).subset{k}.subset_normL2(i)=norm(recstruct(j).subset{k}.subset(i,:)); % named subset_matnorm in other analysis file
            B=tril(ones(recstruct(j).ncells),-1);
            B(B==1)=recstruct(j).subset{k}.subset(i,:);
            recstruct(j).subset{k}.subset_matnonorm(i,:,:)=B+B';
            recstruct(j).subset{k}.subset_matnormL2(i,:,:)=recstruct(j).subset{k}.subset_matnonorm(i,:,:)/recstruct(j).subset{k}.subset_normL2(i);
        end
    end
end

% Second step is to calculate the clustering on thresholded matrix musing thres_weight from previously

for j=1:nrec
    recstruct(j).clusteringperthres_SelectSub_bin=zeros(length(thres_weight),length(laserconditions)); % binary clustering for the selected subset of each laser condition
    recstruct(j).clusteringperthres_SelectSub_SoundP_bin=zeros(length(thres_weight),length(laserconditions));
    recstruct(j).clusteringperthres_SelectSub_wei=zeros(length(thres_weight),length(laserconditions)); % weighted clustering for the selected subset of each laser condition    
    recstruct(j).clusteringperthres_SelectSub_SoundP_wei=zeros(length(thres_weight),length(laserconditions));
    for k=recstruct(j).laserconditions      
        for i=1:length(thres_weight)
            a=squeeze(recstruct(j).subset{k}.subset_matnonorm(recstruct(j).SelectSubset(k),:,:));
            a(a<thres_weight(i)/recstruct(j).nedges)=0;
            
            % clustering_coeff function gives the clustering per node - average
            % over all non-zero nodes? do both?
            b=clustering_coef_wu(a);
            recstruct(j).clusteringperthres_SelectSub_wei(i,k)=mean(b);
            % mean over all SoundP cells
            recstruct(j).clusteringperthres_SelectSub_SoundP_wei(i,k)=mean(b(recstruct(j).activityworkspace.Ind_CellSelection{3}));
            
            a2=a;
            a2(a>=thres_weight(i)/recstruct(j).nedges)=1;
            b2=clustering_coef_bu(a2);
            recstruct(j).clusteringperthres_SelectSub_bin(i,k)=mean(b2);
             % mean over all SoundP cells
            recstruct(j).clusteringperthres_SelectSub_SoundP_bin(i,k)=mean(b2(recstruct(j).activityworkspace.Ind_CellSelection{3}));
        end
    end
end

% third step is to plot

% First line is binary clustering, second line is weighted clustering
figure;
for k=laserconditions   
    subplot(2,4,k), hold on;
    a=[recstruct.clusteringperthres_SelectSub_bin];
    loglog(thres_weight,a(:,k:3:end)','LineWidth',2,'Color',cm(k,:))
    title({[mouseline,' - ',num2str(nrec),' recordings - laser ',num2str(k)], ...
    'Binary clustering of highest subset','for each recording'})
    xlabel('threshold / nedges')
    ylabel('Binary clustering')
    % ylim([1e-5 1])

    subplot(2,4,4), hold on;
    a=[recstruct.clusteringperthres_SelectSub_bin];
    loglog(thres_weight,mean(a(:,k:3:end)'),'LineWidth',2)
    % loglog(thres_weight,Grand_SelectSub_ratioedgesperthres_mean+Grand_SelectSub_ratioedgesperthres_sem,'LineWidth',0.5,'LineStyle','--')
    % loglog(thres_weight,max(Grand_SelectSub_ratioedgesperthres_mean-Grand_SelectSub_ratioedgesperthres_sem,0),'LineWidth',0.5,'LineStyle','--')
    title({[mouseline,' - ',num2str(nrec),' recordings '], ...
    'Binary clustering','from all highest subset per laser condition'})
    xlabel('threshold / nedges')
    ylabel('Mean binary clustering of SelectSub ')
    % ylim([1e-5 1])
    % legend('laser 1', 'laser 2', 'laser 3')

    subplot(2,4,4+k), hold on;
    a=[recstruct.clusteringperthres_SelectSub_wei];
    loglog(thres_weight,a(:,k:3:end)','LineWidth',2,'Color',cm(k,:))
    title({[mouseline,' - ',num2str(nrec),' recordings - laser ',num2str(k)], ...
    'Weighted clustering of highest subset','for each recording'})
    xlabel('threshold / nedges')
    ylabel('Ratio of SelectSub nodes above threshold')
    % ylim([1e-2 1])

    subplot(2,4,8), hold on;
    a=[recstruct.clusteringperthres_SelectSub_wei];
    loglog(thres_weight,mean(a(:,k:3:end)'),'LineWidth',2)
    % loglog(thres_weight,Grand_SelectSub_rationodesperthres_mean+Grand_SelectSub_rationodesperthres_sem,'LineWidth',0.5,'LineStyle','--')
    % loglog(thres_weight,max(Grand_SelectSub_rationodesperthres_mean-Grand_SelectSub_rationodesperthres_sem,0),'LineWidth',0.5,'LineStyle','--')
    title({[mouseline,' - ',num2str(nrec),' recordings '], ...
    'Weighted clustering','from all highest subset per laser condition'})
    xlabel('threshold / nedges')
    ylabel('Ratio of SelectSub nodes above threshold')
    % ylim([1e-3 1])
end

subplot(2,4,4)
legend('laser 1', 'laser 2', 'laser 3')
subplot(2,4,8)
legend('laser 1', 'laser 2', 'laser 3')
  
% Same figure but with mean clustering of SoundP cells
figure;
for k=laserconditions   
    subplot(2,4,k), hold on;
    a=[recstruct.clusteringperthres_SelectSub_SoundP_bin];
    loglog(thres_weight,a(:,k:3:end)','LineWidth',2,'Color',cm(k,:))
    title({[mouseline,' - ',num2str(nrec),' recordings - laser ',num2str(k)], ...
    'Binary clustering of highest subset','Sound+ cells for each recording'})
    xlabel('threshold / nedges')
    ylabel('Binary clustering')
    % ylim([1e-5 1])

    subplot(2,4,4), hold on;
    a=[recstruct.clusteringperthres_SelectSub_SoundP_bin];
    loglog(thres_weight,mean(a(:,k:3:end)'),'LineWidth',2)
    % loglog(thres_weight,Grand_SelectSub_ratioedgesperthres_mean+Grand_SelectSub_ratioedgesperthres_sem,'LineWidth',0.5,'LineStyle','--')
    % loglog(thres_weight,max(Grand_SelectSub_ratioedgesperthres_mean-Grand_SelectSub_ratioedgesperthres_sem,0),'LineWidth',0.5,'LineStyle','--')
    title({[mouseline,' - ',num2str(nrec),' recordings '], ...
    'Binary clustering for Sound+ cells','from all highest subset per laser condition'})
    xlabel('threshold / nedges')
    ylabel('Mean binary clustering of SelectSub ')
    % ylim([1e-5 1])
    % legend('laser 1', 'laser 2', 'laser 3')

    subplot(2,4,4+k), hold on;
    a=[recstruct.clusteringperthres_SelectSub_SoundP_wei];
    loglog(thres_weight,a(:,k:3:end)','LineWidth',2,'Color',cm(k,:))
    title({[mouseline,' - ',num2str(nrec),' recordings - laser ',num2str(k)], ...
    'Weighted clustering of highest subset','Sound+ cells for each recording'})
    xlabel('threshold / nedges')
    ylabel('Ratio of SelectSub nodes above threshold')
    % ylim([1e-2 1])

    subplot(2,4,8), hold on;
    a=[recstruct.clusteringperthres_SelectSub_SoundP_wei];
    loglog(thres_weight,mean(a(:,k:3:end)'),'LineWidth',2)
    % loglog(thres_weight,Grand_SelectSub_rationodesperthres_mean+Grand_SelectSub_rationodesperthres_sem,'LineWidth',0.5,'LineStyle','--')
    % loglog(thres_weight,max(Grand_SelectSub_rationodesperthres_mean-Grand_SelectSub_rationodesperthres_sem,0),'LineWidth',0.5,'LineStyle','--')
    title({[mouseline,' - ',num2str(nrec),' recordings '], ...
    'Weighted clustering of Sound+ cells','from all highest subset per laser condition'})
    xlabel('threshold / nedges')
    ylabel('Ratio of SelectSub nodes above threshold')
    % ylim([1e-3 1])
end

subplot(2,4,4)
legend('laser 1', 'laser 2', 'laser 3')
subplot(2,4,8)
legend('laser 1', 'laser 2', 'laser 3')  
  
 %% 4) connectivity matrices averaged across activity groups
 
for j=1:nrec
    Ind_CellSelection=recstruct(j).activityworkspace.Ind_CellSelection;
    % red cells
    recstruct(j).CellID(10).name='red';
    recstruct(j).CellID(10).index= Ind_CellSelection{2}';
    % sound+ cells -> laserp -> soundn -> lasern
    recstruct(j).CellID(1).name='soundplasern';
    recstruct(j).CellID(1).index=setdiff(intersect(Ind_CellSelection{3},Ind_CellSelection{6}),recstruct(j).CellID(10).index)'; % sound+ and Laser-, not red
    recstruct(j).CellID(2).name='soundplaser0';
    recstruct(j).CellID(2).index=setdiff(setdiff(setdiff(Ind_CellSelection{3},Ind_CellSelection{5}),Ind_CellSelection{6}),recstruct(j).CellID(10).index)'; % sound+ and Laser0, not red
    recstruct(j).CellID(3).name='soundplaserp';
    recstruct(j).CellID(3).index=setdiff(intersect(Ind_CellSelection{3},Ind_CellSelection{5}),recstruct(j).CellID(10).index)'; % sound+ and Laser+, not red
    recstruct(j).CellID(4).name='sound0laserp';
    recstruct(j).CellID(4).index=setdiff(setdiff(setdiff(Ind_CellSelection{5},Ind_CellSelection{3}),Ind_CellSelection{4}),recstruct(j).CellID(10).index)'; % sound0 and Laser+, not red
    recstruct(j).CellID(5).name='soundnlaserp';
    recstruct(j).CellID(5).index=setdiff(intersect(Ind_CellSelection{4},Ind_CellSelection{5}),recstruct(j).CellID(10).index)'; % sound- and Laser+, not red
    recstruct(j).CellID(6).name='soundnlaser0';
    recstruct(j).CellID(6).index=setdiff(setdiff(setdiff(Ind_CellSelection{4},Ind_CellSelection{5}),Ind_CellSelection{6}),recstruct(j).CellID(10).index)'; % sound- and Laser0, not red
    recstruct(j).CellID(7).name='soundnlasern';
    recstruct(j).CellID(7).index=setdiff(intersect(Ind_CellSelection{4},Ind_CellSelection{6}),recstruct(j).CellID(10).index)'; % sound- and Laser-, not red
    recstruct(j).CellID(8).name='sound0lasern';
    recstruct(j).CellID(8).index=setdiff(setdiff(setdiff(Ind_CellSelection{6},Ind_CellSelection{3}),Ind_CellSelection{4}),recstruct(j).CellID(10).index)'; % sound0 and Laser-, not red
    recstruct(j).CellID(9).name='sound0laser0';
    recstruct(j).CellID(9).index=setdiff(setdiff(setdiff(setdiff(setdiff(Ind_CellSelection{1}',Ind_CellSelection{3}),Ind_CellSelection{4}),Ind_CellSelection{5}),Ind_CellSelection{6}),recstruct(j).CellID(10).index)'; % sound0 and Laser0, not red

    recstruct(j).Order=[recstruct(j).CellID.index];
end
 
%% plot the matrices re ordered by activity
% one figure per laser condition, all recordings, without averaging
% using norm L2 to plot matrices

for k=laserconditions
    figure;
    for j=1:nrec
        subplot(4,4,j), hold on;
        imagesc(squeeze(recstruct(j).subset{k}.subset_matnormL2(recstruct(j).SelectSubset(k),recstruct(j).Order,recstruct(j).Order)))    
        % y axis
        line([-1, -1],[1, length([recstruct(j).CellID(1:3).index])],'Color',[1 0.9 0],'LineWidth',2) %sound+ cells
        line([-2, -2],[1+length([recstruct(j).CellID(1:2).index]), 1+length([recstruct(j).CellID(1:5).index])],'Color','r','LineWidth',2 ) % laser+ cells
        line([-1, -1],[1+length([recstruct(j).CellID(1:4).index]), 1+length([recstruct(j).CellID(1:7).index])],'Color','g','LineWidth',2) % sound- cells
        line([-2, -2],[1+length([recstruct(j).CellID(1:6).index]), 1+length([recstruct(j).CellID(1:8).index])],'Color','b','LineWidth',2) % laser- cells
        line([-2, -2],[1, length([recstruct(j).CellID(1).index])],'Color','b','LineWidth',2) % laser- and sound+ cells  first part
        line([-3, -3],[length(recstruct(j).Order)-length([recstruct(j).CellID(10).index]), length(recstruct(j).Order)],'Color','r','LineWidth',3) % red cells
        % x axis
        line([1, length([recstruct(j).CellID(1:3).index])],[-1, -1],'Color',[1 0.9 0],'LineWidth',2) %sound+ cells
        line([1+length([recstruct(j).CellID(1:2).index]), 1+length([recstruct(j).CellID(1:5).index])],[-2, -2],'Color','r','LineWidth',2 ) % laser+ cells
        line([1+length([recstruct(j).CellID(1:4).index]), 1+length([recstruct(j).CellID(1:7).index])],[-1, -1],'Color','g','LineWidth',2) % sound- cells
        line([1+length([recstruct(j).CellID(1:6).index]), 1+length([recstruct(j).CellID(1:8).index])],[-2, -2],'Color','b','LineWidth',2) % laser- cells
        line([1, length([recstruct(j).CellID(1).index])],[-2, -2],'Color','b','LineWidth',2) % laser- and sound+ cells  first part
        line([length(recstruct(j).Order)-length([recstruct(j).CellID(10).index]), length(recstruct(j).Order)],[-3, -3],'Color','r','LineWidth',3) % red cells
        
        title({[mouseline,' - laser ', num2str(k)], ['matrix of rec ', num2str(j)]})
        xlabel('nodes')
        ylabel('nodes')
        axis tight square
        set(gca,'YDir','reverse') 
    end
end

%% average subset matrix by cell activity

for j=1:nrec
    for k=recstruct(j).laserconditions
        recstruct(j).subset{k}.subset_grouped=NaN(size(recstruct(j).subset{k}.subset,1),10,10); % 10 is the number of groups
        for l=1:size(recstruct(j).subset{k}.subset,1)
            for i=1:10   
                if ~isempty(recstruct(j).CellID(i).index)
                    for m=1:10
                        if ~isempty(recstruct(j).CellID(m).index)
                            recstruct(j).subset{k}.subset_grouped(l,i,m)=mean(mean(recstruct(j).subset{k}.subset_matnormL2(l,recstruct(j).CellID(i).index,recstruct(j).CellID(m).index),'omitnan'),'omitnan');                            
                        end
                    end
                end
            end
        end
    end
end

for k=laserconditions
    figure;
    for j=1:nrec
        subplot(4,4,j), hold on;
        imagesc(squeeze(recstruct(j).subset{k}.subset_grouped(recstruct(j).SelectSubset(k),:,:)),[0, 0.025])    
 
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
        
        title({[mouseline,' - laser ', num2str(k)], ['grouped matrix of rec ', num2str(j)]})
        xlabel('nodes')
        ylabel('nodes')
        axis tight square
        set(gca,'YDir','reverse') 
        colorbar
    end
end

%% average across recordings

Grand_SelectedSubset_grouped_all=zeros(length(laserconditions),nrec,10,10);

for j=1:nrec
    for k=recstruct(j).laserconditions
        Grand_SelectedSubset_grouped_all(k,j,:,:)=recstruct(j).subset{k}.subset_grouped(recstruct(j).SelectSubset(k),:,:);
    end
end

figure;
for k=laserconditions
    subplot(1,4,k),
    imagesc(squeeze(nanmean(Grand_SelectedSubset_grouped_all(k,:,:,:),2)),[0, 0.015])    

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
        
    title({[mouseline,' - ',num2str(nrec),' recordings '], ...
        ['mean matrix - laser ', num2str(k)]})
    xlabel('nodes')
    ylabel('nodes')
    axis tight square
    set(gca,'YDir','reverse') 
    colorbar
end

%% calculate efficiency? 