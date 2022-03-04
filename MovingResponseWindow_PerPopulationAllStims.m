% Code for calculating response on optimal response window per
% population, over all stims (one time window per recording)


%%%% Option 4: like option 2 with moving window, but on whole population:
%%%% take time interval that maximizes the number of signif cells. 

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
    av_raster_poststim_time(ii,jj,kk)=nanmean(b2);
    sem_raster_poststim_time(ii,jj,kk)=nanstd(b2)/sqrt(repeats);
    end
end
end

% for each stim, for each window, calculate how many cells have signif response
% and how many of them are pos or neg responses
% for ii=1:length(stimInfo.index) 
%     for kk=1:length(Tdelay_test)
%         a=find(tstats_time(ii,:,kk)==1);
%         Nsignifcells_time(ii,kk)=sum(tstats_time(ii,:,kk));
%         Nsignifcells_pos_time(ii,kk)=length(find(av_raster_poststim_time(ii,a,kk)>0));
%         Nsignifcells_neg_time(ii,kk)=length(find(av_raster_poststim_time(ii,a,kk)<0));    
%     end
% end

%% Option 4b: Over all stim conditions find first response window that maximizes the number of
% different cells with signif responses
%%%% Note: here we're arbitrarily taking the first window with max number
%%%% of cells, but it could be changed to last of middle point.


    for kk=1:length(Tdelay_test)
        Nsignifcells_allstims_time(kk)=length(find(sum(tstats_time(:,:,kk))>=1));
    end
%     a=find(Nsignifcells_allstims_time(:,:)==max(Nsignifcells_allstims_time));
     Tdelay_allstims=Tdelay_test(find(Nsignifcells_allstims_time(:)==max(Nsignifcells_allstims_time(:)),1,'first'));

%     [a b]=find(Nsignifcells_time(:,:)==max(max(Nsignifcells_time)));
%     Tdelay_allstims=Tdelay_test(b(1));


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
title({['mouse, - ',num2str(date),' - ',num2str(size(raster,3)),' cells'],...
    ['Number of signif cells vs onset of time window'],...
     ['to at least one stim condition stim conditions (blue)'],...
    ['Window begin = ',num2str(Tdelay_allstims),'s, Number of signif cells = ',num2str(max(Nsignifcells_allstims_time))],...
    ['Duration window = ',num2str(Tdelay_windur),'s, alpha = ',num2str(alpha)]})

%     ['For each stim (black/pink/red) and over all stim conditions (blue)'],...

for ii=1:length(stimInfo.index)
    % compute stats, average and std with chosen Tdelay
    tstats(ii,:)=tstats_time(ii,:,find(Nsignifcells_allstims_time(:)==max(Nsignifcells_allstims_time(:)),1,'first'));
    pstats(ii,:)=pstats_time(ii,:,find(Nsignifcells_allstims_time(:)==max(Nsignifcells_allstims_time(:)),1,'first'));
    av_raster_poststim(ii,:)=av_raster_poststim_time(ii,:,find(Nsignifcells_allstims_time(:)==max(Nsignifcells_allstims_time(:)),1,'first'));
    sem_raster_poststim(ii,:)=sem_raster_poststim_time(ii,:,find(Nsignifcells_allstims_time(:)==max(Nsignifcells_allstims_time(:)),1,'first'));
end
