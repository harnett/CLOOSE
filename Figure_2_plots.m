%% Part 1 : Single Neuron
close all
clear all

rng(1)

load('Z:\CLOOSE benchmark tests\Various_BCI_approaches\DF_062624_NDNF_07_000_009.mat');
load('Z:\BCI\Data\NDNF_Opto\NDNF_07\062624_NDNF_07\Baseline\baselineMask_2.mat')
load('Z:\CLOOSE benchmark tests\Various_BCI_approaches\ROIS_062624_NDNF_07_000_009.mat')

figure(1)
imshow(imsharpen(imadjust(baseline.RefImage, [0 0.3],[]),'Amount', 0.8 ))
hold on
for i = 1:numel(ROI)
    ROIq(i) = drawpolygon(gca, 'Position',ROI(1,i).Position,...
        'LineWidth', 1, 'Deletable', true, 'Label', num2str(i), ...
        'LabelVisible', 'hover', 'Color', 'y', 'FaceAlpha', 0, 'MarkerSize', 2);
end

oneTrial_length = 435;
pop1 = 1;
pop2 = 2;
targetAngle = 90;
trgP1 = 3;
trgP2 = 1;

DF = DF([trgP1, trgP2], :);

trial_start = randperm(size(DF,2), 200);
df_oneTrial  = cell(1, length(trial_start));
z_zz = cell(1, length(trial_start));
F0_all = cell(1, length(trial_start));
F_2 = repmat(DF(:, :), 1, 2); 
F0_tmp = nan(size(DF(:, :)));

for i = 1:size(DF(:, :), 2)
    F0_tmp(:,i) = prctile(F_2(:,i:i+(oneTrial_length*4)-1), 10, 2);
end
F0_tmp= (repmat(F0_tmp, 1, 2));

for itrial = 1:length(trial_start)
    itrial
    beg_end = [trial_start(itrial)+size(DF(:, :),2)-oneTrial_length; ...
        trial_start(itrial)+size(DF(:, :),2)];
    
    F0_all{itrial} = F0_tmp(:,beg_end(1):beg_end(2));
    df_oneTrial{itrial} =  (F_2(:, trial_start(itrial) : trial_start(itrial) + oneTrial_length));
    if ~isempty(pop2)
        z_zz{itrial}= nanmean(df_oneTrial{itrial}(pop1, :),1) - nanmean(df_oneTrial{itrial}(pop2,:),1);
    else
        z_zz{itrial}= nanmean(df_oneTrial{itrial},1);
    end
end

if ~isempty(pop2)
    dfmx1 = max(cat(2, z_zz{:}));
    dfmx2 = min(cat(2, z_zz{:}));
    threshs = dfmx2: 0.005: dfmx1;  
else
    threshs = 0:0.01:max(DF, [], 'all');
end

baseline.zzz = z_zz;
Day = 1; 

propp = nan(1,length(threshs));
for i = 1:length(threshs)
    total_large = nan(1,length(z_zz));
    for ii =  1:length(z_zz)
        mvMu = movmean(z_zz{ii}, 3);
        %         mvMu = z_zz{ii};
        if sum(mvMu > threshs(i)) > 0
            total_large(ii) = 1;
        else
            total_large(ii) = 0;
        end
    end
    propp(1,i) = sum(total_large)/length(total_large);
    propp(2,i) = threshs(i);
end
baseline.propTrialsTresh = propp;
Thresh.propTrialsTresh = baseline.propTrialsTresh;
baseline.selected_tresh = propp(2,min(find(propp(1,:) <= 0.5)));
[muFit, sigmaFit] = normfit(cat(2, z_zz{:}));
baseline.muFit = muFit;
baseline.sigmaFit = sigmaFit;
Thresh.muFit = muFit;
Thresh.sigmaFit = sigmaFit;
z_sc = (baseline.selected_tresh - muFit)/sigmaFit;
% baseline.selected_range = muFit - (z_sc*sigmaFit);
baseline.selected_range = min(cat(2, z_zz{:}));

cc =  cat(2, baseline.zzz{:});
% Data on the left side of the distribution
d = cc(cc <= nanmedian(cc));
% Mirror the distribution on the right
dd = max(d)-d;
%concatenate the left and right side of the distribution
ddc = cat(2,d, dd);
%fit
[baseline.muHalfFit, baseline.sigmaHalfFit] = normfit(ddc);
% z score
baseline.z_HalfFit = (baseline.selected_tresh - baseline.muHalfFit)/baseline.sigmaHalfFit;


angle_bins  = mod((targetAngle - 90), 360) : 15 : targetAngle;
act_bins  = nan(1, 7);
act_bins(2:7)  = linspace(baseline.selected_range, baseline.selected_tresh, 6); 
act_bins(1) = act_bins(2) - (act_bins(3) - act_bins(2));
act_bins = act_bins - baseline.selected_range;
rng = baseline.selected_range;

x = DF(pop1, :) -  DF(pop2, :);
xx = nan(size(DF, 2), 1);
for n = 1:size(DF, 2)
    if n >= 3
        z = mean(x(n-2:n)) - rng;
    else
        z = x(n) - rng;
    end


    try
        idx = max(find(z >= act_bins)); %#ok<MXFND>
    catch
        idx = 1;
    end

    xx(n) = angle_bins(idx);
end

figure(1)
plot(DF(1, :)./max(DF(1, :)) + 1, 'r')
hold on
plot(xx./max(xx))
hold on
yline(baseline.selected_tresh./max(DF(1, :)) + 1, '--k')
box off
set(gca,'TickDir','out'); % The only other option is 'in'
hold off


DF_avg_succ = nan(200, 75);
DF_avg_notsucc = nan(200, 75);


total_large = nan(1,length(z_zz));
for ii =  1:length(z_zz)
    mvMu = movmean(z_zz{ii}, 3);
    %         mvMu = z_zz{ii};
    if sum(mvMu > baseline.selected_tresh) > 0
        id_s = min(find(mvMu > baseline.selected_tresh));
        try
            DF_avg_succ(ii, :) = mvMu(id_s - 15 : id_s + 59);
        catch
        end

    else
        try
            DF_avg_notsucc(ii, :) = mvMu(end - 74: end);
        catch
        end
    end
end
DF_avg_succ = DF_avg_succ - mean(nanmean(DF_avg_succ(:, 1:5)));
SEM1 = nanstd(DF_avg_succ)./sqrt(sum(~isnan(DF_avg_succ(:, 1))));
figure(100)
plot(nanmean(DF_avg_succ), 'r')
hold on
shadedplot(1:75, nanmean(DF_avg_succ)+ SEM1, nanmean(DF_avg_succ) - SEM1, 'r')
hold on
ylim([-0.2 1.5])
yticks(0:0.5:1.5)
xlim([0 75])
xticks(0:15:75)
box off
set(gca,'TickDir','out'); % The only other option is 'in'
hold off


%% Part 2 : P+ vs P-
close all
clear all
rng(1)
load('Z:\CLOOSE benchmark tests\Various_BCI_approaches\DF_062624_NDNF_07_000_009.mat');

oneTrial_length = 435;
pop1 = 1:4;
pop2 = 5:8;
targetAngle = 90;
trgP1 = [3, 6, 8, 9];
trgP2 = [2, 10, 11, 12];

DF = DF([trgP1, trgP2], :);

trial_start = randperm(size(DF,2), 200);
df_oneTrial  = cell(1, length(trial_start));
z_zz = cell(1, length(trial_start));
fff = cell(1, length(trial_start));
F0_all = cell(1, length(trial_start));
F_2 = repmat(DF(:, :), 1, 2); 
F0_tmp = nan(size(DF(:, :)));

for i = 1:size(DF(:, :), 2)
    F0_tmp(:,i) = prctile(F_2(:,i:i+(oneTrial_length*4)-1), 10, 2);
end
F0_tmp= (repmat(F0_tmp, 1, 2));

for itrial = 1:length(trial_start)
    itrial
    beg_end = [trial_start(itrial)+size(DF(:, :),2)-oneTrial_length; ...
        trial_start(itrial)+size(DF(:, :),2)];
    
    F0_all{itrial} = F0_tmp(:,beg_end(1):beg_end(2));
    df_oneTrial{itrial} =  (F_2(:, trial_start(itrial) : trial_start(itrial) + oneTrial_length));
    if ~isempty(pop2)
        z_zz{itrial}= nanmean(df_oneTrial{itrial}(pop1, :),1) - nanmean(df_oneTrial{itrial}(pop2,:),1);
    else
        z_zz{itrial}= nanmean(df_oneTrial{itrial},1);
    end
    fff{itrial} = cat(1, df_oneTrial{itrial}(pop1, :), df_oneTrial{itrial}(pop2, :))
end

if ~isempty(pop2)
    dfmx1 = max(cat(2, z_zz{:}));
    dfmx2 = min(cat(2, z_zz{:}));
    threshs = dfmx2: 0.005: dfmx1;  
else
    threshs = 0:0.01:max(DF, [], 'all');
end

baseline.zzz = z_zz;
Day = 1; 

propp = nan(1,length(threshs));
for i = 1:length(threshs)
    total_large = nan(1,length(z_zz));
    for ii =  1:length(z_zz)
        mvMu = movmean(z_zz{ii}, 3);
        %         mvMu = z_zz{ii};
        if sum(mvMu > threshs(i)) > 0
            total_large(ii) = 1;
        else
            total_large(ii) = 0;
        end
    end
    propp(1,i) = sum(total_large)/length(total_large);
    propp(2,i) = threshs(i);
end
baseline.propTrialsTresh = propp;
Thresh.propTrialsTresh = baseline.propTrialsTresh;
baseline.selected_tresh = propp(2,min(find(propp(1,:) <= 0.5)));
[muFit, sigmaFit] = normfit(cat(2, z_zz{:}));
baseline.muFit = muFit;
baseline.sigmaFit = sigmaFit;
Thresh.muFit = muFit;
Thresh.sigmaFit = sigmaFit;
z_sc = (baseline.selected_tresh - muFit)/sigmaFit;
% baseline.selected_range = muFit - (z_sc*sigmaFit);
baseline.selected_range = min(cat(2, z_zz{:}));

cc =  cat(2, baseline.zzz{:});
% Data on the left side of the distribution
d = cc(cc <= nanmedian(cc));
% Mirror the distribution on the right
dd = max(d)-d;
%concatenate the left and right side of the distribution
ddc = cat(2,d, dd);
%fit
[baseline.muHalfFit, baseline.sigmaHalfFit] = normfit(ddc);
% z score
baseline.z_HalfFit = (baseline.selected_tresh - baseline.muHalfFit)/baseline.sigmaHalfFit;


angle_bins  = mod((targetAngle - 90), 360) : 15 : targetAngle;
act_bins  = nan(1, 7);
act_bins(2:7)  = linspace(baseline.selected_range, baseline.selected_tresh, 6); 
act_bins(1) = act_bins(2) - (act_bins(3) - act_bins(2));
act_bins = act_bins - baseline.selected_range;
rng = baseline.selected_range;

x = mean(DF(pop1, :)) -  mean(DF(pop2, :));
xx = nan(size(DF, 2), 1);
for n = 1:size(DF, 2)
    if n >= 3
        z = mean(x(n-2:n)) - rng;
    else
        z = x(n) - rng;
    end


    try
        idx(n) = max(find(z >= act_bins)); %#ok<MXFND>
    catch
        idx(n) = 1;
    end

    xx(n) = angle_bins(idx(n));
end


figure(1)
plot(mean(DF(pop1, :))./max(mean(DF(pop1, :))) + 5, 'r')
hold on
plot(mean(DF(pop2, :))./max(mean(DF(pop2, :))) + 4, 'b')
hold on
plot(mean(DF(pop1, :)) - mean(DF(pop2, :))./max(mean(DF(pop1, :)) - mean(DF(pop2, :))) + 2, 'k')
hold on
plot(xx./max(xx))
hold on
yline(baseline.selected_tresh./max(mean(DF(pop1, :)) - mean(DF(pop2, :))) + 2, '--k')
box off
set(gca,'TickDir','out'); % The only other option is 'in'
hold off


DF_avg_succ = nan(200, 75);
DF_avg_notsucc = nan(200, 75);


total_large = nan(1,length(z_zz));
for ii =  1:length(z_zz)
    mvMu = movmean(z_zz{ii}, 3);
    %         mvMu = z_zz{ii};
    if sum(mvMu > baseline.selected_tresh) > 0
        id_s = min(find(mvMu > baseline.selected_tresh));
        try
            DF_avg_succ(ii, :) = mean(fff{ii}(pop1, id_s - 15 : id_s + 59));
            DF_avg_notsucc(ii, :) = mean(fff{ii}(pop2, id_s - 15 : id_s + 59));
        catch
        end
    end
end
% 
% DF_avg_notsucc = DF_avg_notsucc - mean(nanmean(DF_avg_succ(:, 1:5)));
% SEM2 = nanstd(DF_avg_notsucc)./sqrt(sum(~isnan(DF_avg_notsucc(:, 1))));
% DF_avg_succ = DF_avg_succ - mean(nanmean(DF_avg_succ(:, 1:5)));
% SEM1 = nanstd(DF_avg_succ)./sqrt(sum(~isnan(DF_avg_succ(:, 1))));


figure(100)
plot(nanmean(DF_avg_succ), 'r')
hold on
shadedplot(1:75, nanmean(DF_avg_succ)+ SEM1, nanmean(DF_avg_succ) - SEM1, 'r')
hold on
% ylim([-0.2 1.5])
% yticks(0:0.5:1.5)
% xlim([0 75])
% xticks(0:15:75)
plot(nanmean(DF_avg_notsucc), 'b')
hold on
shadedplot(1:75, nanmean(DF_avg_notsucc)+ SEM2, nanmean(DF_avg_notsucc) - SEM2, 'b')
hold on
box off
set(gca,'TickDir','out'); % The only other option is 'in'
hold off



%% Part 3: Population dynamics

close all
clear all

load('Z:\CLOOSE benchmark tests\Various_BCI_approaches\DF_062624_NDNF_07_000_009.mat');
load('Z:\CLOOSE benchmark tests\Various_BCI_approaches\BlineforPopDynamics.mat')

oneTrial_length = 435;
pop1 = 1;
pop2 = 2;
targetAngle = 90;

DF = DF(2:end, :);

trial_start = randperm(size(DF,2), 200);
df_oneTrial  = cell(1, length(trial_start));
z_zz = cell(1, length(trial_start));
F0_all = cell(1, length(trial_start));
F_2 = repmat(DF(:, :), 1, 2); 
F0_tmp = nan(size(DF(:, :)));

DF_T_bline = DF';
[coeff, score, ~, ~, expl] = pca(DF_T_bline);
mydata_mean = mean(DF_T_bline); %Find mean of data (columns)
mydata_mean = repmat(mydata_mean,size(DF_T_bline, 1),1); %Replicate mean vector to matrix for subtraction
my_data_norm = DF_T_bline - mydata_mean; % Normalize data to zero mean y subtraction
scores_b = my_data_norm*coeff; %Manually calculate scores using PCA coeff and normalized data


for itrial = 1:length(trial_start)
    itrial
    beg_end = [trial_start(itrial)+size(DF(:, :),2)-oneTrial_length; ...
        trial_start(itrial)+size(DF(:, :),2)];

    m_d_m = mydata_mean(1:beg_end(2) - beg_end(1) + 1, :);

%     mydata_mean = mean(DF_T_bline); %Find mean of data (columns)
%     mydata_mean = repmat(mydata_mean,size(F_2(:, trial_start(itrial) : trial_start(itrial) + oneTrial_length), 1),1);
    df_oneTrial{itrial} = (F_2(:, trial_start(itrial) : trial_start(itrial) + oneTrial_length))' - m_d_m ;
    z_zz{itrial} = df_oneTrial{itrial}*coeff; %Manually calculate scores using PCA coeff and normalized data
end

dist_PCs = squareform(pdist(score(:, 2)));

dfmx1 = max(dist_PCs, [], 'all');
[x_m,y_m] = find(dist_PCs==dfmx1);

ref_coord_PC = [score(x_m(1), 1), score(x_m(1), 2)];
dfmx2 = 0;

threshs = linspace(dfmx2, dfmx1, 200);

baseline.zzz = z_zz;
Day = 1; 

propp = nan(1,length(threshs));
all_dist_available  = nan(1, length(z_zz));
for i = 1:length(threshs)
    i
    total_large = nan(1,length(z_zz));
    for ii =  1:length(z_zz)

        ct_m = cat(1, ref_coord_PC, z_zz{ii}(:, 1:2));

        mvMu = squareform(abs(pdist(ct_m)));

%         figure(1)
%         plot(score(:, 1), score(:, 2), 'k')
%         hold on
%         plot(z_zz{ii}(:, 1), z_zz{ii}(:, 2), 'g')
%         hold on
%         scatter(ref_coord_PC(1), ref_coord_PC(2), 'r')
%         hold on
%         title(num2str(min(ct_m, [], 'all')));
%         hold off
% 
%         figure(2)
%         histogram(mvMu(1, 2:end))
%         hold off
% 
%         pause()

        %         mvMu = z_zz{ii};
        if i == length(threshs)
            all_dist_available(ii) = max(mvMu(1, 2:end));
        end
        if sum(mvMu(1, 2:end) < threshs(i)) > 0
            total_large(ii) = 1;
        else
            total_large(ii) = 0;
        end
    end
    propp(1,i) = sum(total_large)/length(total_large);
    propp(2,i) = threshs(i);
end
baseline.propTrialsTresh = propp;
Thresh.propTrialsTresh = baseline.propTrialsTresh;
baseline.selected_tresh = propp(2,min(find(propp(1,:) >= 0.3)));
[muFit, sigmaFit] = normfit(cat(2, z_zz{:}));
baseline.muFit = muFit;
baseline.sigmaFit = sigmaFit;
Thresh.muFit = muFit;
Thresh.sigmaFit = sigmaFit;
z_sc = (baseline.selected_tresh - muFit)/sigmaFit;
% baseline.selected_range = muFit - (z_sc*sigmaFit);
baseline.selected_range = max(all_dist_available);

cc =  cat(2, baseline.zzz{:});
% Data on the left side of the distribution
d = cc(cc <= nanmedian(cc));
% Mirror the distribution on the right
dd = max(d)-d;
%concatenate the left and right side of the distribution
ddc = cat(2,d, dd);
%fit
[baseline.muHalfFit, baseline.sigmaHalfFit] = normfit(ddc);
% z score
baseline.z_HalfFit = (baseline.selected_tresh - baseline.muHalfFit)/baseline.sigmaHalfFit;


angle_bins  = mod((targetAngle - 90), 360) : 15 : targetAngle;
% angle_bins = flip(angle_bins);
act_bins  = nan(1, 7);
act_bins  = linspace(dfmx1, baseline.selected_tresh, 7); 

x = nan(size(DF, 2), size(DF, 1));
xx = nan(size(DF, 2), 1);

figure(2)
plot(score(:, 1), score(:, 2), 'Color', [0.7 0.7 0.7])
xlim([min(score(:, 1)) - 0.2, max(score(:, 1))+0.2])
ylim([min(score(:, 2))- 0.2, max(score(:, 2))+ 0.2])
hold on
scatter(ref_coord_PC(1), ref_coord_PC(2),  72, 'filled', 'r')
hold on
viscircles(repmat(ref_coord_PC, 7, 1), act_bins, 'Color','k', 'LineWidth', 1)
box off
set(gca,'TickDir','out'); % The only other option is 'in'
xticks(-1.5:2.5:3.5)
yticks(-1.5:2:2.5)

for ii = [1, 5, 13]% 1:length(z_zz)
    ii
    
    
    figure(2)
    if ii == 1
        plot(baseline.zzz{ii}(:, 1), baseline.zzz{ii}(:, 2), 'Color', [0.1 0 0])
    elseif ii == 5
        plot(baseline.zzz{ii}(:, 1), baseline.zzz{ii}(:, 2), 'Color',[0.5 0 0])
    elseif ii == 13
        plot(baseline.zzz{ii}(:, 1), baseline.zzz{ii}(:, 2), 'Color',[1 0 0])
    end
    hold on
%     title(num2str(min(ct_m, [], 'all')));
    hold on
%     pause()
end

c = 0;
for i = [1,2,18,29,30]
    figure(1)
    max(DF(i, :))
    c=c+1;
    plot(DF(i, :)./max(DF(i, :)) + c, 'k')
    hold on

end
set(gca,'TickDir','out'); % The only other option is 'in'
hold off
box off


for n = 2800:size(DF, 2)

    df_tmp = DF(:, n)';

    x(n, :) = (DF(:, n)'- mydata_mean(1, :))*coeff;

    dst = abs(pdist(cat(1,ref_coord_PC, x(n, 1:2))));
    z = dst;
    try
        idx(n) = max(find(z <= act_bins)); %#ok<MXFND>
    catch
        idx(n) = 1;
    end

    xx(n) = angle_bins(idx(n));
    
%     figure(1)
%     plot(xx)
%     hold on
% 
%     figure(2)
%     hold on
%     scatter(x(n, 1), x(n, 2), 'r')
%     hold on
% 
% %     figure(3)
% %     hold on
% %     plot(n, dst, 'r')
% %     hold on
% 
%     pause(0.001)


end

%save baseline for reproducibility
% save('Z:\CLOOSE benchmark tests\Various_BCI_approaches\BlineforPopDynamics.mat', 'baseline')
