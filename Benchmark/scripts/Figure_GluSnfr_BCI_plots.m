%% Part 1 : Single Neuron
close all
clear all

rng(1)

% load('Z:\CLOOSE_benchmark_tests\Spine_imaging\DF_gcampspinebci04_000_006.mat');
% load('Z:\BCI\glusnfr_spine_bci\glusnfrbci4\Baseline\baselineMask_1.mat')
% load('Z:\CLOOSE_benchmark_tests\Spine_imaging\ROIS_gcampspinebci04_000_006.mat')

load('Z:\CLOOSE_benchmark_tests\Spine_imaging\DF_gcampspinebci04_000_006.mat');
load('Y:\users\zilan662\Glusnfr_spine_metadata\glusnfr_spine_bci\061324_glusnfrspinebci5\Baseline\baselineMask_1.mat')
load('Z:\CLOOSE_benchmark_tests\Spine_imaging\ROIS_gcampspinebci04_000_006.mat')

% figure(1)
% imshow(imsharpen(imadjust(baseline.RefImage, [0 0.8],[]),'Amount', 0.8 ))
% hold on
% for i = 1:numel(ROI)
%     ROIq(i) = drawpolygon(gca, 'Position',ROI(1,i).Position,...
%         'LineWidth', 1, 'Deletable', true, 'Label', num2str(i), ...
%         'LabelVisible', 'hover', 'Color', 'y', 'FaceAlpha', 0, 'MarkerSize', 2);
% end

% [b_par, stat] = robustfit(DF(3, :), DF(1, :));
% [b_par2, stat2] = robustfit(DF(3, :), DF(2, :));
% 
% subDF1 = stat.resid;
% subDF2 = stat2.resid;
% 
% scatter(subDF2, stat2.resid)
% 
% figure(3)
% plot(stat.resid)

figure(100)
plot(DF(1, :)./max(DF(1, :)) + 2, 'b')
hold on
plot(DF(2, :)./max(DF(2, :)) + 3, 'r')
% hold on
% plot(stat.resid./max(stat.resid), 'b')
% hold on
% plot(stat2.resid./max(stat2.resid)+ 1, 'r')

hold off


oneTrial_length = 435;
pop1 = 1;
pop2 = 2;
targetAngle = 90;
trgP1 = 2;
trgP2 = 1;

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
    fff{itrial} = cat(1, df_oneTrial{itrial}(pop1, :), df_oneTrial{itrial}(pop2, :));
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

x = mean(DF(pop1, :), 1) -  mean(DF(pop2, :), 1);
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
plot(DF(pop1, 4000:5300)./max(DF(pop1, :)) + 4, 'r')
hold on
plot(DF(pop2, 4000:5300)./max(DF(pop2, :)) + 3, 'b')
hold on
plot((DF(pop1, 4000:5300) - DF(pop2, 4000:5300))./(max(DF(pop1, :) - DF(pop2, :))) + 2, 'k')
hold on
plot(xx(4000:5300)./max(xx))
hold on
yline(baseline.selected_tresh./max(DF(pop1, :) - DF(pop2, :)) + 2, '--k')
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
            DF_avg_succ(ii, :) = mean(fff{ii}(pop1, id_s - 15 : id_s + 59), 1);
            DF_avg_notsucc(ii, :) = mean(fff{ii}(pop2, id_s - 15 : id_s + 59), 1);
        catch
        end
    end
end
% 
% DF_avg_notsucc = DF_avg_notsucc - mean(nanmean(DF_avg_succ(:, 1:5)));
SEM2 = nanstd(DF_avg_notsucc)./sqrt(sum(~isnan(DF_avg_notsucc(:, 1))));
% DF_avg_succ = DF_avg_succ - mean(nanmean(DF_avg_succ(:, 1:5)));

SEM1 = nanstd(DF_avg_succ)./sqrt(sum(~isnan(DF_avg_succ(:, 1))));

figure(100)
plot(nanmean(DF_avg_succ(:, 1:60)), 'r')
hold on
shadedplot(1:60, nanmean(DF_avg_succ(:, 1:60))+ SEM1(:, 1:60), nanmean(DF_avg_succ(:, 1:60)) - SEM1(:, 1:60), 'r')
hold on
% ylim([-0.2 1.5])
% yticks(0:0.5:1.5)
% xlim([0 75])
% xticks(0:15:75)
plot(nanmean(DF_avg_notsucc(:, 1:60)), 'b')
hold on
shadedplot(1:60, nanmean(DF_avg_notsucc(:, 1:60))+ SEM2(:, 1:60), nanmean(DF_avg_notsucc(:, 1:60)) - SEM2(:, 1:60), 'b')
hold on
box off
set(gca,'TickDir','out'); % The only other option is 'in'
hold off
