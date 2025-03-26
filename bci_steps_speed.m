close all
clear all

roi_10 = load('Z:\CLOOSE benchmark tests\testing_bci_bench_marktest\Baseline_10\baselineMask_1.mat');
roi_100 = load('Z:\CLOOSE benchmark tests\testing_bci_bench_marktest\Baseline_100\baselineMask_1.mat');


imgpro_tt = roi_10.baseline.timeStamp_imgprocess(2:end).*1000;  
mc_tt = roi_10.baseline.timeStamp_mc_corr(2:end).*1000;   
Fextr_tt = roi_10.baseline.timeStamp_F_extr(2:end).*1000;    
plt_tt = roi_10.baseline.timeStamp_plot(2:end).*1000;    
stimpres_tt = roi_10.baseline.timeStamp_stim_pres(2:end).*1000;  

imgpro_tt2 = roi_100.baseline.timeStamp_imgprocess(2:end).*1000;  
mc_tt2 = roi_100.baseline.timeStamp_mc_corr(2:end).*1000;   
Fextr_tt2 = roi_100.baseline.timeStamp_F_extr(2:end).*1000;    
plt_tt2 = roi_100.baseline.timeStamp_plot(2:end).*1000;    
stimpres_tt2 = roi_100.baseline.timeStamp_stim_pres(2:end).*1000;  


% figure(1)
% histogram(log10(imgpro_tt))
% 
% figure(2)
% histogram(log10(mc_tt))
% 
% figure(3)
% histogram(log10(Fextr_tt))
% 
% figure(4)
% histogram(log10(plt_tt))
% 
% figure(5)
% histogram(log10(stimpres_tt))


figure(11)
histogram((imgpro_tt), 50, 'FaceColor', 'k')
box off
set(gca,'TickDir','out'); % The only other option is 'in'
xlim([0 3.5])
ylim([0 7000])
yticks(0:3500:7000)
hold off

figure(110)
histogram((imgpro_tt), 50, 'FaceColor', 'k')
box off
set(gca,'TickDir','out'); % The only other option is 'in'
xlim([1 3.5])
ylim([0 70])
yticks(0:35:70)
hold off



figure(22)
histogram((mc_tt), 50, 'FaceColor', 'k')
box off
set(gca,'TickDir','out'); % The only other option is 'in'
xlim([0 25])
ylim([0 8000])
yticks(0:4000:8000)
hold off

figure(220)
histogram((mc_tt), 50, 'FaceColor', 'k')
box off
set(gca,'TickDir','out'); % The only other option is 'in'
xlim([5 25])
ylim([0 202])
yticks(0:100:201)
hold off


figure(33)
histogram((Fextr_tt)./10,25, 'FaceColor', 'k')
box off
set(gca,'TickDir','out'); % The only other option is 'in'
xlim([0 0.75])
xticks(0:0.25:0.75)
ylim([0 9000])
yticks(0:4500:9000)
hold off

figure(330)
histogram((Fextr_tt)./10,25, 'FaceColor', 'k')
box off
set(gca,'TickDir','out'); % The only other option is 'in'
xlim([0.25 0.75])
xticks(0.25:0.25:0.75)
ylim([0 120])
yticks(0:60:120)
hold off


figure(44)
histogram((plt_tt), 50, 'FaceColor', 'k')
box off
set(gca,'TickDir','out'); % The only other option is 'in'
xlim([0 12])
xticks(0:6:12)
ylim([0 6000])
yticks(0:3000:6000)
hold off

figure(440)
histogram((plt_tt), 50, 'FaceColor', 'k')
box off
set(gca,'TickDir','out'); % The only other option is 'in'
xlim([2.5 12.5])
xticks(2.5:5:12.5)
ylim([0 1000])
yticks(0:500:1000)
hold off


figure(55)
histogram((stimpres_tt), 'FaceColor', 'k')
box off
set(gca,'TickDir','out'); % The only other option is 'in'
xlim([0 25])
xticks(0:5:25)
ylim([0 10000])
yticks(0:5000:10000)
hold off


figure(550)
histogram((stimpres_tt), 'FaceColor', 'k')
box off
set(gca,'TickDir','out'); % The only other option is 'in'
xlim([4.5 25])
xticks(5:5:25)
ylim([0 30])
yticks(0:15:30)
hold off



mu_imgpro_tt = mean(roi_10.baseline.timeStamp_imgprocess(2:end).*1000);  
mu_mc_tt = mean(roi_10.baseline.timeStamp_mc_corr(2:end).*1000);   
mu_Fextr_tt = mean(roi_10.baseline.timeStamp_F_extr(2:end).*100);    
mu_plt_tt = mean(roi_10.baseline.timeStamp_plot(2:end).*1000);    
mu_stimpres_tt = mean(roi_10.baseline.timeStamp_stim_pres(2:end).*1000);  

all_mu = cat(1,imgpro_tt,  mc_tt, Fextr_tt./10, plt_tt, stimpres_tt);

all_mu = mean(all_mu, 2);

SEMs = (std(all_mu, [], 2))./sqrt(size(all_mu, 2));

figure(1)
bar(all_mu, 'FaceColor', 'k') 
hold on
errorbar(all_mu, SEMs, '.', 'Color', 'k') 
box off
set(gca,'TickDir','out'); % The only other option is 'in'
ylim([0 5])
yticks([0:2.5:5])
hold off



% 
% 
% figure(10)
% histogram(log10(imgpro_tt2))
% 
% figure(20)
% histogram(log10(mc_tt2))
% 
% figure(30)
% histogram(log10(Fextr_tt2))
% 
% figure(40)
% histogram(log10(plt_tt2))
% 
% figure(50)
% histogram(log10(stimpres_tt2))
% 
% 
% figure(110)
% histogram((imgpro_tt2))
% 
% figure(220)
% histogram((mc_tt2))
% 
% figure(330)
% histogram((Fextr_tt2))
% 
% figure(440)
% histogram((plt_tt2))
% 
% figure(550)
% histogram((stimpres_tt2))
% 
% 
