close all
clear all

tt1 = load('Z:\CLOOSE benchmark tests\testing_bci_bench_marktest\halfimage_uint8_external.mat');
tt2 = load('Z:\CLOOSE benchmark tests\testing_bci_bench_marktest\fullimage_uint8_external.mat');
tt3 = load('Z:\CLOOSE benchmark tests\testing_bci_bench_marktest\halfimage_uint16_external.mat');
tt4 = load('Z:\CLOOSE benchmark tests\testing_bci_bench_marktest\fullimage_uint16_external.mat');

tt5 = load('Z:\CLOOSE benchmark tests\testing_bci_bench_marktest\halfimage_uint8_internal.mat');
tt6 = load('Z:\CLOOSE benchmark tests\testing_bci_bench_marktest\fullimage_uint8_internal.mat');
tt7 = load('Z:\CLOOSE benchmark tests\testing_bci_bench_marktest\halfimage_uint16_internal.mat');
tt8 = load('Z:\CLOOSE benchmark tests\testing_bci_bench_marktest\fullimage_uint16_internal.mat');

tts = (cat(1,tt5.tt(2:end),tt6.tt(2:end),tt7.tt(2:end),tt8.tt(2:end), tt1.tt(2:end), tt2.tt(2:end),tt3.tt(2:end),tt4.tt(2:end))).*1000; 

all_mu = mean(tts, 2);
all_med = median(tts, 2);

SEMs = (std(tts, [], 2))./sqrt(size(tts, 2));

figure(1)
bar(all_mu, 'FaceColor', 'w') 
hold on
errorbar(all_mu, SEMs, '.', 'Color', 'k') 
box off
set(gca,'TickDir','out'); % The only other option is 'in'
hold off

x1 = uint16(10);
whos x1

x2 = uint8(10);
whos x2

figure(2)
violin(tts') 
box off
set(gca,'TickDir','out'); % The only other option is 'in'
hold off

xpx_f = 512;
xpx_h = 256;
ypx = 796;
bit8 = 1;
bit16 = 2;
bytes = [xpx_h*ypx*bit8; xpx_h*ypx*bit16; xpx_f*ypx*bit8; xpx_f*ypx*bit16; ...
    xpx_h*ypx*bit8; xpx_h*ypx*bit16; xpx_f*ypx*bit8; xpx_f*ypx*bit16];

figure(3)
scatter(bytes(1:4), all_mu(1:4), 'filled', 'r')
hold on
scatter(bytes(5:8), all_mu(5:8),'filled', 'b')
hold on
coefficients = polyfit(bytes(1:4), all_mu(1:4), 1);
xFit = linspace(min(bytes(1:4)), max(bytes(1:4)), 1000);
yFit = polyval(coefficients , xFit);
plot(bytes(1:4), all_mu(1:4), 'r.', 'MarkerSize', 15); % Plot training data.
hold on; % Set hold on so the next plot does not blow away the one we just drew.
plot(xFit, yFit, 'r-', 'LineWidth', 1); % Plot fitted line.
hold on
coefficients = polyfit(bytes(5:8), all_mu(5:8), 1);
xFit = linspace(min(bytes(5:8)), max(bytes(5:8)), 1000);
yFit = polyval(coefficients , xFit);
plot(bytes(5:8), all_mu(5:8), 'b.', 'MarkerSize', 15); % Plot training data.
hold on; % Set hold on so the next plot does not blow away the one we just drew.
plot(xFit, yFit, 'b-', 'LineWidth', 1); % Plot fitted line.
xlim([0 900000])
xticks(0:450000:900000)
ylim([0 8])
yticks(0:4:8)
box off
set(gca,'TickDir','out'); % The only other option is 'in'
hold off


