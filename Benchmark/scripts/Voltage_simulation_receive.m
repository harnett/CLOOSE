%% Linlin Data

close all
clear all

filename = 'C:\Users\Harnettlab\OneDrive\Desktop\Voltage_Data\Sq_camera.bin';
xpixels = 96;
ypixels = 200;
mov = readBinMov(filename, xpixels, ypixels);
fft_tts  = fft2(mean(mov(:, :, 1:50), 3)); 
clear mov

IPv4 = '127.0.0.1';
port = 30001;          % Port number used for communication
buffsz = 20;    

rois = load('C:\Users\Harnettlab\OneDrive\Desktop\Voltage_Data\ROI.mat');
roiMask = cell(numel(rois.RoiInfo),1);
for i = 1%1:numel(rois.RoiInfo)
    roiMask{i, 1} = poly2mask(rois.RoiInfo(1,i).Position(:, 1),rois.RoiInfo(1,i).Position(:, 2),xpixels,ypixels);
end
nroi = 1; %size(roiMask, 1);

fs = 1000;  % Sampling frequency in Hz
cutoff = 50; % High-pass cutoff frequency in Hz

% Design a high-pass filter using a Butterworth filter
filterOrder = 4; % Choose a reasonable filter order
d = designfilt('highpassiir', 'FilterOrder', filterOrder, ...
               'HalfPowerFrequency', cutoff, 'SampleRate', fs, ...
               'DesignMethod', 'butter');
% Apply the filter using filtfilt for zero-phase filtering

noise_level = [8.4568174464555; 6.42921065016092; 6.3283443718719; 5.45695172669764; 4.94371695815561]; 
% BCI Part
pl = 0;
roi_one_frame = nan(nroi, 30000);

% Establish loop speed
avg_f = 1;
avg_mc  = 50;
if avg_f ~= 1
    avg_loop = avg_f;
elseif avg_f == 1
    avg_loop = avg_mc;
end
raw_Data_stack = uint16(zeros(xpixels, ypixels, avg_loop));
f_rate = nan(nroi, 30000/avg_loop);
filtered_signal = nan(nroi, 30000);

s = nan(1, 30000);
s_img_pro = nan(1, 30000);
s_tranf = nan(1, 30000);
t_filt = nan(1, 30000/avg_loop);
t_peaks = nan(5, 30000/avg_loop);
s_mc = nan(1, 30000/avg_loop);
s_F_ext = nan(1, 30000/avg_loop);
 
mc = 1;

warning('off', 'all')

tcpipServer = tcpip(IPv4, port, 'NetworkRole', 'client', ...
    'InputBufferSize', (xpixels*ypixels*buffsz), 'Terminator', 'CR/LF');

% Open the TCP/IP connection.
fopen(tcpipServer);


while pl < 30000

    pl = pl +1;
    if pl == 1
        timerVal = tic;
    end
    s(pl) = toc(timerVal);
    
    s_tranf_tic = tic;
    vect_img = uint16(fread(tcpipServer, xpixels*ypixels, 'uint16'));
    s_tranf(pl) = toc(s_tranf_tic);

    t_img_pro = tic;
    rawData_tmp = reshape(vect_img, xpixels, ypixels);
    s_img_pro(pl) = toc(t_img_pro);
    
    rm = rem(pl, avg_loop);

    if rm ~= 0
        raw_Data_stack(:, :, rm) = rawData_tmp;
        s_img_pro(pl) = toc(t_img_pro);
        continue
    elseif rm == 0
        raw_Data_stack(:, :, avg_loop) = rawData_tmp;
        raw_Data_onePl = mean(raw_Data_stack, 3);
        s_img_pro(pl) = toc(t_img_pro);

        % Motion Correction
        s_mc_tic  = tic;
        if mc == 1
            output = dftregistration(fft_tts, fft2(raw_Data_onePl));
            xsh = output(1, 4);
            ysh = output(1, 3);
        end
        



        if avg_f ~= 1
            % Just extraxct fluorescence from average image
            raw_Data_onePl_mc = imtranslate(raw_Data_onePl, [xsh ysh]);
            for i = 1:nroi
                roi_one_frame(i,pl) = double(mean(raw_Data_onePl_mc(roiMask{i, 1})));
            end
        elseif avg_f == 1
            % Apply MC to the entire stack and extract
            raw_Data_stack_mc = imtranslate(raw_Data_stack, [xsh ysh]);
            s_mc(1, pl/avg_loop) = toc(s_mc_tic);

            s_F_ext_tic = tic;
            reshapedStack = reshape(raw_Data_stack_mc, [], size(raw_Data_stack_mc, 3));
            for i = 1:nroi
                maskVector = roiMask{i, 1}(:);
                maskedPixels = reshapedStack(maskVector, :);
                roi_one_frame(i,pl - avg_loop + 1: pl) = mean(maskedPixels, 1);
            end
            s_F_ext(1, pl/avg_loop) = toc(s_F_ext_tic);


            idx = pl - avg_loop + 1: pl;
            ticFilt = tic;
            filtered_signal(1:nroi,idx) = filtfilt(d, roi_one_frame(:,idx)')';
            t_filt(:, pl/avg_loop) = toc(ticFilt);
            tic_peaks = tic ;
            peaks = cell(nroi, 1);
            for iflt = 1:size(filtered_signal, 1)
                if any(filtered_signal(iflt,idx) > noise_level(iflt)*4)
                    [peaks{iflt, 1}, ~] = findpeaks(filtered_signal(iflt, idx), 'MinPeakHeight', noise_level(iflt)*4);
                end
                t_peaks(iflt, pl/avg_loop) = toc(tic_peaks);
            end 
%             [peaks, ~] = arrayfun(@(i) findpeaks(filtered_signal(i, idx), 'MinPeakHeight', noise_level(i)*4), ...
%                          1:size(filtered_signal, 1), 'UniformOutput', false);
            f_rate(:, pl/avg_loop) = cellfun(@length, peaks);
            
        end
        raw_Data_stack = uint16(zeros(xpixels, ypixels, avg_loop));
    end
end
toc(timerVal)

all_mu = cat(1, nanmean(s_tranf), nanmean(s_img_pro),nanmean(s_mc), nanmean(s_F_ext), nanmean(t_filt),nanmean(t_peaks(1, :)), nanmean(diff(s)))*1000;

sem_diff = nanstd(diff(s)) / sqrt(sum(~isnan(diff(s)))) * 1000;
sem_s_img_pro = nanstd(s_img_pro) / sqrt(sum(~isnan(s_img_pro))) * 1000;
sem_s_tranf = nanstd(s_tranf) / sqrt(sum(~isnan(s_tranf))) * 1000;
sem_t_filt = nanstd(t_filt) / sqrt(sum(~isnan(t_filt))) * 1000;
sem_t_peaks = nanstd(t_peaks(1, :)) / sqrt(sum(~isnan(t_peaks(1, :)))) * 1000;
sem_s_mc = nanstd(s_mc) / sqrt(sum(~isnan(s_mc))) * 1000;
sem_s_F_ext = nanstd(s_F_ext) / sqrt(sum(~isnan(s_F_ext))) * 1000;

% Create SEM vector
all_sem = cat(1, sem_s_tranf, sem_s_img_pro,  ...
    sem_s_mc, sem_s_F_ext, sem_t_filt, sem_t_peaks, sem_diff);

% Labels for each bar (adjust as needed)
labels = {'Transfer', 'Img pro',  'MC', 'F Ext', 'Filter', 'Pk Detect', 'diff(s)'};

% Create the bar plot
figure;
bar_handle = bar(all_mu, 'FaceColor','w');  % Optional: set a nice color

% Hold the plot for adding error bars
hold on;

% Add error bars
errorbar(1:length(all_mu), all_mu, all_sem, 'k.', 'LineWidth', 1.5);

% Customize the plot
set(gca, 'XTick', 1:length(labels), 'XTickLabel', labels, 'XTickLabelRotation', 45, 'TickDir', 'out');
ylabel('Time (ms)');
title('Mean values with SEM error bars');
grid off;
box off;
