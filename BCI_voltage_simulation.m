%% Linlin Data

close all
clear all

filename = 'Z:\CLOOSE_benchmark_tests\Voltage_Imaging\Sq_camera.bin';
nrow = 96;
ncol = 200;
mov = readBinMov(filename, nrow, ncol);
figure(100); imshow(imadjust(uint16(mean(mov(:, :, 1:100), 3))));
rois = load('Z:\CLOOSE_benchmark_tests\Voltage_Imaging\ROI.mat');
roiMask = cell(numel(rois.RoiInfo),1);
for i = 1:numel(rois.RoiInfo)
    drawpolygon('Position', rois.RoiInfo(1,i).Position,...
        'LineWidth', 1, 'Deletable', true, 'Label', num2str(i), ...
        'LabelVisible', 'hover', 'Color', 'y', 'FaceAlpha', 0);
    roiMask{i, 1} = poly2mask(rois.RoiInfo(1,i).Position(:, 1),rois.RoiInfo(1,i).Position(:, 2),nrow,ncol);
end
nroi = size(roiMask, 1);

s = nan(1, size(mov, 3));
fft_tts  = fft2(mean(mov(:, :, 1:50), 3)); figure(1); imshow(imadjust(uint16(mean(mov(:, :, 1:50), 3))))
mc = 1;

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
roi_one_frame = nan(numel(rois.RoiInfo), size(mov, 3));

% Establish loop speed
avg_f = 1;
avg_mc  = 50;
if avg_f ~= 1
    avg_loop = avg_f;
elseif avg_f == 1
    avg_loop = avg_mc;
end
raw_Data_stack = uint16(zeros(nrow, ncol, avg_loop));
f_rate = nan(numel(rois.RoiInfo), size(mov, 3)/avg_loop);
filtered_signal = nan(numel(rois.RoiInfo), size(mov, 3));
warning('off', 'all')

while pl < size(mov, 3)

    pl = pl +1;
    if pl == 1
        timerVal = tic;
    end
    s(pl) = toc(timerVal);

    rawData_tmp = mov(:, :, pl);

    rm = rem(pl, avg_loop);

    if rm ~= 0
        raw_Data_stack(:, :, rm) = rawData_tmp;
        continue
    elseif rm == 0
        raw_Data_stack(:, :, avg_loop) = rawData_tmp;
        raw_Data_onePl = mean(raw_Data_stack, 3);

        % Motion Correction
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
            reshapedStack = reshape(raw_Data_stack_mc, [], size(raw_Data_stack_mc, 3));

            for i = 1:nroi
                maskVector = roiMask{i, 1}(:);
                maskedPixels = reshapedStack(maskVector, :);
                roi_one_frame(i,pl - avg_loop + 1: pl) = mean(maskedPixels, 1);
            end
            idx = pl - avg_loop + 1: pl;
            filtered_signal(1:nroi,idx) = filtfilt(d, roi_one_frame(:,idx)')';
            peaks = cell(nroi, 1);
            for iflt = 1:size(filtered_signal, 1)
                if any(filtered_signal(iflt,idx) > noise_level(iflt)*4)
                    [peaks{iflt, 1}, ~] = findpeaks(filtered_signal(iflt, idx), 'MinPeakHeight', noise_level(iflt)*4);
                end
            end 
%             [peaks, ~] = arrayfun(@(i) findpeaks(filtered_signal(i, idx), 'MinPeakHeight', noise_level(i)*4), ...
%                          1:size(filtered_signal, 1), 'UniformOutput', false);
            f_rate(:, pl/avg_loop) = cellfun(@length, peaks);
        end
        raw_Data_stack = uint16(zeros(nrow, ncol, avg_loop));
    end
end
toc(timerVal)

figure(200)
iflt = 4;
bla = roi_one_frame(iflt, 10:end);
bla2 = filtered_signal(iflt, 10:end);

norm_data = (bla - min(bla)) / ( max(bla) - min(bla) );
norm_data2 = (bla2 - min(bla2)) / ( max(bla2) - min(bla2) );

plot(norm_data(6150:6450) + 1)
hold on
plot(norm_data2(6150:6450))
hold off


figure(2500)
findpeaks(filtered_signal(iflt, 10:end), 'MinPeakHeight', noise_level(iflt)*4);
hold off


figure(400)
iflt = 1;
bla = roi_one_frame(iflt, 10:end);
bla2 = filtered_signal(iflt, 10:end);

norm_data = (bla - min(bla)) / ( max(bla) - min(bla) );
norm_data2 = (bla2 - min(bla2)) / ( max(bla2) - min(bla2) );

plot(norm_data(200:500) + 1)
hold on
plot(norm_data2(200:500))
hold off


figure(4500)
findpeaks(filtered_signal(iflt, 10:end), 'MinPeakHeight', noise_level(iflt)*4);
hold off

figure(250)
plot(roi_one_frame(1, 200:250))
hold on
plot(roi_one_frame(4, 200:250))
hold off

figure(300)
findpeaks(filtered_signal(1, 10:end))

figure(400)
noise_level = std(abs(diff(filtered_signal(2, 10:end))));
plot(filtered_signal(1, 10:end))
hold on 
yline(noise_level*5)


figure(500)
noise_level = std(abs(diff(filtered_signal(3, 10:end))));
plot(roi_one_frame(3, 10:end))
hold on
findpeaks(filtered_signal(3, 10:end), 'MinPeakHeight', noise_level*4)
hold off
% 
% 
% noise_level = std(abs(diff(filtered_signal(:, 10:end), 1, 2)), [], 2);

for ifig = 6:10
    figure(ifig); imshow(imadjust(mov(:, :, ifig)))
end

pause()