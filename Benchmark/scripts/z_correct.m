%% z-correct
close all
clear all

frames_tmp = sbxgrabframe('Z:\Scanbox\Data\NDNF\NDNF_06\070124_NDNF_06\070124_NDNF_06_000_000\070124_NDNF_06_000_000', 1, 1000);
frames_tmp = squeeze(frames_tmp); 
load('Z:\BCI\Data\NDNF_Opto\NDNF_06\070124_NDNF_06\Baseline_noopto\baselineMask_1.mat')

% frames_rando = sbxgrabframe('Z:\Scanbox\Data\NDNF\NDNF_04\062224_NDNF_04\062224_NDNF_04_000_003\062224_NDNF_04_000_003', 1, 100);
% frames_rando = squeeze(frames_rando); 

fft_tts = fft2((baseline.RefImage(:, 100:end)));
ref_img = baseline.RefImage(:, 100:end);

output = nan(size(frames_tmp, 3), 4);
for i = 1:size(frames_tmp, 3)

    reg = dftregistration(fft_tts, fft2(frames_tmp(:, 100:end, i)));
    newImg = imtranslate(frames_tmp(:, 100:end, i), reg(3:4));

    tic
    output(i,1) = immse(newImg(100:300, 100:300), ref_img(100:300, 100:300));
    toc

%     tic
%     output(i,1) = multissim(newImg(100:300, 100:300), ref_img(100:300, 100:300));
%     toc

%     tic
%     output(i,:) = ssim(frames_tmp(:, 100:end, i), baseline.RefImage(:, 100:end));
%     toc

%     tic
%     C = corrcoef(fft_tts, newImg);
%     output(i,:) = C(2);
%     toc

%     tic
%     output(i,:) = dftregistration(fft_tts, fft2(frames_tmp(:, 100:end, i)));
%     toc


%     tic
%     [output(i,1),output(i,2)] = fftalign(baseline.RefImage, frames_tmp(:, :, i));
%     toc

%     tic
%     output(i,1) = max(xcorr(reshape(baseline.RefImage(:, 100:end), [], 1), reshape(frames_tmp(:, 100:end, i), [], 1)));
%     toc

%     pause()
end
%frames = frames_tmp(:,92:end, 1:2:end);
figure(1)
plot(output(1:2:end, 1), 'r')
hold on
plot(output(2:2:end, 1), 'b')
hold off

figure(2)
plot(output(2:2:end, 1))

figure(3)
plot(output(1:2:end, 2))

figure(4)
plot(output(2:2:end, 2))

% figure(5)
% plot(output(1:2:end, 3))
% 
% figure(6)
% plot(output(2:2:end, 3))