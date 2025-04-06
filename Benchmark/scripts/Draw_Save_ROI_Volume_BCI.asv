close all
clear all

load('Z:\BCI\Data\NDNF_Opto\NDNF_06\062424_NDNF_06\Baseline\baselineMask_2.mat')
dband = 80;
dband2 = 40;
frames_tmp = sbxgrabframe('Z:\Scanbox\Data\NDNF\NDNF_06\062424_NDNF_06\062424_NDNF_06_000_006\062424_NDNF_06_000_006', 1, 24000);
frames_tmp = squeeze(frames_tmp); frames_pl1 = frames_tmp(:,:, 1:2:end); frames_pl0 = frames_tmp(:,:, 2:2:end);

for ifr = 1:size(frames_pl1, 3)
    img__1 = imtranslate(frames_pl1(:, :, ifr), baseline.ImTranslateXY(:, ifr));
    frames_pl1(:, :, ifr) = img__1;

    img__0 = imtranslate(frames_pl0(:, :, ifr), baseline.ImTranslateXY(:, ifr));
    frames_pl0(:, :, ifr) = img__0;
end


figure(1); imshow(imadjust(imsharpen(uint16(mean(frames_pl1(dband2:end, dband:end, 1:50), 3)))));
% h = gca;
% for i = 1:10
%     ROI(i) = drawpolygon(h, 'LineWidth', 1, 'Deletable', true, ...
%         'Label', num2str(i), 'LabelVisible', 'hover', 'Color', 'y', ...
%         'FaceAlpha', 0);
% end
% Roi_plane = repmat(1, numel(ROI), 1);
% nR = numel(ROI);

figure(2); imshow(imadjust(imsharpen(uint16(mean(frames_pl0(dband2:end, dband:end, 1:50), 3)))));
% h = gca;
% for i = 1:10
%     nR = nR + 1;
%     ROI(nR) = drawpolygon(h, 'LineWidth', 1, 'Deletable', true, ...
%         'Label', num2str(i), 'LabelVisible', 'hover', 'Color', 'y', ...
%         'FaceAlpha', 0);
% end
% Roi_plane = cat(1, Roi_plane, repmat(2, numel(ROI), 1));

load('Z:\CLOOSE_benchmark_tests\Volume_bci\ROIS_062424_NDNF_06_000_006.mat')
load('Z:\CLOOSE_benchmark_tests\Volume_bci\ROI_plane_062424_NDNF_06_000_006.mat')
for i = 1:numel(ROI)
    if Roi_plane(i) == 1
        h = figure(1);
    elseif Roi_plane(i) == 2
        h = figure(2);
    end
    whatever = drawpolygon(gca, 'Position', cat(2, ROI(1,i).Position(:, 1) -80, ROI(1,i).Position(:, 2) -40),...
        'LineWidth', 1, 'Deletable', true, ...
        'Label', num2str(i), 'LabelVisible', 'hover', 'Color', 'y', ...
        'FaceAlpha', 0);
end

% save('Z:\CLOOSE_benchmark_tests\Volume_bci\ROIS_062424_NDNF_06_000_006.mat', 'ROI');
% save('Z:\CLOOSE_benchmark_tests\Volume_bci\ROI_plane_062424_NDNF_06_000_006.mat', 'Roi_plane');


pixelsX = 796;
LinesY = 512;
for iroi = 1:length(ROI)
    tb = nan(LinesY, pixelsX);
    for ix = 1:LinesY
        yq = repmat(ix, pixelsX,1);
        xq = [1:pixelsX]';
        tb(ix,1:pixelsX) = inpolygon(xq, yq, ROI(1,iroi).Position(:,1),...
        ROI(1,iroi).Position(:,2));
    end
    mask{iroi,1} = logical(tb);
end
roiMask = mask;

roi_one_frame = nan(size(roiMask, 1), 12000);
for i = 1:size(roiMask, 1)
    if Roi_plane(i) == 1
        for n = 1:size(frames_pl1, 3)
            frm_tmp  = frames_pl1(:, :, n);
            roi_one_frame(i,n) = double(mean(frm_tmp(roiMask{i, 1})));
        end
    elseif Roi_plane(i) == 2
        for n = 1:size(frames_pl1, 3)
            frm_tmp  = frames_pl0(:, :, n);
            roi_one_frame(i,n) = double(mean(frm_tmp(roiMask{i, 1})));
        end
    end
end
% 
% figure(100)
% plot(DF(2, :))

DF = calculateDf_smart(movmean(roi_one_frame, 10, 2),1800,10);
save('Z:\CLOOSE_benchmark_tests\Volume_bci\DF_062424_NDNF_06_000_006.mat', 'DF');
