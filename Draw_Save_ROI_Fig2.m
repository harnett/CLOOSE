close all
clear all

load('Z:\BCI\Data\NDNF_Opto\NDNF_07\062624_NDNF_07\Baseline_noopto\baselineMask_1.mat')

frames_tmp = sbxgrabframe('Z:\Scanbox\Data\NDNF\NDNF_07\062624_NDNF_07\062624_NDNF_07_000_009\062624_NDNF_07_000_009', 1, 24000);
frames_tmp = squeeze(frames_tmp); frames = frames_tmp(:,:, 1:2:end);



% writerObj = VideoWriter('spineAct.avi');
% writerObj.FrameRate = 60;
% open(writerObj);
for ifr = 1:size(frames, 3)
    img__ = imtranslate(frames(:, :, ifr), baseline.ImTranslateXY(:, ifr));
    frames(:, :, ifr) = img__;
%     writeVideo(writerObj, uint8(frames(:,:,ifr)));
end
% close(writerObj);
% implay('spineAct.avi');

figure(1); imshow(imadjust(baseline.RefImage));

h = gca;

for i = 1:31
    ROI(i) = drawpolygon(h, 'LineWidth', 1, 'Deletable', true, ...
        'Label', num2str(i), 'LabelVisible', 'hover', 'Color', 'y', ...
        'FaceAlpha', 0);
end

save('Z:\CLOOSE benchmark tests\Various_BCI_approaches\ROIS_062624_NDNF_07_000_009.mat', 'ROI');


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
    for n = 1:size(frames, 3)
        frm_tmp  = frames(:, :, n);
        roi_one_frame(i,n) = double(mean(frm_tmp(roiMask{i, 1})));
    end
end
% 
% figure(100)
% plot(DF(2, :))

DF = calculateDf_smart(movmean(roi_one_frame, 10, 2),1800,10);
save('Z:\CLOOSE benchmark tests\Various_BCI_approaches\DF_062624_NDNF_07_000_009.mat', 'DF');
