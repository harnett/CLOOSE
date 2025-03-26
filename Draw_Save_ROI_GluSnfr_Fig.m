close all
clear all

glusnf = 1;
load('Y:\users\zilan662\Glusnfr_spine_metadata\glusnfr_spine_bci\061324_glusnfrspinebci5\Baseline\baselineMask_1.mat')

frames_tmp_ = sbxgrabframe('Y:\users\zilan662\Glusnfr_spine\glusnfrspinebci05\061324_glusnfrspinebci05\061324_glusnfrspinebci05_000_000\061324_glusnfrspinebci05_000_000', 1, 15000);
frames_tmp_ = squeeze(frames_tmp_); 

d_band = 80;
fold_fct = 6;
LinesY = size(frames_tmp_, 1);
pixelsX = 796;
fold = LinesY/fold_fct;
split_fr = nan(fold_fct, 2);
for i = 1:fold_fct
    split_fr(i, 1:2) = [fold*i - fold + 1, fold*i];
end

frames = uint16(zeros(fold, size(frames_tmp_, 2),size(frames_tmp_, 3)));
r_tmp_special = uint16(zeros(fold, 796, 12));
for ifr = 1:size(frames, 3)
    raw_Data = frames_tmp_(:, :, ifr);
    r_tmp = uint16(zeros(fold, 796, fold_fct));
    for ifc = 1:fold_fct %#ok<BDSCI>
        fuck_index = (fold*(ifc-1))+ 1 : fold*ifc;
        r_tmp(:, :, ifc) = raw_Data(fuck_index, :);
        if ifr == 2 
            figure(ifc + 10)
            imshow(imadjust(imsharpen(r_tmp(:, d_band:end, ifc), Radius=200,Amount=2)));
            r_tmp_special(:, :, ifc) = raw_Data(fuck_index, :);
        end
        if ifr == 3
            figure(ifc + 100)
            imshow(imadjust(imsharpen(r_tmp(:, d_band:end, ifc), Radius=200,Amount=2)));
            r_tmp_special(:, :, ifc+6) = raw_Data(fuck_index, :);
        end
    end
    frames(:, :, ifr) = uint16(mean(r_tmp, 3));
end
figure(700)
imshow(imadjust(imsharpen(uint16(mean(r_tmp_special(:, 80:end, 4:11), 3)), Radius=200,Amount=2)))

for ifr = 1:size(frames, 3)
    img__ = imtranslate(frames(:, :, ifr), baseline.ImTranslateXY(:, ifr));
    frames(:, :, ifr) = img__;
end

figure(1); imshow(imadjust(imsharpen(baseline.RefImage(:, d_band:end), Radius=4,Amount=1)));
% 
% for i = 1:10
%     figure(i)
%     imshow(imadjust(imsharpen(frames(:, d_band:end, i), Radius=25,Amount=1)));
% end

% pause()

h = gca;

for i = 1:3
    ROI(i) = drawpolygon(h, 'LineWidth', 1, 'Deletable', true, ...
        'Label', num2str(i), 'LabelVisible', 'hover', 'Color', 'y', ...
        'FaceAlpha', 0);
end

save('Z:\CLOOSE_benchmark_tests\Spine_imaging\ROIS_gcampspinebci04_000_006.mat', 'ROI');


for iroi = 1:length(ROI)
    tb = nan(fold, pixelsX);
    for ix = 1:fold
        yq = repmat(ix, pixelsX,1);
        xq = [1:pixelsX]';
        tb(ix,1:pixelsX) = inpolygon(xq, yq, ROI(1,iroi).Position(:,1),...
        ROI(1,iroi).Position(:,2));
    end
    mask{iroi,1} = logical(tb);
end
roiMask = mask;

roi_one_frame = nan(size(roiMask, 1), 15000);
for i = 1:size(roiMask, 1)
    for n = 1:size(frames, 3)
        frm_tmp  = frames(:, :, n);
        roi_one_frame(i,n) = double(mean(frm_tmp(roiMask{i, 1})));
    end
end
% 
% figure(100)
% % plot(DF(2, :))
% [~, X_sub] = robustfit(roi_one_frame(3, :), roi_one_frame(1, :));
% [~, Y_sub] = robustfit(roi_one_frame(3, :), roi_one_frame(2, :));
% roi_one_frame(1, :) = X_sub.resid;
% roi_one_frame(2, :) = Y_sub.resid;

DF = calculateDf_smart(movmean(roi_one_frame, 10, 2),1800,10);
save('Z:\CLOOSE_benchmark_tests\Spine_imaging\DF_gcampspinebci04_000_006.mat', 'DF');
