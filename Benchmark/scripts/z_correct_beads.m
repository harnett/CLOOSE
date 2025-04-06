%% z-correct
close all
clear all

rng(10)

fl = 'Z:\CLOOSE_benchmark_tests\z_correction\test_beads_z_correction\test_beads_z_correction_000_002\test_beads_z_correction_000_002';

dead_band = 80;
cXY = [30, 210; 70, 198];
z_nplanes = 21;

frames_tmp = sbxgrabframe(fl, 1, -1);
frames_tmp = squeeze(frames_tmp); 
frames_tmp = frames_tmp (:, dead_band:end, :);
% Plot for sanity check
% for ifr = 1:size(frames_tmp, 3)
%     figure(1)
%     imshow(imadjust(frames_tmp(:, :, ifr)))
%     pause()
% end
ref_img = uint16(mean(frames_tmp(:, :, 11:z_nplanes:210), 3));

ref_img_all = uint16(zeros(size(frames_tmp, 1), size(frames_tmp, 2), z_nplanes));
for ipl = 1:z_nplanes
    
    if ipl == 21
        ref_img_all(:, :, ipl)  = uint16(mean(frames_tmp(:, :, 22:21:210), 3));
    else
        ref_img_all(:, :, ipl)  = uint16(mean(frames_tmp(:, :, ipl:21:210), 3));
    end

    figure(ipl)
    imshow((ref_img_all(:, :, ipl)))

end


figure(100)
imshow(ref_img)
hold off

fft_tts = fft2(ref_img);

rand_vect = randsample(size(frames_tmp, 3),size(frames_tmp, 3));
OrigPlane = rem(rand_vect, 21);

output = nan(size(frames_tmp, 3), 21);

times_processing = nan(size(frames_tmp, 3), 21);

output_rand = nan(size(frames_tmp, 3), 2);
for k = 1:size(frames_tmp, 3)

%     if rem(k, 21) == 0 || rem(k, 21) > 16 
%         continue
%     end
% 
    reg = dftregistration(fft_tts(cXY(2,1):cXY(2, 2), cXY(1,1):cXY(1, 2)), fft2(frames_tmp(cXY(2,1):cXY(2, 2), cXY(1,1):cXY(1, 2), k)));
%     newImg = imtranslate(frames_tmp(:, :, k), [reg(4), reg(3)]);
    newImg = frames_tmp(:, :, k);
    
    for iz = 1:z_nplanes
        tic
%         A = newImg(cXY(2,1):cXY(2, 2), cXY(1,1):cXY(1, 2));
%         B = ref_img_all(cXY(2,1):cXY(2, 2), cXY(1,1):cXY(1, 2), iz);
%         output(k,iz) = 1- corr(double(reshape(A, [], 1)), double(reshape(B, [], 1)));

        output(k,iz) = immse(newImg(cXY(2,1):cXY(2, 2), cXY(1,1):cXY(1, 2)), ref_img_all(cXY(2,1):cXY(2, 2), cXY(1,1):cXY(1, 2), iz));
        CMMSE = immse(newImg(cXY(2,1):cXY(2, 2), cXY(1,1):cXY(1, 2)), ref_img_all(cXY(2,1):cXY(2, 2), cXY(1,1):cXY(1, 2), iz));
        output(k,iz) = CMMSE./mean(newImg(cXY(2,1):cXY(2, 2), cXY(1,1):cXY(1, 2)), 'all');


%         output(k,iz) = psnr(newImg(cXY(2,1):cXY(2, 2), cXY(1,1):cXY(1, 2)), ref_img_all(cXY(2,1):cXY(2, 2), cXY(1,1):cXY(1, 2), iz), "DataFormat","SSC");

%         tn = dftregistration(newImg(cXY(2,1):cXY(2, 2), cXY(1,1):cXY(1, 2)), ref_img_all(cXY(2,1):cXY(2, 2), cXY(1,1):cXY(1, 2), iz));
%         output(k,iz) = tn(1);

        times_processing(k,iz) = toc; 
        
%         figure(iz)
%         imshow(imabsdiff(newImg(cXY(2,1):cXY(2, 2), cXY(1,1):cXY(1, 2)), ref_img_all(cXY(2,1):cXY(2, 2), cXY(1,1):cXY(1, 2), iz)))
%         
    end
% 
%     figure(100)
%     plot([output(k, 21),output(k, 1:20)], 'k')
%     title(num2str(k))
%     box off
%     xlim([0 22]) 
%     xticks(1:5:21)
%     set(gca,'TickDir','out'); % The only other option is 'in'
%     xline(11, 'r')
% 
% 
%     hold off
% 
% % 
%     pause()

%     rnd_idx = rand_vect(k);
%     reg = dftregistration(fft_tts, fft2(frames_tmp(:, :, rnd_idx)));
%     newImg = imtranslate(frames_tmp(:, 100:end, rnd_idx), reg(3:4));
%     
%     tic
%     output_rand(k,1) = immse(newImg(cXY(2,1):cXY(2, 2), cXY(1,1):cXY(1, 2)), ref_img(cXY(2,1):cXY(2, 2), cXY(1,1):cXY(1, 2)));
%     output_rand(k,2) = toc;

%     tic
%     output(i,1) = multissim(newImg(100:300, 100:300), ref_img(100:300, 100:300));
%     toc

%     tic
%     output(i,:) =  (frames_tmp(:, 100:end, i), baseline.RefImage(:, 100:end));
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

[mnn, idx] = min(output, [], 2);
cxc = repmat([1:20, 0]', 1000, 1); cxc = cxc(1:numel(idx));

cxc2 = find(cxc <22); 
cxc3 = cxc(cxc2);

figure(500)
% histogram(idx - rem(1:size(frames_tmp, 3), 21));
% histogram(idx(cxc2) - cxc(cxc2), [-8.5:1:8.5], 'Normalization', 'probability')
histogram(rem(idx(cxc2), 21) - cxc(cxc2),  'Normalization', 'probability', 'FaceColor', 'k')
box off
set(gca,'TickDir','out'); % The only other option is 'in'
ylim([0 0.6])
yticks(0:0.3:0.6)
hold off


figure(501)
plot(rem(idx(cxc2), 21) - cxc(cxc2))

figure(750)
tm = nansum(times_processing,  2);
histogram(tm(cxc2).*1000, 0:0.25:7,  'Normalization', 'probability', 'FaceColor', 'k')
box off
set(gca,'TickDir','out'); % The only other option is 'in'
ylim([0 0.7])
yticks(0:0.35:0.7)
% xlim([0 5])
% xticks(0:1:5)
hold off

figure(1000)
plot(idx(cxc2) - cxc(cxc2))