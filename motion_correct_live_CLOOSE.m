%% MOTION CORRECTION SCRIPT
close all
clear all

% fname  = 'Z:\Scanbox\Data\rbp6_3\20210707_rbp6_3\20210707_rbp6_3_000_004'; 

Directory = 'Z:\Scanbox\Data';
cases = {'rbp6_3'; 'rbp10'; 'rbp11'; 'rbp12'; 'rbp13'; 'rbp16'};

% MC_values = cell(6, 6);
% 
% x_y_coord = cell(6, 6);

load('Z:\CLOOSE benchmark tests\MC\MC_values_new_method.mat')
% load('Z:\CLOOSE benchmark tests\MC\x_y_coord.mat')

for icase  = 1%5:6%numel(cases)
    listingDays  = dir([Directory '\' char(cases(icase))]);
    listingDays  = listingDays(~ismember({listingDays.name},{'.','..', 'matchedRoi', 'hide'}));
    for idays = 5%numel(listingDays)
        listingRecs = dir(fullfile([Directory '\' char(cases(icase)) '\'...
            char(listingDays(idays).name)], '*.sbx'));
        listingRecsT = struct2table(listingRecs);
        sortedT = sortrows(listingRecsT, 'name');
        sortedRecs = table2struct(sortedT);
        sortedRecs=sortedRecs(~ismember({sortedRecs.name},{'.','..'}));
        
        dldr = [Directory '\' char(cases(icase)) '\' char(listingDays(idays).name) '\'];
        
        disp(dldr);
        
        
        frames_tmp = sbxgrabframe([dldr sortedRecs(1).name(1:end-4)], 1, -1);
        frames_tmp = squeeze(frames_tmp); frames = frames_tmp(:,92:end, 1:2:end);
        
        load([dldr 'suite2p\plane1\Fall.mat'])
        fixed = ops.refImg;

        imshow(imadjust(fixed))

        pause()

        answer = questdlg('Use this image for Motion Correction?', ...
            'Options', ...
            'Yes','No', 'Yes');
        switch answer
            case 'Yes'
                figure(1)
                imshow(imadjust(fixed))

                [xcoord, ycoord] = ginput(4);
                x_y_coord{icase, idays} =  [xcoord, ycoord] ;
                vx = nan(4, 1);
                vy = nan(4, 1);
                for ixx = 1:size(xcoord)
                    vx(ixx, [1, 2]) = [round(xcoord(ixx) - 99), round(xcoord(ixx) + 99)];
                    vy(ixx, [1, 2]) = [round(ycoord(ixx) - 64), round(ycoord(ixx) + 63)];
                end

                if any(vx(:, 1) < 100)
                    idx = find(vx(:, 1) < 100);
                    vx(idx, 1) = 100; vx(idx, 2) = 298;
                end
                if any(vx(:, 2) > 794)
                    idx = find(vx(:, 2) > 794);
                    vx(idx, 1) = 596; vx(idx, 2) = 794;
                end
                if any(vy(:, 1) < 1)
                    idx = find(vy(:, 1) < 1);
                    vy(idx, 1) = 1; vy(idx, 2) = 128;
                end
                if any(vy(:, 2) > 510)
                    idx = find(vy(:, 2) > 510);
                    vy(idx, 1) = 383; vy(idx, 2) = 510;
                end

        end
        close all

%         vy = round(size(fixed,2)/4);
%         vx = round(size(fixed,1)/4);

%         fft_tts = nan(vx, vy, 4);
%         vx1 = [vx vx*2-1; vx vx*2-1; vx*2 vx*3-1; vx*2 vx*3-1];
%         vy1 = [vy vy*2-1; vy*2 vy*3-1; vy vy*2-1; vy*2 vy*3-1];

        fft_tts = nan(vy(1, 2) - vy(1, 1) + 1, vx(1, 2) - vx(1, 1) + 1, 3);
        % vx1 = [vx vx*2-1; vx vx*2-1; vx*2 vx*3-1; vx*2 vx*3-1];
        % vy1 = [vy vy*2-1; vy*2 vy*3-1; vy vy*2-1; vy*2 vy*3-1];

        vx1 = vy; %[vx vx*3-1; vx vx*3-1];
        vy1 = vx; %[vy vy*2-1; vy*2 vy*3-1];
        for i = 1:size(vx1,1)
            fft_tts(:,:,i) = fft2(fixed(vx1(i,1):vx1(i,2), vy1(i,1):vy1(i,2)));
        end


        output = nan(4, 4, size(frames, 3));
        xsh = nan(1, size(frames, 3));
        ysh = nan(1, size(frames, 3));
        for kk = 1:size(frames, 3)
            for i = 1:size(fft_tts, 3)
%                 output(i,:, kk) = dftregistration(fft_tts(:,:,i), fft2(double(frames(vx1(i,1):vx1(i,2), vy1(i,1):vy1(i,2), kk))));
                figure(i)
                imshow(imadjust(fixed(vx1(i,1):vx1(i,2), vy1(i,1):vy1(i,2), 1)))
            end
            xsh(kk) = round(nanmedian(output(:,4, kk)));
            ysh(kk) = round(nanmedian(output(:,3, kk)));
        end
        MC_values{icase, idays} = [xsh; ysh; ops.xoff(1:numel(xsh)); ops.yoff(1:numel(xsh))];

        pause()
    end
end

save('Z:\CLOOSE benchmark tests\MC\MC_values_new_method.mat', 'MC_values', '-v7.3')
save('Z:\CLOOSE benchmark tests\MC\x_y_coord.mat', 'x_y_coord', '-v7.3')