close all
clear all

load('Z:\CLOOSE_benchmark_tests\MC\MC_values_new_method.mat')


dispX = cell(6,6);
dispY = cell(6, 6);

xc = nan(6);
yc = nan(6);

for k = 1:6
    for kk = 1:5
        k
        kk

        dispX{k, kk} = MC_values{k,kk}(1, :) -  (-MC_values{k,kk}(3, :)); 
        
        figure(1)
        plot(dispX{k, kk})
        hold off

        dispY{k, kk} = MC_values{k,kk}(2, :) - (-MC_values{k,kk}(4, :));

        figure(2)
        plot(dispY{k, kk})
        hold off

        cR = corrcoef(double(MC_values{k,kk}(1, :)), double(-(MC_values{k,kk}(3, :))));
        cR = cR(2);
        title(cR)
        xc(k,kk) = cR;

        cR = corrcoef(double(MC_values{k,kk}(2, :)), double(-(MC_values{k,kk}(4, :))));
        cR = cR(2);
        title(cR)
        yc(k,kk) = cR;


        figure(1)
        scatter(MC_values{k,kk}(1, :), -(MC_values{k,kk}(3, :)), 'r')
        cR = corrcoef(double(MC_values{k,kk}(1, :)), double(-(MC_values{k,kk}(3, :))));
        cR = cR(2);
        title(cR)
        xc(k,kk) = cR;
        mr = refline([1,0]);
        mr.Color = 'r';
        hold off
        PxlShX{k,kk} = double(MC_values{k,kk}(1, :)) - double(-(MC_values{k,kk}(3, :)));
        UmX(k,kk) = nanmean(abs(PxlShX{k,kk})); stdX(k,kk) = nanstd(abs(PxlShX{k,kk}));
        UmX_99prctile(k,kk) = prctile(abs(PxlShX{k,kk}), 99);
        
        mtxX = nan(max(MC_values{k,kk}(3, :))- min(MC_values{k,kk}(3, :)), max(MC_values{k,kk}(1, :))- min(MC_values{k,kk}(1, :)));
        c = 0;
        for ix = double(min(MC_values{k,kk}(3, :)): max(MC_values{k,kk}(3, :)))
            c = c + 1;
            cc = 0;
            for ixx = double(min(MC_values{k,kk}(1, :)): max(MC_values{k,kk}(1, :)))
                cc = cc + 1;
                mtxX(cc, c) = sum(MC_values{k,kk}(3, :) == ix & MC_values{k,kk}(1, :) == ixx) / sum(MC_values{k,kk}(3, :) == ix); 
            end
        end
        figure(7)
        heatmap(mtxX, 'CellLabelColor','none')
        colormap parula;
        title([min(MC_values{k,kk}(3, :)); max(MC_values{k,kk}(3, :))])
        
        
        figure(3)
        histogram(abs(PxlShX{k,kk}))
        hold off
        
        figure(2)
        scatter(MC_values{k,kk}(2, :), -(MC_values{k,kk}(4, :)), 'g')
%         heatscatter(MC_values{k,kk}(2, :), -(MC_values{k,kk}(4, :)))
        cR = corrcoef(double(MC_values{k,kk}(2, :)), double(-(MC_values{k,kk}(4, :))));
        cR = cR(2);
        title(cR)
        yc(k,kk) = cR;
        mr = refline([1,0]);
        mr.Color = 'r';
        hold off
        PxlShY{k,kk} = double(MC_values{k,kk}(2, :)) - double(-(MC_values{k,kk}(4, :)));
        UmY(k,kk) = nanmean(abs(PxlShY{k,kk})); stdY(k,kk) = nanstd(abs(PxlShY{k,kk}));
        UmY_99prctile(k,kk) = prctile(abs(PxlShY{k,kk}), 99);
        
        figure(4)
        histogram(abs(PxlShY{k,kk}))
        hold off
        
        % figure(5)
        % violin(PxlShY{k,kk}')
        % hold off
        % 
        % figure(6)
        % violin(PxlShX{k,kk}')
        % hold off
        
        
        mtxY = nan(max(MC_values{k,kk}(4, :))- min(MC_values{k,kk}(4, :)), max(MC_values{k,kk}(2, :))- min(MC_values{k,kk}(2, :)));
        c = 0;
        for ix = double(min(MC_values{k,kk}(4, :)): max(MC_values{k,kk}(4, :)))
            c = c + 1;
            cc = 0;
            for ixx = double(min(MC_values{k,kk}(2, :)): max(MC_values{k,kk}(2, :)))
                cc = cc + 1;
                mtxY(cc, c) = sum(MC_values{k,kk}(4, :) == ix & MC_values{k,kk}(2, :) == ixx) / sum(MC_values{k,kk}(4, :) == ix); 
            end
        end
        figure(8)
        heatmap(mtxY, 'CellLabelColor','none')
        colormap parula;
        title([min(MC_values{k,kk}(4, :)); max(MC_values{k,kk}(4, :))])
        
%         pause()

    end
end



CI99X = UmX + 2.57*stdX;
CI99Y = UmY + 2.57*stdY;

figure(10)
boxplot([reshape(xc, [], 1), reshape(yc, [], 1)])
box off
hold off

figure(11)
violin([reshape(UmX_99prctile, [], 1), reshape(UmY_99prctile, [], 1)], 'facecolor', 'w', 'edgecolor', 'k')
box off
hold off

figure(12)
histogram(reshape(UmX_99prctile, [], 1),'Normalization', 'probability', 'facecolor', 'w', 'edgecolor', 'k')
box off
xticks(0:10)
xlim([0 5.75])
yticks(0:0.2:1)
set(gca,'TickDir','out'); % The only other option is 'in'
box off
hold off

figure(13)
histogram(reshape(UmY_99prctile, [], 1),'Normalization', 'probability', 'facecolor', 'w', 'edgecolor', 'k')
xticks(0:10)
box off
xlim([0 5.75])
yticks(0:0.2:1)
set(gca,'TickDir','out'); % The only other option is 'in'
hold off

umX = mean(UmX, "all");
umY = mean(UmY, "all");

SEM1 = (std(UmX, [], 'all'))./sqrt(numel(umX, 2));
SEM2 = (std(UmY, [], 'all'))./sqrt(numel(umY, 2));

figure(1)
bar([umX, umY], 'FaceColor', 'w') 
hold on
errorbar([umX, umY], [SEM1, SEM2] , '.', 'Color', 'k') 
box off
set(gca,'TickDir','out'); % The only other option is 'in'
ylim([0 0.25])
yticks([0:0.25:0.25])
hold off


