tt2 = tt;
tt2(1) = nan;
tt2 = tt2*1000;
edges = linspace(prctile(tt2, 1), prctile(tt2, 99), 100);

figure(1)
histogram(tt2, edges, FaceColor='k')
box off
set(gca,'TickDir','out'); % The only other option is 'in'
ylabel('Count')
xlabel('ms')
title('uint16 512x796 internal')
yticks([0:500:1000])
ylim([0 1000])                          
xticks([0:4:8])
xlim([0 8])

