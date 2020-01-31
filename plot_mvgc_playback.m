function plot_mvgc_playback(FF_all,FF_playback,t)
[~,idx] = inRange(t,[-4 -1.5]);
cla
hold on
colors = {'r','k'};
k = 1;
for ff = {FF_all{1},FF_playback}
    ff = ff{1};
    mu = nanmean(ff(:,:,idx),'all');
    sigma = nanstd(ff(:,:,idx),[],'all');
    boundedline(t,squeeze(nanmean((ff-mu)/sigma,[1 2])),squeeze(nanstd((ff-mu)/sigma,[],[1 2]))/sqrt(size(ff,1)),colors{k},'alpha')
    k = k + 1;
end

set(gca,'FontSize',25)
xlabel('Time (s)')
ylabel('GC (baseline normalized)')
h = findobj(gca,'Type','Line');
legend(h,{'Communication','Playback'});
legend box off
axis square
xlim([-2.2 2.2])

end

