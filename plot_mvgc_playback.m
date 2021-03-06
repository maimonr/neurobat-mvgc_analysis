function plot_mvgc_playback(mvgcStruct,FF_playback)
[~,idx] = inRange(mvgcStruct.time,[-4 -1.5]);
cla
hold on
colors = {'r','k'};
k = 1;
nBat = length(mvgcStruct.FF);
nT = length(mvgcStruct(1).time);
ff_inter_listener = nan(nBat,1,nT);
exp_bat_k = 1;
for bat_k = 1:length(mvgcStruct.FF)
    ff = mvgcStruct.FF{bat_k};
    nBat = size(ff,1);
    non_bat_idx = setdiff(1:nBat,exp_bat_k);
    ff = ff(non_bat_idx,non_bat_idx,:);
    ff_inter_listener(bat_k,1,:) = squeeze(nanmean(ff,[1 2]))';
    exp_bat_k = exp_bat_k + 1;
    if exp_bat_k == 4
        exp_bat_k = 1;
    end
end
for ff = {ff_inter_listener,FF_playback}
    ff = ff{1};
    mu = nanmean(ff(:,:,idx),'all');
    sigma = 1; %nanstd(ff(:,:,idx),[],'all');
    boundedline(mvgcStruct.time,squeeze(nanmean((ff-mu)/sigma,[1 2])),squeeze(nanstd((ff-mu)/sigma,[],[1 2]))/sqrt(size(ff,1)),colors{k},'alpha')
    k = k + 1;
end

set(gca,'FontSize',25)
xlabel('Time (s)')
ylabel('GC (baseline subtracted)')
h = findobj(gca,'Type','Line');
legend(h,{'Communication','Playback'});
legend box off
axis square
xlim([-2.2 2.2])

end

