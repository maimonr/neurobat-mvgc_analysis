function plot_dGC_trial_shuffle(mvgcData)

ff_self_rep = nan(size(mvgcStruct_trial_shuffle.FF_self));
timeWin = [-0.5 0.5];
[~,idx] = inRange(mvgcStruct_trial_shuffle.time,timeWin);
nBoot = size(ff_self_rep,3);
nBat = size(mvgcData.mvgcStruct_comb.FF_self,3);
for bat_k = 1:nBat
    ff_self_rep(:,:,:,bat_k) = repmat(mvgcData.mvgcStruct_comb.FF_self(:,:,bat_k),[1 1 nBoot]);
end
dGC = (ff_self_rep(:,idx,:,:) - mvgcStruct_trial_shuffle.FF_self(:,idx,:,:));

cla
hold on
histogram(dGC,'Normalization','probability','FaceColor','k','FaceAlpha',1)
plot([0 0],get(gca,'YLim'),'r--');
p = 1 - sum(dGC>0,'all')/numel(dGC);
t = text(-0.01,0.035,sprintf('p = %0.2f',p),'FontSize',25);
box off
axis square
ylabel('Probability')
xlabel('\DeltaGC (actual - trial shuffle)')
h = gca;
h.YTick = [];
xlim([-0.015 0.015])
h.XTick = [-0.015 0 0.015];
