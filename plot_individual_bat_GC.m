function plot_individual_bat_GC(mvgcData,plotType)


[~,idx] = inRange(mvgcData.mvgcStruct_comb.time,[-1 1]);
[~,baseIdx] = inRange(mvgcData.mvgcStruct_comb.time,[-4 -1]);
mu = squeeze(nanmean(mvgcData.mvgcStruct_comb.FF_self(:,baseIdx,:,:),[2 3]));

switch plotType
    case 'batDist'
        titleStrs = {'Listener -> Caller','Caller -> Listener'};
        for source_k = 1:2
            subplot(1,2,source_k)
            boxplot(squeeze(nanmean(mvgcData.mvgcStruct_comb.FF_self(source_k,idx,:,:),2)) - mu(source_k,:))
            h = gca;
            h.XTickLabel = strsplit(num2str(mvgcData.mvgcStruct_comb.batNums));
            ylim([-7e-3 14e-3])
            set(gca,'FontSize',25)
            h.XTickLabelRotation = 25;
            box off
            axis square
            title(titleStrs{source_k});
        end
        
    case 'linearFit'
        
        x = squeeze(nanmean(mvgcData.mvgcStruct_comb.FF_self(1,idx,:,:),[2 3]))' - mu(1,:);
        y = squeeze(nanmean(mvgcData.mvgcStruct_comb.FF_self(2,idx,:,:),[2 3]))' - mu(2,:);
        L = fitlm(x,y);
        L.plot
        xlabel('Listener -> Caller')
        ylabel('Caller -> Listener')
        set(gca,'FontSize',25)
        axis square
        legend off
        title('')
        box off
        text(0,8e-3,sprintf('p = %0.2f \nR^{2} = %0.2f',L.Coefficients.pValue(1),L.Rsquared.Ordinary),'FontSize',25)
end

end