function plot_fc_operant_gc(FF_operant_call,trialInfo_operant,operant_t,mvgcStruct_comb)
% cla
k = 1;
baselineT = -1.5;
colors = {'r','k'};
nBat = 4;
% [operant_caller_idxs,operant_listener_idxs,operant_caller_bat_num,operant_listener_bat_num] = deal(zeros(1,length(FF_operant_call)));

% for exp_k = 1:length(FF_operant_call)
%     x = cellfun(@(x) x.used_bat_nums,trialInfo_operant{exp_k},'un',0);
%     x = [x{:}];
%     batNums = str2double(trialInfo_operant{exp_k}{1}.batNums);
%     
%     operant_caller_bat_num(exp_k) = mode(x);
%     operant_listener_bat_num(exp_k) = setdiff(batNums,operant_caller_bat_num(exp_k));
%     operant_caller_idxs(exp_k) = find(batNums == operant_caller_bat_num(exp_k));
%     operant_listener_idxs(exp_k) = setdiff(1:length(batNums),operant_caller_idxs(exp_k));
% end
% 
% for operant_idx = {operant_listener_idxs,operant_caller_idxs}
%     X = cellfun(@(x,idx) reshape(x(:,idx,:),[size(x,1) size(x,3)]),FF_operant_call,num2cell(operant_idx{1}),'un',0);
%     X = cat(1,X{:});
%     mu = nanmean(X(:,operant_t<baselineT,:),'all');
%     sigma = 1; %nanstd(X(:,operant_t<baselineT,:),[],'all');
%     %     X = nanmean(X,1);
%     XNorm = (X - mu)/sigma;
%     %     X_conf_int = nan(2,length(operant_t));
%     %     [X_conf_int(1,:),X_conf_int(2,:)] = mvgc_confint(0.05,X,8,50,1e3,4,4,0,'chi2');
%     %     X_conf_int = (X_conf_int - mu)/sigma;
%     %     boundedline(operant_t,XNorm',X_conf_int'/sqrt(2),[colors{k} '--'],'alpha')
%     boundedline(operant_t,nanmean(XNorm,1),nanstd(XNorm,[],1)/sqrt(length(operant_listener_idxs)),[colors{k} '--'],'alpha')
%     k = k + 1;
% end

% Element 1 of 1st dimension of FF_self (source_target_k = 1) is the
% influence of all other bats _ON_ a given bat (target). Element 2 of 1st dimension of
% FF_self (source_target_k = 2) is the influence _OF_ a given bat (source)
% on all other bats.
% target_source_batNums = {operant_listener_bat_num,operant_caller_bat_num};
for target_source_k = 1:2
%     batIdx = ismember(mvgcStruct_comb.batNums,target_source_batNums{target_source_k});
    FF_current = mvgcStruct_comb.FF_self(target_source_k,:,:);
    mu = squeeze(nanmean(FF_current(:,mvgcStruct_comb.time<baselineT,:,:),'all'));
    sigma = 1;% nanstd(FF_current(:,mvgcStruct_comb.time<baselineT,:,:),[],'all');
    nBat = size(FF_current,ndims(FF_current));
    FFNorm = (FF_current - mu)/sigma;
    boundedline(mvgcStruct_comb.time,squeeze(nanmean(FFNorm,[3 4])),squeeze(nanstd(FFNorm,[],[3 4]))./sqrt(nBat),colors{target_source_k},'alpha')
end
h = findobj(gca,'Type','Line');
legend(flipud(h),{'Operant Listener -> Caller','Operant Caller -> Listener','FC Listener -> Caller','FC Caller -> Listener'})
xlim([-2.2 2.2])
axis square
set(gca,'FontSize',25)
legend box off
xlabel('Time (s)')
ylabel('G-causality (baseline normalized)')