function plot_mvgc_by_operant_pair(trialInfo,mvgcStruct)

nPair = length(trialInfo);
operantPair = nan(nPair,2);
operantVocalizer = nan(nPair,1);
for pair_k = 1:nPair
    operantPair(pair_k,:) = cellfun(@str2num,trialInfo{pair_k}{1}.batNums);
    x = cellfun(@(x) x.used_bat_nums,trialInfo{pair_k},'un',0);
    x = [x{:}];
    operantVocalizer(pair_k) = mode(x);
end

operantPair = {operantPair,operantPair'};
operantVocalizer = {operantVocalizer,flipud(operantVocalizer)};
nPair_type = length(operantPair);

expIdx = strcmp({mvgcStruct.expType},'operant');
nT = length(mvgcStruct(expIdx).time);
FF_by_operant_pair = zeros(2,nT,nPair,nPair_type);

for pair_k = 1:nPair
    for pair_type_k = 1:nPair
        included_bat_idx = mvgcStruct(expIdx).batNums == operantVocalizer{pair_type_k}(pair_k);
        ff_tmp = mvgcStruct(expIdx).FF{included_bat_idx};
        
        calling_bat_idx =  mvgcStruct(expIdx).exp_bat_nums == operantVocalizer{pair_type_k}(pair_k);
        listening_bat_idx =  mvgcStruct(expIdx).exp_bat_nums == setdiff(operantPair{pair_type_k}(pair_k,:),operantVocalizer{pair_type_k}(pair_k));
        
        FF_by_operant_pair(1,:,pair_k,pair_type_k) = ff_tmp(calling_bat_idx,listening_bat_idx,:);
        FF_by_operant_pair(2,:,pair_k,pair_type_k) = ff_tmp(listening_bat_idx,calling_bat_idx,:);
    end
end

[~,idx] = inRange(mvgcStruct(1).time,[-Inf -1.5]);
colors = {'r','k';'g','b'};
mu = squeeze(mean(FF_by_operant_pair(:,idx,:,:),[2 3]));
sigma = squeeze(std(FF_by_operant_pair(:,idx,:,:),[],[2 3]));

for k = 1:2
    for j = 1:2
        X = (FF_by_operant_pair(j,:,:,k) - mu(j,k))/sigma(j,k);
        boundedline(t,squeeze(mean(X,3)),std(X,[],3)/sqrt(nPair),colors{j,k},'alpha')
    end
end