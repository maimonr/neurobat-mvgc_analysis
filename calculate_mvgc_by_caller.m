function [mvgcStruct, mvgcStruct_comb] = calculate_mvgc_by_caller(lfpPower,trialInfo,varargin)

pnames = {'timeWin','stepSize','selectBatType','modelOrder'};
dflts  = {[-4 4],0.25,'individual','AIC'};
[timeWin,stepSize,select_bat_type,modelOrder] = internal.stats.parseArgs(pnames,dflts,varargin{:});

mvgcStruct = struct('FF',[],'FF_conf',[],'selectBat',[],'expType',[],'time',[]);
selectBat_strs = 'includeBat';
expStrs = {'adult','operant'};

for exp_k = 1:length(lfpPower)
    idx = find(~cellfun(@isempty,trialInfo{exp_k}),1);
    exp_bat_nums = cellfun(@str2num,trialInfo{exp_k}{idx}.batNums);
    batNums = exp_bat_nums;
    switch select_bat_type
        case 'nonfamiliar'
            batNums = exp_bat_nums;
            trialInfo_tmp = [trialInfo{exp_k}{:}];
            all_used_bat_nums = [trialInfo_tmp.used_bat_nums];
            batNums = setdiff(all_used_bat_nums,batNums);
            batNums = batNums(~isnan(batNums));
            batCount = zeros(1,length(batNums));
            for select_bat_k = 1:length(batNums)
                batCount(select_bat_k) = sum(all_used_bat_nums == batNums(select_bat_k));
            end
            batNums = batNums(batCount > 50);
    end
    nBat = length(batNums);
    
    mvgcStruct(exp_k).expType = expStrs{exp_k};
    mvgcStruct(exp_k).selectBat = selectBat_strs;
    [mvgcStruct(exp_k).FF,mvgcStruct(exp_k).FF_conf] = deal(cell(1,nBat));
    
    for bat_k = 1:length(batNums)
        [mvgcStruct(exp_k).FF{bat_k},t,mvgcStruct(exp_k).FF_conf{bat_k}] = calculate_mvgc(lfpPower{exp_k},trialInfo{exp_k},...
            'trialSelection',selectBat_strs,'selectBat',batNums(bat_k),...
            'timeWin',timeWin,'stepSize',stepSize,'modelOrder',modelOrder);
    end
    mvgcStruct(exp_k).time = t;
    mvgcStruct(exp_k).batNums = batNums;
    mvgcStruct(exp_k).exp_bat_nums = exp_bat_nums;
end

%%
for exp_k = 1:length(mvgcStruct)
    nBat = length(mvgcStruct(exp_k).batNums);
    [FF_self,FF_self_conf] = deal(cell(1,nBat));
    
    for bat_k = 1:nBat
        ff_tmp = mvgcStruct(exp_k).FF{bat_k};
        ff_conf_tmp = mvgcStruct(exp_k).FF_conf{bat_k};
        switch select_bat_type
            case 'individual'
                FF_self{bat_k} = cat(1,ff_tmp(bat_k,:,:),permute(ff_tmp(:,bat_k,:),[2 1 3]));
                FF_self{bat_k} = FF_self{bat_k}(:,setdiff(1:nBat,bat_k),:);
                FF_self{bat_k} = squeeze(mean(FF_self{bat_k},2));
                
                FF_self_conf{bat_k} = cat(1,ff_conf_tmp(bat_k,:,:,:),permute(ff_conf_tmp(:,bat_k,:,:),[2 1 3 4]));
                FF_self_conf{bat_k} = FF_self_conf{bat_k}(:,setdiff(1:nBat,bat_k),:,:);
                FF_self_conf{bat_k} = squeeze(mean(FF_self_conf{bat_k},2));
            case 'nonfamiliar'
                FF_self{bat_k} = squeeze(nanmean(ff_tmp,[1 2]))';
                FF_self_conf{bat_k} = nanmean(ff_conf_tmp,[1 2]);
                FF_self_conf{bat_k} = reshape(FF_self_conf{bat_k},[1 size(FF_self_conf{bat_k},3) 2]);
        end
    end
    
    FF_self = cat(3,FF_self{:});
    FF_self_conf = cat(4,FF_self_conf{:});
    FF_self_conf = permute(FF_self_conf,[1 2 4 3]);
    
    [mvgcStruct(exp_k).FF_self,mvgcStruct(exp_k).FF_self_conf] = deal(FF_self,FF_self_conf);
end

%%

t_idx{1} = true(1,length(mvgcStruct(1).time));
t_idx{2} = false(1,length(mvgcStruct(2).time));
t_idx{2}(5:22) = true;
t = mvgcStruct(1).time;
batNums = [mvgcStruct(1).batNums,mvgcStruct(2).batNums];
exp_bat_nums = [mvgcStruct(1).exp_bat_nums,mvgcStruct(2).exp_bat_nums];

[FF,FF_conf,FF_self,FF_self_conf] = deal(cell(1,length(expStrs)));

for exp_k = 1:length(expStrs)
    
    FF{exp_k} = mvgcStruct(exp_k).FF;
    FF{exp_k} = cellfun(@(ff) ff(:,:,t_idx{exp_k}),FF{exp_k},'un',0);
    
    FF_conf{exp_k} = mvgcStruct(exp_k).FF_conf;
    FF_conf{exp_k} = cellfun(@(ff) ff(:,:,t_idx{exp_k},:),FF_conf{exp_k},'un',0);
    
    FF_self{exp_k} = mvgcStruct(exp_k).FF_self(:,t_idx{exp_k},:);
    FF_self_conf{exp_k} = mvgcStruct(exp_k).FF_self_conf(:,t_idx{exp_k},:,:);
    
end

FF = [FF{:}];
FF_conf = [FF_conf{:}];

FF_self = cat(3,FF_self{:});
FF_self_conf = cat(3,FF_self_conf{:});

mvgcStruct_comb = struct('FF',{FF},'FF_conf',{FF_conf},'FF_self',FF_self,'FF_self_conf',FF_self_conf,'selectBat',selectBat_strs,'time',t,'batNums',batNums,'exp_bat_nums',exp_bat_nums);

end
