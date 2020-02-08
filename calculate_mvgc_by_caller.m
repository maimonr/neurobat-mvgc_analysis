function mvgcStruct = calculate_mvgc_by_caller(lfpPower,trialInfo,varargin)

pnames = {'timeWin','stepSize','selectBatType','modelOrder','minCalls','acmaxlags','subSample','nSub'};
dflts  = {[-4 4],0.25,'individual','AIC',50,100,'none',10};
[timeWin,stepSize,select_bat_type,modelOrder,minCalls,acmaxlags,subSample,nSub] = internal.stats.parseArgs(pnames,dflts,varargin{:});

mvgcStruct = struct('FF',[],'FF_conf',[],'selectBat',[],'expType',[],'time',[]);
selectBat_strs = 'includeBat';
expStrs = {'adult','operant'};
tic;
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
            batNums = batNums(batCount > minCalls);
    end
    nBat = length(batNums);
    
    mvgcStruct(exp_k).expType = expStrs{exp_k};
    mvgcStruct(exp_k).selectBat = selectBat_strs;
    [mvgcStruct(exp_k).FF,mvgcStruct(exp_k).FF_conf] = deal(cell(1,nBat));
    
    for bat_k = 1:length(batNums)
        [mvgcStruct(exp_k).FF{bat_k},t,mvgcStruct(exp_k).FF_conf{bat_k}] = calculate_mvgc(lfpPower{exp_k},trialInfo{exp_k},...
            'trialSelection',selectBat_strs,'selectBat',batNums(bat_k),...
            'timeWin',timeWin,'stepSize',stepSize,'modelOrder',modelOrder,...
            'acmaxlags',acmaxlags,'subSample',subSample,'nSub',nSub);
        sprintf('%d bats out of %d done, %d/%d exps., %d sec elapsed',bat_k,length(batNums),exp_k,length(lfpPower),toc);
    end
    mvgcStruct(exp_k).time = t;
    mvgcStruct(exp_k).batNums = batNums;
    mvgcStruct(exp_k).exp_bat_nums = exp_bat_nums;
end
end
