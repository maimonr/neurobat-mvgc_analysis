function mvgc_over_time = calculate_mvgc_over_time(lfpPower,trialInfo,varargin)
pnames = {'timeWin','stepSize','modelOrder','binningMethod','minCalls'};
dflts  = {[-4 4],0.25,8,'equalCalls',50};
[timeWin,stepSize,modelOrder,binningMethod,minCalls] = internal.stats.parseArgs(pnames,dflts,varargin{:});
nExp = length(lfpPower);
for exp_k = 1:nExp
    X = cat(3,lfpPower{exp_k}{:});
    trialInfo_tmp = [trialInfo{exp_k}{:}];
    used_date_idx = ~cellfun(@isempty,lfpPower{exp_k});
    expDates = cellfun(@(expDate,lfppower) repmat(expDate,1,size(lfppower,3)),{trialInfo_tmp.expDate},lfpPower{exp_k}(used_date_idx),'un',0);
    expDates = [expDates{:}];
    all_used_bat_nums = [trialInfo_tmp.used_bat_nums];
    assert(issorted(expDates))
    
    switch binningMethod
        
        case 'equalCalls'
            n_date_bins = 10;
            dateIdx = round(linspace(1,length(expDates)+1,n_date_bins));
        case 'sessions'
            session_exp_dates = [trialInfo_tmp.expDate];
            assert(issorted(session_exp_dates))
            nCalls_by_session = zeros(1,length(session_exp_dates));
            for session_k = 1:length(session_exp_dates)
                nCalls_by_session(session_k) = sum(expDates == session_exp_dates(session_k));
            end
            usableSessions = nCalls_by_session > minCalls;
            usable_exp_dates = session_exp_dates(usableSessions);
            n_date_bins = length(usable_exp_dates);
            
            usableIdx = ismember(expDates,usable_exp_dates);
            
            expDates = expDates(usableIdx);
            X = X(:,:,usableIdx);
            all_used_bat_nums = all_used_bat_nums(usableIdx);
            
            dateIdx = nan(1,n_date_bins);
            for session_k = 1:n_date_bins
                dateIdx(session_k) = find(expDates == usable_exp_dates(session_k),1,'first');
            end
    end
    
    trialInfo_tmp = trialInfo_tmp(1);
    trialInfo_tmp.used_bat_nums = [];
    trialInfo_tmp.callID = [];
    
    %%
    
    FF = cell(1,n_date_bins-1);
    for k = 1:n_date_bins-1
        XX = X(:,:,dateIdx(k):dateIdx(k+1)-1);
        trialInfo_tmp.used_bat_nums = all_used_bat_nums(dateIdx(k):dateIdx(k+1)-1);
        [FF{k},t] = calculate_mvgc(XX,trialInfo_tmp,'timeWin',timeWin,'stepSize',stepSize,'modelOrder',modelOrder);
    end
    
    batNums = cellfun(@str2num,trialInfo_tmp.batNums);
    
    mvgc_over_time(exp_k) = struct('FF',{FF},'time',t,'expDates',expDates,'dateIdx',dateIdx,'all_used_bat_nums',all_used_bat_nums,'batNums',batNums,'expType','adult');
end