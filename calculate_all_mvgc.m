

mOrder = 8;
timeWin = [-4 4];
stepSize = 0.25;

if ~exist('lfpPower','var')
    [lfpPower, trialInfo] = load_ps_data_mvgc('call');
end

if  ~exist('mvgcStruct','var')
    [mvgcStruct, mvgcStruct_comb] = calculate_mvgc_by_caller(lfpPower,trialInfo,'timeWin',timeWin,'stepSize',stepSize,'modelOrder',mOrder,'selectBatType','individual');
end

if ~exist('mvgc_over_time','var')
    mvgc_over_time = calculate_mvgc_over_time(lfpPower,trialInfo);
end

if ~exist('mvgcStruct_nonfamiliar','var')
    [mvgcStruct_nonfamiliar, mvgcStruct_comb_nonfamiliar] = calculate_mvgc_by_caller(lfpPower,trialInfo,'selectBatType','nonfamiliar','timeWin',timeWin,'stepSize',stepSize,'modelOrder',mOrder);
end

if ~exist('FF_sim','var')
    [FF_sim, FF_sim_boot] = calculate_mvgc_similarity(mvgc_over_time);
end

if ~exist('FF_all','var')
    FF_all = cell(1,2);
    for exp_k = 1:2
        [FF_all{exp_k},t] = calculate_mvgc(lfpPower{exp_k},trialInfo{exp_k},'timeWin',[-4 4],'stepSize',0.25,'modelOrder',8);
    end
end

if ~exist('lfpPower_operant','var')
    [lfpPower_operant, trialInfo_operant] = load_ps_data_mvgc('operant');
end

if ~exist('FF_operant_call','var')
    FF_operant_call = cell(1,2);
    for pair_k = 1:2
        FF_operant_call{pair_k} = calculate_mvgc(lfpPower_operant{pair_k},trialInfo_operant{pair_k},'stepSize',stepSize,'timeWin',timeWin,'modelOrder',mOrder);
    end
end

[lfpPower_playback, trialInfo_playback] = load_ps_data_mvgc('playback');
[FF_playback,t] = calculate_mvgc(lfpPower_playback{1},trialInfo_playback{1},'stepSize',stepSize,'timeWin',timeWin,'modelOrder',mOrder);
