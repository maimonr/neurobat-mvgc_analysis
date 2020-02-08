function [mvgcStruct, mvgcStruct_comb, mvgcStruct_nonfamiliar, mvgcStruct_comb_nonfamiliar, FF_sim, FF_sim_boot, FF_playback] = calculate_all_mvgc(mvgcData,varargin)

pnames = {'calculationSelection','outDir'};
dflts  = {'all',[]};
[calculationSelection,outDir] = internal.stats.parseArgs(pnames,dflts,varargin{:});

[mvgcStruct, mvgcStruct_comb, mvgcStruct_nonfamiliar, mvgcStruct_comb_nonfamiliar, FF_sim, FF_sim_boot, FF_playback] = deal(NaN);

if ~iscell(calculationSelection) && strcmp(calculationSelection,'all')
    calculationSelection = {'comm','over_time','familiar','similarity','allBats','operant','playback'};
end

saveData = ~isempty(outDir);
if saveData
    dateStr = datestr(datetime,'mmddyyyy');
    fName = fullfile(outDir,['mvgcData_' dateStr '.mat']);
    if ~exist(fName,'file')
        runTime = datetime;
        save(fName,'runTime');
    end
end

mOrder = 8;
timeWin = [-3 3];
stepSize = 0.5;
acmaxlags = 100;
subSample = 'boot';
nSub = 1e3;

if ismember('comm',calculationSelection)
    mvgcStruct = ...
        calculate_mvgc_by_caller({mvgcData.comm.lfpPower},{mvgcData.comm.trialInfo}...
        ,'timeWin',timeWin,'stepSize',stepSize,'modelOrder',mOrder,...
        'selectBatType','individual','acmaxlags',acmaxlags,'subSample',subSample,...
        'nSub',nSub);
    if saveData
        save(fName,'-append','mvgcStruct')
    end
    
    mvgcStruct_comb = get_mvgc_self(mvgcStruct,'individual');
    
    if saveData
        save(fName,'-append','mvgcStruct_comb')
    end
    
    sprintf('Finished calculating Comm. mvgc\n')
end

if ismember('over_time',calculationSelection)
    mvgc_over_time = calculate_mvgc_over_time({mvgcData.comm.lfpPower},{mvgcData.comm.trialInfo}...
        ,'timeWin',timeWin,'stepSize',stepSize,'modelOrder',mOrder,...
        'acmaxlags',acmaxlags,'subSample',subSample,'nSub',nSub);
    sprintf('Finished calculating mvgc over time\n')
    if saveData
        save(fName,'-append','mvgc_over_time')
    end
end

if ismember('familiar',calculationSelection)
    mvgcStruct_nonfamiliar =...
        calculate_mvgc_by_caller({mvgcData.comm.lfpPower},{mvgcData.comm.trialInfo},...
        'selectBatType','nonfamiliar','timeWin',timeWin,'stepSize',stepSize,...
        'modelOrder',mOrder,'acmaxlags',acmaxlags,'subSample',subSample,'nSub',nSub);
    mvgcStruct_comb_nonfamiliar = get_mvgc_self(mvgcStruct_nonfamiliar,'nonfamiliar');
    
    sprintf('Finished calculating familiar mvgc\n')
    if saveData
        save(fName,'-append','mvgcStruct_nonfamiliar', 'mvgcStruct_comb_nonfamiliar')
    end
end

if ismember('similarity',calculationSelection)
    [FF_sim, FF_sim_boot] = calculate_mvgc_similarity(mvgc_over_time);
    sprintf('Finished calculating mvgc similarity\n')
    if saveData
        save(fName,'-append','FF_sim', 'FF_sim_boot')
    end
    
end

if ismember('allBats',calculationSelection)
    FF_all = cell(1,2);
    for exp_k = 1:2
        FF_all{exp_k} = calculate_mvgc(mvgcData.comm(exp_k).lfpPower,mvgcData.comm(exp_k).trialInfo,...
            'timeWin',timeWin,'stepSize',stepSize,'modelOrder',mOrder,...
            'acmaxlags',acmaxlags,'subSample',subSample,'nSub',nSub);
    end
    sprintf('Finished calculating mvgc all bats\n')
    if saveData
        save(fName,'-append','FF_all')
    end
end

if ismember('operant',calculationSelection)
    FF_operant_call = cell(1,2);
    for pair_k = 1:2
        FF_operant_call{pair_k} = ...
            calculate_mvgc(mvgcData.operant(pair_k).lfpPower,mvgcData.operant(pair_k).trialInfo,...
            'stepSize',stepSize,'timeWin',timeWin,'modelOrder',mOrder,...
            'acmaxlags',acmaxlags,'subSample',subSample,'nSub',nSub);
    end
    sprintf('Finished calculating operant mvgc\n')
    if saveData
        save(fName,'-append','FF_operant_call')
    end
end

if ismember('playback',calculationSelection)
    [FF_playback,~] = ...
        calculate_mvgc(mvgcData.playback.lfpPower,mvgcData.playback.trialInfo,...
        'stepSize',stepSize,'timeWin',timeWin,'modelOrder',mOrder,...
        'acmaxlags',acmaxlags,'subSample',subSample,'nSub',nSub);
    sprintf('Finished calculating playback mvgc\n')
    if saveData
        save(fName,'-append','FF_playback')
    end
end
end
