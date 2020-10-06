function [lfpPower, trialInfo] = load_ps_data_mvgc(baseDir,expType,callType,used_exp_dates)

excl_bat_nums = {'71360','11682'};

switch callType
    
    case 'call'
        fName_str = '*call_trig_ps_corr.mat';
    case 'operant'
        fName_str = {'*call_trig_operant_box_1_ps_corr.mat','*call_trig_operant_box_2_ps_corr.mat'};
    case 'playback'
        fName_str = {'*playback_ps_corr.mat'};
end


lfpDir = fullfile(baseDir, [expType '_recording'],'data_analysis_results\lfp_data_analysis');
call_trig_ps_fnames = dir(fullfile(lfpDir,fName_str));
eData = ephysData(expType);
batNums = setdiff(eData.batNums,excl_bat_nums);

if ~isempty(used_exp_dates)
   expDates = arrayfun(@(x) datetime(x.name(1:8),'InputFormat','yyyyMMdd'),call_trig_ps_fnames);
   dateIdx = ismember(expDates,used_exp_dates);
   call_trig_ps_fnames = call_trig_ps_fnames(dateIdx);
end

nFiles = length(call_trig_ps_fnames);
[lfpPower, trialInfo] = deal(cell(1,nFiles));

for f_k = 1:nFiles
    m = matfile(fullfile(call_trig_ps_fnames(f_k).folder,call_trig_ps_fnames(f_k).name));
    expParams = m.expParams;
    
    switch callType
        
        case 'call'
            usable = isfield(expParams,'ps_time') && all(ismember(batNums,cellflat(expParams.batNums))) && ismember('cross_brain_corr',fieldnames(m));
        case 'operant'
            usable = isfield(expParams,'ps_time') && all(ismember(cellflat(expParams.batNums),batNums)) && ismember('cross_brain_corr',fieldnames(m));
        case 'playback'
            usable = isfield(expParams,'ps_time') && all(ismember(cellflat(expParams.batNums),batNums)) && ismember('cross_brain_corr',fieldnames(m));
    end
    
    if usable
        call_trig_ps = load(fullfile(call_trig_ps_fnames(f_k).folder,call_trig_ps_fnames(f_k).name),'ps','n_call_artifact_times','specParams','expParams');
        [lfpPower{f_k},trialInfo{f_k}] = prepare_call_ps_for_mvgc(call_trig_ps,'eData',eData,'includedBats',batNums,'averageType','tetrode');
        trialInfo{f_k}.expDate = repmat(datetime(call_trig_ps_fnames(f_k).name(1:8),'InputFormat','yyyyMMdd'),1,size(lfpPower{f_k},3));
    else
        disp('skipping experiment')
    end
end
