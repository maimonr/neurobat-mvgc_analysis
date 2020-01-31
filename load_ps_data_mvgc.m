function [lfpPower, trialInfo] = load_ps_data_mvgc(callType)

switch callType
    
    case 'call'        
        expStrs = {'adult','adult_operant'};
        fName_str = repmat({'*call_trig_ps_corr.mat'},1,2);
    case 'operant'
        expStrs = repmat({'adult_operant'},1,2);
        fName_str = {'*call_trig_operant_box_1_ps_corr.mat','*call_trig_operant_box_2_ps_corr.mat'};
    case 'playback'     
        expStrs = {'adult'};
        fName_str = {'*playback_ps_corr.mat'};
end

[lfpPower,trialInfo] = deal(cell(1,length(expStrs)));
for exp_k = 1:length(expStrs)
    
    lfpDir = ['E:\ephys\' expStrs{exp_k} '_recording\data_analysis_results\lfp_data_analysis'];
    call_trig_ps_fnames = dir(fullfile(lfpDir,fName_str{exp_k}));
    eData = ephysData(expStrs{exp_k});
    batNums = setdiff(eData.batNums,'71360');
    nFiles = length(call_trig_ps_fnames);
    [lfpPower{exp_k}, trialInfo{exp_k}] = deal(cell(1,nFiles));
    
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
            [lfpPower{exp_k}{f_k},trialInfo{exp_k}{f_k}] = prepare_call_ps_for_mvgc(call_trig_ps,'eData',eData,'includedBats',batNums,'averageType','tetrode');
            trialInfo{exp_k}{f_k}.expDate = repmat(datetime(call_trig_ps_fnames(f_k).name(1:8),'InputFormat','yyyyMMdd'),1,size(lfpPower{exp_k}{f_k},3));
        else
            disp('skipping experiment')
        end
    end
    
end