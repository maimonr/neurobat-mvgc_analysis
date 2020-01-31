function [lfp_power, trialInfo] = prepare_call_ps_for_mvgc(call_trig_ps,varargin)

pnames = {'averageType','f_bin','eData','includedBats'};
dflts  = {'all',[70 150],[],[]};
[averageType,f_bin,eData,includedBats] = internal.stats.parseArgs(pnames,dflts,varargin{:});

recordedBats = cellflat(call_trig_ps.expParams.batNums);

if ~isempty(includedBats)
    batIdx = ismember(recordedBats,includedBats);
else
    batIdx = true(1,length(recordedBats));
end

recordedBats = recordedBats(batIdx);

nBat = sum(batIdx);
nTrial = size(call_trig_ps.ps,2);
nChannel = size(call_trig_ps.ps,3);
nT = size(call_trig_ps.ps,4);

t = call_trig_ps.expParams.ps_time';
playbackFlag = false;
if ~isfield(call_trig_ps.expParams,'included_bat_nums')
    playbackFlag = true;
else
    batNums = cellfun(@(x) mode(cellfun(@str2num, setdiff(cellflat(x),'unidentified'))),call_trig_ps.expParams.included_bat_nums);
end

lfp_call_time_s = abs(diff(call_trig_ps.expParams.call_t_win));
call_t_length = call_trig_ps.expParams.fs*lfp_call_time_s;
max_n_artifact = call_t_length*call_trig_ps.expParams.max_artifact_frac;
artifact_trial_idx = call_trig_ps.n_call_artifact_times>max_n_artifact;
artifact_removed_ps = call_trig_ps.ps(batIdx,:,:,:,:);

for bat_k = 1:nBat
    
    for trial_k = 1:nTrial
        for ch_k = 1:nChannel
            if artifact_trial_idx(bat_k,trial_k,ch_k)
                artifact_removed_ps(bat_k,trial_k,ch_k,:,:) = NaN;
            end
        end
    end
end

lfp_power_all = get_f_bin_lfp_power(artifact_removed_ps,call_trig_ps.specParams.freqs,f_bin);

switch averageType
    case 'all'
        lfp_power_all_avg = zeros(nBat,nT,nTrial);
        for bat_k = 1:nBat
            lfp_power_all_avg(bat_k,:,:) = squeeze(nanmean(lfp_power_all(bat_k,:,:,:),3))';
        end
        lfp_power = lfp_power_all_avg;
    case 'tetrode'
        nTt = 4;
        nChannel_per_tt = 4;
        lfp_power_all_avg = zeros(nBat*nTt,nT,nTrial);
        all_channel_k = 1;
        for bat_k = 1:nBat
            usedChannels = eData.activeChannels{strcmp(eData.batNums,recordedBats{bat_k})};
            tetrodeChannels = reshape(0:(nTt*nChannel_per_tt)-1,nTt,[]);
            
            for tt_k = 1:nTt
                channelIdx = ismember(usedChannels,tetrodeChannels(:,tt_k));
                current_lfp_power = squeeze(nanmean(lfp_power_all(bat_k,:,channelIdx,:),3))';
                lfp_power_all_avg(all_channel_k,:,:) = current_lfp_power;
                all_channel_k = all_channel_k + 1;
            end
        end
        lfp_power = lfp_power_all_avg;
    case 'none'
        nChannel = size(call_trig_ps.ps,3);
        lfp_power = zeros(nBat*nChannel,nT,nTrial);
        all_channel_k = 1;
        for bat_k = 1:nBat
            for channel_k = 1:nChannel
                current_lfp_power = squeeze(lfp_power_all(bat_k,:,channel_k,:))';
                lfp_power(all_channel_k,:,:) = current_lfp_power;
                all_channel_k = all_channel_k + 1;
            end
        end
end

used_trials = squeeze(~any(isnan(lfp_power),[1 2]));
lfp_power = lfp_power(:,:,used_trials);

trialInfo.batNums = recordedBats;
trialInfo.time = t;

if ~playbackFlag
    trialInfo.used_bat_nums = batNums(used_trials);
    trialInfo.callID = call_trig_ps.expParams.used_call_IDs(used_trials);
end


end