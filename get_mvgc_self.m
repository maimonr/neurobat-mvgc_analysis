function mvgcStruct_comb = get_mvgc_self(mvgcStruct,select_bat_type)

for exp_k = 1:length(mvgcStruct)
    nBat = length(mvgcStruct(exp_k).batNums);
    [FF_self,FF_self_conf] = deal(cell(1,nBat));
    
    for bat_k = 1:nBat
        ff_tmp = mvgcStruct(exp_k).FF{bat_k};
        ff_conf_tmp = mvgcStruct(exp_k).FF_conf{bat_k};
        switch select_bat_type
            case 'individual'
                FF_self{bat_k} = cat(1,ff_tmp(bat_k,:,:,:),permute(ff_tmp(:,bat_k,:,:),[2 1 3 4]));
                FF_self{bat_k} = FF_self{bat_k}(:,setdiff(1:nBat,bat_k),:,:);
                FF_self{bat_k} = squeeze(mean(FF_self{bat_k},2));
                
                FF_self_conf{bat_k} = cat(1,ff_conf_tmp(bat_k,:,:,:),permute(ff_conf_tmp(:,bat_k,:,:),[2 1 3 4]));
                FF_self_conf{bat_k} = FF_self_conf{bat_k}(:,setdiff(1:nBat,bat_k),:,:);
                FF_self_conf{bat_k} = squeeze(mean(FF_self_conf{bat_k},2));
            case 'nonfamiliar'
                FF_self{bat_k} = squeeze(nanmean(ff_tmp,[1 2]));
                if size(FF_self{bat_k}) == 1
                    FF_self{bat_k} = FF_self{bat_k}';
                end
                FF_self_conf{bat_k} = nanmean(ff_conf_tmp,[1 2]);
                FF_self_conf{bat_k} = reshape(FF_self_conf{bat_k},[1 size(FF_self_conf{bat_k},3) 2]);
        end
    end
    
    FF_self = cat(ndims(FF_self{1})+1,FF_self{:});
    FF_self_conf = cat(4,FF_self_conf{:});
    FF_self_conf = permute(FF_self_conf,[1 2 4 3]);
    
    [mvgcStruct(exp_k).FF_self,mvgcStruct(exp_k).FF_self_conf] = deal(FF_self,FF_self_conf);
end

%%
batNums = [mvgcStruct(1).batNums,mvgcStruct(2).batNums];
exp_bat_nums = [mvgcStruct(1).exp_bat_nums,mvgcStruct(2).exp_bat_nums];

[FF,FF_conf,FF_self,FF_self_conf] = deal(cell(1,length(mvgcStruct)));

if length(unique(cellfun(@length,{mvgcStruct.time}))) > 1
    [~,diff_t_idx(2)] = min(cellfun(@length,{mvgcStruct.time}));
    diff_t_idx(1) = setdiff(1:2,diff_t_idx(2));
    
    t_idx{diff_t_idx(1)} = true(1,length(mvgcStruct(diff_t_idx(1)).time));
    t_idx{diff_t_idx(2)} = false(1,length(mvgcStruct(diff_t_idx(2)).time));
    t = mvgcStruct(diff_t_idx(1)).time;
    [~,overlap_t_idx(1)] = min(abs(mvgcStruct(diff_t_idx(2)).time - mvgcStruct(diff_t_idx(1)).time(1)));
    [~,overlap_t_idx(2)] = min(abs(mvgcStruct(diff_t_idx(2)).time - mvgcStruct(diff_t_idx(1)).time(end)));
    t_idx{2}(overlap_t_idx(1):overlap_t_idx(2)) = true;
    
    for exp_k = 1:length(mvgcStruct)
        
        FF{exp_k} = mvgcStruct(exp_k).FF;
        FF{exp_k} = cellfun(@(ff) ff(:,:,t_idx{exp_k},:),FF{exp_k},'un',0);
        
        FF_conf{exp_k} = mvgcStruct(exp_k).FF_conf;
        FF_conf{exp_k} = cellfun(@(ff) ff(:,:,t_idx{exp_k},:),FF_conf{exp_k},'un',0);
        
        FF_self{exp_k} = mvgcStruct(exp_k).FF_self(:,t_idx{exp_k},:,:);
        FF_self_conf{exp_k} = mvgcStruct(exp_k).FF_self_conf(:,t_idx{exp_k},:,:);
        
    end
    
else
    for exp_k = 1:length(mvgcStruct)
        FF{exp_k} = mvgcStruct(exp_k).FF;
        FF_conf{exp_k} = mvgcStruct(exp_k).FF_conf;
        FF_self{exp_k} = mvgcStruct(exp_k).FF_self;
        FF_self_conf{exp_k} = mvgcStruct(exp_k).FF_self_conf;
    end
    t = mvgcStruct(1).time;
end

FF = [FF{:}];
FF_conf = [FF_conf{:}];

FF_self = cat(ndims(FF_self{1}),FF_self{:});
FF_self_conf = cat(3,FF_self_conf{:});

mvgcStruct_comb = struct('FF',{FF},'FF_conf',{FF_conf},'FF_self',FF_self,'FF_self_conf',FF_self_conf,'time',t,'batNums',batNums,'exp_bat_nums',exp_bat_nums);
end
