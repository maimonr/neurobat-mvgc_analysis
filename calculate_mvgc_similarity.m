function [FF_sim, FF_sim_boot, FF_over_time, FF_over_time_boot] = calculate_mvgc_similarity(mvgc_over_time)

if length(unique(cellfun(@length,{mvgc_over_time.time}))) > 1
    [~,short_t_idx] = min(cellfun(@length,{mvgc_over_time.time}));
    t_idx = cell(1,length(mvgc_over_time));
    t_bounds = mvgc_over_time(short_t_idx).time([1 end]) + 0.2*[-1 1];
    for k = 1:length(mvgc_over_time)
        if k ~= short_t_idx
            [~,t_idx{k}] = inRange(mvgc_over_time(k).time,t_bounds);
        else
            t_idx{k} = true(1,length(mvgc_over_time(k).time));
        end
    end
else
    t_idx = {true(1,unique(cellfun(@length,{mvgc_over_time.time})))};
end
t = mvgc_over_time(1).time;
nBoot = 1e3;

ff_tmp = arrayfun(@(exp,idx) cellfun(@(x) x(:,:,idx{1}),exp.FF,'un',0),mvgc_over_time,t_idx,'un',0);
idx = cellfun(@(ff) ~eye(size(ff{1},1)),ff_tmp,'un',0);
nPairs = sum(cellfun(@(x) sum(x,'all'),idx));

nDate = length(mvgc_over_time(1).dateIdx)-1;

FF_sim = nan(length(t),nDate);
FF_sim_boot = nan(length(t),nDate,nBoot);

FF_over_time = nan(length(t),nPairs,nDate);
FF_over_time_boot = nan(length(t),nPairs,nDate,nBoot);

for t_k = 1:length(t)
    a = cellfun(@(exp) cellfun(@(ff) ff(:,:,t_k),exp,'un',0), ff_tmp,'un',0);
    a = cellfun(@(x,idx) cellfun(@(a) a(idx),x,'un',0),a,idx,'un',0);
    a = vertcat(a{:});
    a = cell2mat(a);
    FF_over_time(t_k,:,:) = a;
    
    R = corr(a);
    R(logical(eye(size(R,1)))) = NaN;
    FF_sim(t_k,:) = nanmean(R,1);
    
    
    for boot_k = 1:nBoot
        aPerm = nan(size(a));
        for date_k = 1:size(a,2)
           permIdx = randperm(size(a,1));
           aPerm(:,date_k) = a(permIdx,date_k);
        end
        FF_over_time_boot(t_k,:,:,boot_k) = aPerm;
        R = corr(aPerm);
        R(logical(eye(size(R,1)))) = NaN;
        FF_sim_boot(t_k,:,boot_k) = nanmean(R,1);
    end
    
end
