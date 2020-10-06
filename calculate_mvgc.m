function [FF, t, confint, pval, sig] = calculate_mvgc(lfp_power,trialInfo,varargin)

pnames = {'trialSelection','selectBat','winSize','stepSize','timeWin','modelOrder','averageType','subSample','nSub','acmaxlags'};
dflts  = {'all',[],1,0.25,[-Inf,Inf],'AIC','tetrode','none',10,100};
[trialSelection,selectBat,winSize,stepSize,timeWin,mOrder_flag,averageType,subsampleFlag,nSub,acmaxlags] = internal.stats.parseArgs(pnames,dflts,varargin{:});

if iscell(lfp_power)
    X = cat(3,lfp_power{:});
    trialInfo = [trialInfo{:}];
else
    X = lfp_power;
end
sig_t = trialInfo(1).time;
fs = round(1/mean(diff(sig_t)));
[sig_t,t_idx] = inRange(sig_t,timeWin);

switch trialSelection
    case 'all'
        select_trial_idx = true(1,size(X,3));
    case 'includeBat'
        all_used_bat_nums = [trialInfo.used_bat_nums];
        select_trial_idx = ismember(all_used_bat_nums,selectBat);
    case 'excludeBat'
        all_used_bat_nums = [trialInfo.used_bat_nums];
        select_trial_idx = ~ismember(all_used_bat_nums,selectBat);
end

X_full = X(:,:,select_trial_idx);
X = X(:,t_idx,select_trial_idx);

assert(all(cellfun(@(x) isequal(trialInfo(1).batNums,x),{trialInfo.batNums})));
nBat = length(trialInfo(1).batNums);

switch averageType
    
    case 'all'
        channelIdx = (1:nBat)';
        
    case 'tetrode'
        nTT_per_bat = 4;
        channelIdx = reshape((1:nBat*nTT_per_bat)',nTT_per_bat,nBat);
        
    case 'none'
        nChannel_per_bat = 16;
        channelIdx = reshape((1:nBat*nChannel_per_bat)',nChannel_per_bat,nBat);
end

%%
tstat     = 'chi2';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'Bonferroni';  % multiple hypothesis test correction (see routine 'significance')

ntrials   = size(X,3);     % number of trials
nobs      = size(X,2);   % number of observations per trial

regmode   = '';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = '';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

momax     = 20;     % maximum model order for model order estimation

wind = round(winSize*fs);
ev = round(stepSize*fs);

if ischar(mOrder_flag)
    
    switch mOrder_flag
        case 'AIC'
            [~,~,mOrder] = tsdata_to_infocrit(X_full,momax,icregmode,false);
        case 'BIC'
            [~,~,~,mOrder] = tsdata_to_infocrit(X_full,momax,icregmode,false);
    end
    
elseif isnumeric(mOrder_flag)
    
    mOrder = mOrder_flag;
    
end

%%

wnobs = mOrder+wind;   % number of observations in "vertical slice"

sliding_win_idx = slidingWin(nobs,wnobs,wnobs-ev);
nWin = size(sliding_win_idx,1);

% loop through evaluation points
[pval, sig] = deal(nan(nBat,nBat,nWin));
confint = nan(nBat,nBat,nWin,2);
if strcmp(subsampleFlag,'none')
    FF = nan(nBat,nBat,nWin);
else
    FF = nan(nBat,nBat,nWin,nSub);
end

t = nan(1,nWin);

for win_k = 1:nWin
    t(win_k) = mean(sig_t(sliding_win_idx(win_k,:)));
    Xwin = X(:,sliding_win_idx(win_k,:),:);
    
    if strcmp(subsampleFlag,'none')
        G = get_autocov(Xwin,mOrder,regmode,acmaxlags);
        FF_win = nan(nBat);
    else
        G = NaN;
        FF_win = nan(nBat,nBat,nSub);
    end
    
    for channel_k1 = 1:nBat
        idx1 = channelIdx(:,channel_k1);
        for channel_k2 = 1:nBat
            if channel_k1 == channel_k2
                continue
            end
            idx2 = channelIdx(:,channel_k2);
            
            switch subsampleFlag
                
                case 'none'
                    
                    FF_win(channel_k1,channel_k2,:) = get_FF(G,idx1,idx2);
                    
                case 'perm'
                    
                    FP = permtest_tsdata_to_mvgc(Xwin,idx1,idx2,mOrder,5,nSub);
                    
                case 'boot'
                    
                    FF_win(channel_k1,channel_k2,:) = bootstrap_tsdata_to_mvgc(Xwin,idx1,idx2,mOrder,nSub,acmaxlags);
                    
                case 'empirical'
                    
%                     FE = empirical_var_to_mvgc(A,SIG,nobs,ntrials,idx1,idx2,0,10);
                    
                case 'trialShuffle'
                    parfor boot_k = 1:nSub
                        Xboot = Xwin;
                        permIdx = randperm(ntrials);
                        Xboot(idx1,:,:) = Xboot(idx1,:,permIdx);
                        G = get_autocov(Xboot,mOrder,regmode,acmaxlags);
                        FF_win(channel_k1,channel_k2,boot_k) = get_FF(G,idx1,idx2);
                    end
                    
            end
            
        end
    end
    
    FF(:,:,win_k,:) = FF_win;
    fprintf('%d windows out of %d done, %s\n',win_k,nWin,datestr(datetime));
%     if strcmp(subsampleFlag,'none')
%         pval(:,:,win_k) = mvgc_pval(FF_win,mOrder,nobs,ntrials,1,1,nBat-2,tstat);
%         [confint(:,:,win_k,1),confint(:,:,win_k,2)] = mvgc_confint(alpha,FF_win,mOrder,nobs,ntrials,1,1,nBat-2,tstat);
%         sig(:,:,win_k)  = significance(pval(:,:,win_k),alpha,mhtc);
%     end
end

end

function FF = get_FF(G,idx1,idx2)

FF = autocov_to_mvgc(G,idx1,idx2);

if isbad(FF,false)
    fprintf(2,' *** skipping - GC calculation failed\n');
    return
end

end

function [G,A,SIG,info] = get_autocov(X,mOrder,regmode,acmaxlags,varargin)
if ~isempty(varargin)
    verbose = varargin{1};
else
    verbose = false;
end
[A,SIG] = tsdata_to_var(X,mOrder,regmode);
if isbad(A)
    fprintf(2,' *** skipping - VAR estimation failed\n');
    return
end

[G,info] = var_to_autocov(A,SIG,acmaxlags);
if verbose && info.error
    fprintf(2,' *** skipping - bad VAR (%s)\n',info.errmsg);
    return
end
if verbose && info.aclags < info.acminlags % warn if number of autocov lags is too small (not a show-stopper)
    fprintf(2,' *** WARNING: minimum %d lags required (decay factor = %e)\n',info.acminlags,realpow(info.rho,info.aclags));
end
end
