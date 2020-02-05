dFF_boot = 100*abs(diff(FF_over_time_boot,[],3)./FF_over_time_boot(:,:,1:end-1,:));
shuffle_conf = squeeze(quantile(mean(dFF_boot,[1 2]),[0.05 0.95],4));
shuffle_avg = squeeze(mean(dFF_boot,[1 2 4]));

dFF = 100*abs(diff(FF_over_time,[],3)./FF_over_time(:,:,1:end-1));

date_t = 0:size(FF_over_time,3)-2;
boundedline(linspace(1,100,length(date_t)),squeeze(mean(dFF,[1 2])),squeeze(std(dFF,[],[1 2]))/sqrt(size(FF_over_time,2)*size(FF_over_time,1)),'r','alpha')
boundedline(linspace(1,100,length(date_t)),shuffle_avg,abs(shuffle_conf - repmat(shuffle_avg,1,2)),'k','alpha')
