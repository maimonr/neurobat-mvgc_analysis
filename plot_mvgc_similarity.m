function plot_mvgc_similarity(FF_sim,FF_sim_boot)
alpha = 0.05/size(FF_sim,2);
boundedline(linspace(0,100,size(FF_sim,2)),mean(FF_sim_boot,[1 3]),quantile(FF_sim_boot,1-alpha,[1 3])' - mean(FF_sim_boot,[1 3])','k','alpha')
boundedline(linspace(0,100,size(FF_sim,2)),mean(FF_sim),std(FF_sim)./sqrt(size(FF_sim,1)),'r-o','alpha')
xlabel('% of experimental timeline')
set(gca,'FontSize',25)
h = findobj(gca,'Type','Line');
legend(h,'Shuffled (corrected alpha = 0.05)','Actual (SEM)')
legend box off
legend(h,'Shuffled (corrected alpha = 0.05)','Actual')
ylabel('Average Similarity')