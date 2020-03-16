function plot_mvgc_similarity(FF_sim,FF_sim_boot)
histogram(FF_sim_boot(:),bins,'Normalization','probability','FaceColor','k','FaceAlpha',0.5)
hold on
histogram(FF_sim(:),bins,'Normalization','probability','FaceColor','r','FaceAlpha',0.5)
axis square
xlabel('Connectivity matrix similarity (correlation)')
ylabel('Probability')
set(gca,'FontSize',25)
box off
legend('Suffle','Actual')
legend box off