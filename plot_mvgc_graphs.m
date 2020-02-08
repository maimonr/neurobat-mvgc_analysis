function plot_mvgc_graphs(FF_all,plotType,t)

t_bins = -2:0.5:1;
t_bins(2,:) = t_bins(1,:) + 1;
nT = size(t_bins,2);

edgeScale = 1.1;
nC = 255;
c = [0.75*ones(nC,1) repmat(linspace(0.9,0,nC)',1,2)];

cla

switch plotType
    case 'graph'
        
        G = cell(2,nT);
        maxW = zeros(1,2);
        for exp_k = 1:2
            for t_k = 1:nT
                [~,idx] = inRange(t{exp_k},t_bins(:,t_k));
                G{exp_k,t_k} = digraph(mean(FF_all{exp_k}(:,:,idx),3),'omitselfloops');
            end
            maxW(exp_k) = max(cellfun(@(x) max(x.Edges.Weight),G(exp_k,:)));
        end
        
        for exp_k = 1:2
            for t_k = 1:nT
                subplot(2,nT,nT*(exp_k-1) + t_k)
                edgeWeights = edgeScale *G{exp_k,t_k}.Edges.Weight/maxW(exp_k);
                plot(G{exp_k,t_k},'LineWidth',edgeWeights,'EdgeCData',edgeWeights,'Layout','circle','NodeColor','k','EdgeAlpha',1,'ArrowSize',12)
                axis square
                set(gca,'CLim',edgeScale *[0.1 1])
                if exp_k == 1
                    title(sprintf('%1.1f s',mean(t_bins(:,t_k))));
                end
            end
        end
    case 'matrix'
        for exp_k = 1:2
            maxW = max(FF_all{exp_k},[],'all');
            for t_k = 1:nT
                subplot(2,nT,nT*(exp_k-1) + t_k)
                [~,idx] = inRange(t{exp_k},t_bins(:,t_k));
                plot_pw(edgeScale*mean(FF_all{exp_k}(:,:,idx),3)/maxW)
                axis square
                set(gca,'CLim',[0 1])
                if exp_k == 1
                    title(sprintf('%1.1f s',mean(t_bins(:,t_k))));
                end
                h = gca;
                h.XTick = [];
                h.YTick = [];
            end
        end
        
end

colormap(c);