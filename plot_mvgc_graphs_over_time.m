function plot_mvgc_graphs_over_time(mvgc_over_time,plotType)

edgeScale = 5;
nC = 255;
c = [0.75*ones(nC,1) repmat(linspace(0.9,0,nC)',1,2)];

nT = length(mvgc_over_time(1).FF);

switch plotType
    
    case 'graph'
        maxW = zeros(1,2);
        G = cell(2,nT);
        for exp_k = 1:2
            for t_k = 1:nT
                [~,idx] = min(abs(mvgc_over_time(exp_k).time));
                G{exp_k,t_k} = digraph(mvgc_over_time(exp_k).FF{t_k}(:,:,idx),'omitselfloops');
            end
            maxW(exp_k) = max(cellfun(@(x) max(x.Edges.Weight),G(exp_k,:)));
        end
        
        cla
        for exp_k = 1:2
            for t_k = 1:nT
                subplot(2,nT,nT*(exp_k-1) + t_k)
                edgeWeights = edgeScale *G{exp_k,t_k}.Edges.Weight/maxW(exp_k);
                plot(G{exp_k,t_k},'LineWidth',edgeWeights,'EdgeCData',edgeWeights,'Layout','circle','NodeColor','k','EdgeAlpha',1,'ArrowSize',12)
                axis square
                set(gca,'CLim',edgeScale*[0.1 1])
            end
        end
        
        colormap(c);
        
    case 'matrix'
        
        cla
        for exp_k = 1:2
            [~,idx] = min(abs(mvgc_over_time(exp_k).time));
%             maxW = max(cellfun(@(x) max(x(:,:,idx),[],'all'),mvgc_over_time(exp_k).FF));
            for t_k = 1:nT
                maxW = max(mvgc_over_time(exp_k).FF{t_k}(:,:,idx),[],'all');
                subplot(2,nT,nT*(exp_k-1) + t_k)
                plot_pw(edgeScale*mvgc_over_time(exp_k).FF{t_k}(:,:,idx)/maxW)
                axis square
                set(gca,'CLim',edgeScale*[0 1])
            end
        end
        
        colormap(c);
        
end