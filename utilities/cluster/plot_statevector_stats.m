function G = plot_statevector_stats(k, F, TM, MDT, NT,group)
%%%%%%%%%%%%%%%
% k is the number of clusters
% F frequency of occurrance (subjects * k)
% TM transition matrix (subjects * k * k)
% MDT Mean dwell time (subjects * k)
% NT number of transitions (subjects * 1)
% group (subjects * 1 integers with group IDs)

%%%%%%%%%%%%%
G = figure('color', 'w', 'Position',  [ 200         614        1011         446]);
load coldhot_white_128
if ~exist('group','var') || isempty(group)
    grp = ''; N = 0;
else
    grp = unique(group);N = length(grp)-1;
end

if isvector(F) 
    SS = 1; % single subject 
else
    SS = 0; % group
end

subplot(2,2+N,1)
if SS
    plot(1:k, F, 'k');
else
   
    if isempty(grp)
        plot_with_ste(gca, 1:k, F, []);%, 'k', [.8 .8 .8]);
    else
        lclr = [0 0 0;0 0 1;1 0 0;0 1 0];
        aclr = [.8 .8 .8;0 1 1;1 0 1;1 1 0];
        
        for ii = 1:length(grp)
            
            plot_with_ste(gca, 1:k, F(group==grp(ii),:), [], lclr(ii,:), aclr(ii,:));hold on;
        end
    end
end
box off; set(gca, 'TickDir', 'out')
ylabel('Frequency')
xlabel('State (cluster index)')


subplot(2,2+N,3+N)
if SS
    plot(1:k, MDT, 'k');
else
    if isempty(grp)
        plot_with_ste(gca, 1:k, MDT, []);%, 'k', [.8 .8 .8]);
    else
        lclr = [ 0  0  0; 0 0 1; 1 0 0; 0 1 0];
        aclr = [.8 .8 .8; 0 1 1; 1 0 1; 1 1 0];
        
        for ii = 1:length(grp)
            plot_with_ste(gca, 1:k, MDT(group==grp(ii),:), [], lclr(ii,:), aclr(ii,:));hold on;
        end
    end
    
end
box off; set(gca, 'TickDir', 'out')
ylabel('Mean dwell time (windows)')
xlabel('State (cluster index)')

CM = colormap('hot');

if SS
    subplot(2,2+N,[2,4])
    imagesc(TM, [0 1]); 
    %colormap('gray')
    C = colorbar;
    set(get(C, 'YLabel'), 'String', 'Probability')
    axis square
    set(gca, 'XTick', 1:k, 'YTick', 1:k)
    xlabel('State at t')
    ylabel('State at t-1')

    title(sprintf('Number of transitions: %d', NT))

    


    
else
    if isempty(grp)
        subplot(2,2+N,[2,4])
        %imagesc(squeeze(mean(TM)), [0 1]);
        imagesc(-log10(squeeze(mean(TM))+1e-5));%ED
        colormap(CM(length(CM)/8 +1:end,:))
        C = colorbar;
        set(get(C, 'YLabel'), 'String', '-log_1_0(Probability)')
        %set(get(C, 'YLabel'), 'String', 'Probability')
        axis square
        set(gca, 'XTick', 1:k, 'YTick', 1:k)
        xlabel('State at t')
        ylabel('State at t-1')
        
    else
        for ii = 1:length(grp)
            mxv(ii) = min(min(log10(squeeze(mean(TM(group==grp(ii),:,:)))+1e-5)));
        end
        mxa = max(mxv);
        
        for ii = 1:length(grp)
            subplot(2,2+N,[ii+1,ii+N+3])
            %imagesc(squeeze(mean(TM(group==grp(ii),:,:))), [0 1]);
            %imagesc(log10(squeeze(mean(TM(group==grp(ii),:,:)))), [mxa 0]);%ED
            imagesc(log10(squeeze(mean(TM(group==grp(ii),:,:)))), [floor(mxa) 0])
            title(sprintf('Number of transitions: %0.1f +/- %0.1f', mean(NT(group==grp(ii),:)), std(NT(group==grp(ii),:))))
            
            %set(get(C, 'YLabel'), 'String', '-log_1_0(Probability)')
            colormap(CM);
            
            C = colorbar;
            ylm = get(C,'YLim');
            %set(C,'YLim', [floor(ylm(1)) 0]);
            %set(C,'YTick', floor(ylm(1)):1:0,'YTickLabel',10.^(floor(ylm(1)):1:0))
            set(C,'YLim', [floor(mxa) 0]);
            set(C,'YTick', floor(mxa):1:0,'YTickLabel',10.^(floor(mxa):1:0))
            
            set(get(C, 'YLabel'), 'String', 'log_1_0(Probability)')
            axis square
            set(gca, 'XTick', 1:k, 'YTick', 1:k)
            xlabel('State at t')
            ylabel('State at t-1')
        end
        
    end
    
end
%colormap('gray')
%colormap(hot);
%C = colorbar;
%set(get(C, 'YLabel'), 'String', '-log_1_0(Probability)')
%set(get(C, 'YLabel'), 'String', 'log_1_0(Probability)')
%axis square
%set(gca, 'XTick', 1:k, 'YTick', 1:k)
%xlabel('State at t')
%ylabel('State at t-1')
if SS
    title(sprintf('Number of transitions: %d', NT))
else
    title(sprintf('Number of transitions: %0.1f +/- %0.1f', mean(NT), std(NT)))
end


