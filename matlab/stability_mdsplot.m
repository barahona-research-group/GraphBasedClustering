function stability_mdsplot(c_all,C1,N1,Time,n_L)
% visualise the partitions by MDS with VI distance  
n_Time = length(Time);

% unique the partitions
partitions = cell(n_Time,1);
for i = 1:length(Time)
    %vis = vmat(1+n_L*(i-1):n_L+n_L*(i-1),1+n_L*(i-1):n_L+n_L*(i-1));
    [~,vis] = varinfo(c_all(:,1+n_L*(i-1):n_L+n_L*(i-1))');
    partitions{i} = cell(rank(vis),1);
    %partitions{i}{1} = [1];
    num_k = 0;
    p = zeros(n_L,1);
    for k = 1:n_L
        ks = find(vis(k,:) == 0);
        ks = setdiff(ks,k:n_L);
        if isempty(ks)
            num_k = num_k+1;
            partitions{i}{num_k} = [k];
            p(k) = num_k;
        else
            partitions{i}{p(ks(1))} = [partitions{i}{p(ks(1))},k];
        end
    end
end

partitions_num = cell(n_Time,1);
for i = 1:n_Time
    j = [];
    l = [];
    for k = 1:length(partitions{i})
        j = [j,partitions{i}{k}(1)];
        l = [l,length(partitions{i}{k})];
    end
    partitions{i} = j;
    partitions_num{i} = l;
end

idx = [];
for i = 1:n_Time
    idx = [idx,partitions{i}+n_L*(i-1)];
end

%mds_x = cmdscale(vmat(idx,idx),1);
[~,vmat] = varinfo(c_all(:,idx)');
mds_x = cmdscale(vmat,1);

[~,vmat1] = varinfo([c_all(:,idx),C1]');

figure,
hold on
colormap(jet(round(n_L*4/3)))
cmap = colormap;
cmap = cmap(round(n_L/6):round(n_L/6)+n_L,:);
%colormap(jet(round(n_L*2)))
%cmap = colormap;
%cmap = cmap(round(n_L):n_L+n_L,:);

colormap(cmap)
idx0 = 0;
for i = 1:n_Time
    for k = 1:length(partitions{i})
        if vmat1(length(idx)+i,idx0+k)~=0
            p_1 = plot(Time(i),-mds_x(idx0+k),'s',...
            'Color',cmap(partitions_num{i}(k),:),'MarkerSize',4);
        %else
         %   plot(Time(i),-mds_x(idx0+k),'^','MarkerFaceColor',cmap(partitions_num{i}(k),:),...
          %  'Color',cmap(partitions_num{i}(k),:),'MarkerSize',4);
        end
    end
    idx0 = idx0 + length(partitions{i});
end
idx0 = 0;
for i = 1:n_Time
    for k = 1:length(partitions{i})
        if vmat1(length(idx)+i,idx0+k)==0
            p_2 = plot(Time(i),-mds_x(idx0+k),'o','MarkerFaceColor',cmap(partitions_num{i}(k),:),...
            'Color',cmap(partitions_num{i}(k),:),'MarkerSize',5);
        end
    end
    idx0 = idx0 + length(partitions{i});
end
set(gca,'XScale','log')
xlabel('Markov Time')
ylabel('MDS_1')
yyaxis right
semilogx(Time,N1,'LineWidth',1,'Color',[0 0.7 0]),ylim([max(1,min(N1)-1),max(N1)+1])
ylabel('Number of clusters','Color',[0 0.7 0])
ax = gca;
ax.YAxis(2).Color = [0 0.7 0];
hold off
cb = colorbar;
cb.Location = 'north';
cb.Position = [0.66,0.72,0.2,0.02];
cb.Label.String = 'Appearing frequency';
legend([p_1 p_2],'Louvain partition','MS partition')

