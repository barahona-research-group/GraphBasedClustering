%% test the graph constructions
nmis = [];
aris = [];
distancetype = 'euclidean';

%% Can change graph and parameters
networktype = 'rmst'; %graph 
weighted = 0;
p = [2,4,8]; %parameters
VI = [];
for k = 1:3
    para = p(k);
    D = squareform(pdist(data,distancetype));
    
    G = constructNetworkStructure(data',D,networktype,para);
    
    A = double(G);
    
    Time = 10.^[1:0.05:3];
    [Stb, N1, VI, C1] = stability(A,Time,'full','plot','p');
    
    nmit = N1;
    arit = N1;
    for i = 1:length(Stb)
        nmit(i) = nmi(C1(:,i),c_g);
        arit(i) = adjrand(C1(:,i)+1,c_g);
    end
    nmis(k) = max(nmit);
    aris(k) = max(arit);
    
end
disp(nmis);
disp(aris);

%% Plot the MDS figure
x = load('Data/synthetic.mat');
data = x.x;
n_g = 9;
D = squareform(pdist(data,'euclidean'));
G = constructNetworkStructure(data',D,'cknn',12);
A = double(G); %%8.6 9.9 15.9 23.8 26.6
Time = 10.^[0:0.02:3];
[Stb, N1, VI1, C1,c_all] = stability(A,Time,'full','plot','p','L',500,'M',500);
%[~,vmat] = varinfo(c_all');

n_Time = length(Time);
n_L = 500;

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
colormap(jet(n_L+125))
cmap = colormap;
cmap = cmap(51:550,:);
colormap(cmap)
idx0 = 0;
for i = 1:n_Time
    for k = 1:length(partitions{i})
        if vmat1(length(idx)+i,idx0+k)~=0
            plot(Time(i),-mds_x(idx0+k),'s',...
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
            plot(Time(i),-mds_x(idx0+k),'o','MarkerFaceColor',cmap(partitions_num{i}(k),:),...
            'Color',cmap(partitions_num{i}(k),:),'MarkerSize',5);
        end
    end
    idx0 = idx0 + length(partitions{i});
end
set(gca,'XScale','log')

yyaxis right
semilogx(Time,N1,'LineWidth',1),ylim([2,14])

mds_x1 = cmdscale(vmat,1);
mds_x = mds_x1(idx);

vis = reshape(mds_x,[100,151]);
N_max = max(c_all);
Ns = reshape(N_max,[100,151]);
figure,
hold on
colormap(jet(16))
cmap = colormap;
for k = 1:100   
        semilogx(Time,vis(k,:),'*','Color',cmap(Ns(k,:),:));
end