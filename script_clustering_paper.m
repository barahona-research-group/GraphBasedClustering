%% Iris
x = readtable('Data/bezdekIris.data.txt');
data = x{:,1:4};
c_g = ones(150,1);c_g(51:100) = 2;c_g(101:150) = 3;
n_g = 3;
%% Glass
x = readtable('Data/glass.data.txt');
data = x{:,2:10};
c_g = x{:,11};
c_g(c_g==7) = 4;
data = zscore(data);
n_g = 6;
%% Wine
x = readtable('Data/wine.data.txt');
data = x{:,2:14};
c_g = x{:,1};
data = zscore(data);
n_g = 3;
%% WDBC
x = readtable('Data/wdbc.data.txt');
data = x{:,3:32};
[~,~,c_g] = unique(x{:,2});
data = zscore(data);
n_g = 2;
%% Control chart
x = readtable('Data/synthetic_control.data.txt');
data = x{:,:};
c_g = repmat(1:6,100,1); c_g = c_g(:);
% data = zscore(data);
n_g = 6;
%% Parkinsons
x = readtable('Data/parkinsons.data.txt');
data = x{:,[2:17,19:24]};
c_g = x{:,18} + 1;
data = zscore(data);
n_g = 2;
%% vertebral
x = readtable('Data/column_3C.dat');
data = x{:,1:6};
[~,~,c_g] = unique(x{:,7});
% data = zscore(data);
n_g = 3;
%% breast tissue
x = readtable('Data/breast-tissue.txt');
data = x{:,3:11};
[~,~,c_g] = unique(x{:,2});
data = zscore(data);
n_g = 6;
%% seeds
x = readtable('Data/seeds_dataset.txt');
data = x{:,1:7};
c_g = x{:,8};
data = zscore(data);
n_g = 3;

%% segmentation
x = readtable('Data/segmentation.data.txt');
x1 = readtable('Data/segmentation.test.txt');
x = [x;x1];
data = x{:,2:end};
[~,~,c_g] = unique(x{:,1});
data = zscore(data);
n_g = 7;

%% yeast
x = readtable('Data/yeast.data.txt','Delimiter',' ','MultipleDelimsAsOne',1);
data = x{:,2:9};
[~,~,c_g] = unique(x{:,10});
data = zscore(data);
n_g = 10;

%% Synthetic data
x = load('Data/synthetic.mat');
data = x.x;

%% use the number of clusters in ground truth
k = 7;
D = squareform(pdist(data));
G = constructNetworkStructure(data', D,'cknn',k);
A = double(G);
Time = 10.^[1:0.02:3];
[Stb, N, VI, C] = stability(A,Time,'full','p','plot','v');

is = find(N==n_g);
if ~isempty(is)
    disp('Found partitions with N_ground_truth groups')
    C_tmp = C(:,is);
    nmi_tmp = is;
    ari_tmp = is;
    for i = 1:length(is)
        nmi_tmp(i) = nmi(C_tmp(:,i),c_g);
        ari_tmp(i) = adjrand(C_tmp(:,i),c_g);
    end
else
    k = min(abs(N-n_g));
    is = find(N==(n_g-k) | N==(n_g+k));
    C_tmp = C(:,is);
    for i = 1:length(is)
        nmi_tmp(i) = nmi(C_tmp(:,i),c_g);
        ari_tmp(i) = adjrand(C_tmp(:,i),c_g);
    end
end
nmis = max(nmi_tmp);
aris = max(ari_tmp);

%% analysis Markov Stability to choose the number of clusters
k = 7;
D = squareform(pdist(data));
G = constructNetworkStructure(data', D,'cknn',k);
A = double(G);
Time = 10.^[0:0.1:3];
[Stb, N, VI, C] = stability(A,Time,'full','plot');

% choose a scale from the plot, set the index to be i.
i = 31;
C_out = C(:,i);
nmi(C_out,c_g)
adjrand(C_out,c_g)
purity(C_out,c_g)
