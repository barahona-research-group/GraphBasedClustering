%% compare with other methods.

nmit = [];
arit = [];
prty = [];
% compare given the number of clusters.

%% kmeans
[c_kmeans,C] = kmeans(data,n_g,'Replicates',50);
%% mixture models
obj = fitgmdist(data,n_g,'Replicates',10,'RegularizationValue',0.0001,'Start','RandSample');

c_gmm = cluster(obj,data);
%% hierarchical

Z = linkage(data,'complete','euclidean');
c_hcluster = cluster(Z,'maxclust',n_g);
%% ncut and NJW spectral clustering
neighbor_num = 12;
D = squareform(pdist(data,'euclidean'));
[D_LS,A_LS,LS] = scale_dist(D,floor(neighbor_num/2));

ZERO_DIAG = ~eye(size(data,1));
A_LS = A_LS.*ZERO_DIAG;

% spectral clustering with weighted matrix
n_c = n_g;
[c_ncut,x] = ncutW(A_LS,n_c);
c_ncut = transformHMatrixtoPartitionVector(c_ncut);
[c,x] = gcut(A_LS,n_c);
c_sc = c_ncut;
for i = 1:length(c)
    c_sc(c{i}) = i;
end

% spectral clustering with Cknn graph
G = constructNetworkStructure(data',D,'cknn',7);
A = double(G);
[c_ncut1,x] = ncutW(A,n_c);
c_ncut1 = transformHMatrixtoPartitionVector(c_ncut1);
[c1,x] = gcut(A,n_c);
c_sc1 = c_ncut1;
for i = 1:length(c1)
    c_sc1(c1{i}) = i;
end

%%

nmit(1) = nmi(c_kmeans,c_g);
nmit(2) = nmi(c_gmm,c_g);
nmit(3) = nmi(c_hcluster,c_g);
nmit(4) = nmi(c_ncut,c_g);
nmit(5) = nmi(c_sc,c_g);

arit(1) = adjrand(c_kmeans,c_g);
arit(2) = adjrand(c_gmm,c_g);
arit(3) = adjrand(c_hcluster,c_g);
arit(4) = adjrand(c_ncut,c_g);
arit(5) = adjrand(c_sc,c_g);

prty(1) = purity(c_kmeans,c_g);
prty(2) = purity(c_gmm,c_g);
prty(3) = purity(c_hcluster,c_g);
prty(4) = purity(c_ncut,c_g);
prty(5) = purity(c_sc,c_g);
