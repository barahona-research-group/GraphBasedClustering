# GraphBasedClustering

Multiscale graph-based clustering via Markov Stability
================
Zijing Liu

Introduction
------------

This contains the MATLAB codes for the paper "Graph-based data clustering via multiscale community detection" by Zijing Liu and Mauricio Barahona, published in Appl Netw Sci 5, 3 (2020).

Starting from data points, described as feature vectors, the method produces different geometric graphs and applies multiscale community detection (Markov Stability) to find graph partitions at different levels of resolution, which correspond to clusterings into different numbers of clusters. 

The graph-based clustering via Markov Stability uses the code in https://wwwf.imperial.ac.uk/~mpbara/Partition_Stability/ , also deposited in https://github.com/michaelschaub/PartitionStability

For an illustration, have a look at the notebook [MarkovStabilityClustering.ipynb](https://github.com/barahona-research-group/GraphBasedClustering/blob/master/MarkovStabilityClustering.ipynb)
*  script_clustering_paper.m - the example file for running the framework, same as the notebook. 
*  test_graph_build.m - test different graph constructions.
*  othercompare.m - clustering using other methods.
*  matlab/ - Matlab codes for graph construction, Markov Stability and two spectral clustering methods.
*  Data/ - public and generated datasets.


Citation
--------
Liu, Z., Barahona, M. Graph-based data clustering via multiscale community detection. Appl Netw Sci 5, 3 (2020)
https://doi.org/10.1007/s41109-019-0248-7
