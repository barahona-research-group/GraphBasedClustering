function [S,N,C_new,VI] = stability_postprocess_data(S,N,C,VI,Time,A,params)
%stability_postprocess -- postprocessing of stability results to obtain
%smoothed stability curves and best possible partitions from the obtained
%results.
% Input:    S - path to matlab file in which stability results are
%                      stored eg 'matlab.mat' (created by save).
%           N - vector of the number of communities
%           C - partition matrix. Each column is a partition vector.
%           VI - vector of VI
%           Time - vector of Markov time
%           A - adjacency matrix of the original graph
%           params - struct storing the different options of stability
%requested variable names: S N C VI Time
%the effect of the function is to change the values in S, N, C_new and
%store into eg matlab_PP.mat (which can be recovered as
%load('matlab_PP.mat') )

% load stability files
% load(mat_file);

%TODO add other options for stability postprocessing

%% JEAN CHARLES -- you might want to set the default options differently.
%% for undirected graphs with normalized stability:
%% params.directed = false
%% params.linearise = false

if nargin<7
    params.directed = false;
    params.teleport_tau = 0.15;
    params.laplacian = 'normalized';
    params.linearised = false;
    params.precomputed = false;
end

% copy old results
C_new = C;
% then get set of all old found partitions
C = unique(C','rows');
C = C';
% make numbering start from 1
if min(C)==0
    C= C+1;
end
% get number of nodes
params.NbNodes = length(A);


%%%
% SWITCH CASES
if strcmpi(params.laplacian, 'normalized')
    if ~params.linearised
        StabilityFunction = @stability_FNL;
    else
        StabilityFunction = @stability_LNL;
    end
elseif strcmpi(params.laplacian, 'combinatorial')
    if ~params.linearised
        StabilityFunction = @stability_FCL;
    else
        error('linearised case not defined yet')
    end
else
    error('other cases not defined yet')
end



for j = 1:length(Time)
    fprintf('Post-processing Markov time: %f\t(%i of %i)\n', Time(j), ...
            j, length(Time))

    % get best stability from the found partitions
    [Stemp, Ntemp, Ctemp, VAROUT] = StabilityFunction(A, Time(j),C, params);
    if isfield(VAROUT,'precomputed')
        params.precomputed = VAROUT.precomputed;
        params.pi = VAROUT.pi;
        params.P = VAROUT.P;
    end

    % store it in respective fields..
    S(j) = Stemp;
    C_new(:,j) = Ctemp;
    N(j) = Ntemp;

    % store under old name, but append _PP.mat
%     save([mat_file(1:end-4) '_PP.mat'],'S','N','VI','C_new','Time');
end


end

%%%%%%%%%%%%%%%%%%%
% Different stability variants below.

function [S, N, Cbest, VAROUT] = stability_LNL(Graph, time, C, PARAMS)

%%%%
% FIRST PART: GENERATE "DYNAMICS" MATRIX

VAROUT =[]; % init varying outputs

% "transition matrix" and pi unknown
if PARAMS.precomputed == false
    % directed case: M_ij >> from i to j
    if PARAMS.directed == true
        error('not implemented yet')
        % undirected case
    else

        PI=sparse((diag(sum(Graph)))/sum(sum(Graph)));  %diag matrix with stat distr
        VAROUT.pi = diag(PI);   % store results for future use
        PARAMS.pi = diag(PI);

        solution=sparse(Graph/sum(Graph(:)));

    end

    % stationary distribution and "transition matrix" have been computed before
else
    if PARAMS.directed == true
        error('not implemented yet')
    else
        solution = Graph/sum(Graph(:));
    end
end

stabilities = zeros(1,size(C,2));
for i=1:size(C,2)
    H = transformPartitionVectorToHMatrix(C(:,i));
    stabilities(i) = (1-time) + trace(H' *(solution*time)*H) - (PARAMS.pi' *H)*(H'*PARAMS.pi);
end

index = find(stabilities==max(stabilities),1);

S = stabilities(index);
Cbest = C(:,index);
N = max(Cbest);

end



function [S, N, Cbest, VAROUT] = stability_FNL(Graph, time, C, PARAMS)

%%%%
% FIRST PART: GENERATE "DYNAMICS" MATRIX

VAROUT =[]; % init varying outputs

% "transition matrix" and pi unknown
if PARAMS.precomputed == false
    % directed case: M_ij >> from i to j
    if PARAMS.directed == true
        dout = sum(Graph,2);
        dangling = (dout==0);
        dout(dangling) = 1;
        Dout = sparse(diag(dout));
        clear dout;
        M = (1-PARAMS.teleport_tau)*(Dout\Graph); % deterministic part of transition
        % teleportation according to arXiv:0812.1770
        M =	M + diag(PARAMS.teleport_tau + dangling.*(1-PARAMS.teleport_tau))...
            * ones(PARAMS.NbNodes)/PARAMS.NbNodes;

        clear Dout dangling
        [v lambda_all] = eigs(M'); % largest eigenvalue of transition matrix corresponds to stat.distribution.
        lambda = max(diag(lambda_all));
        v = v(:,diag(lambda_all) == lambda);
        v = abs(v);              % make sure eigenvector is positive
        clear lambda;
        % store results for future use
        VAROUT.precomputed = true;
        VAROUT.pi = v/sum(v);
        PARAMS.pi = VAROUT.pi;
        VAROUT.P = M;

        % now compute exponential transition matrix
        solution = diag(v/sum(v))*expm(time* (M-eye(size(M))) );
        clear M v;
        % symmetrize solution
        solution = (solution+solution')/2;


        % undirected case
    else
        % Generate the matrix exponential
        trans=sparse(diag(    (sum(Graph)).^(-1)     ) * Graph);  %(stochastic) transition matrix
        Lap=sparse(trans-eye(PARAMS.NbNodes));
        % store results for future use
        VAROUT.precomputed = true;
        VAROUT.P = trans;

        clear trans;
        exponential=expm(time.*Lap);
        clear Lap;

        PI=sparse((diag(sum(Graph)))/sum(sum(Graph)));  %diag matrix with stat distr
        VAROUT.pi = diag(PI);   % store results for future use
        PARAMS.pi = diag(PI);

        solution=PI*exponential;
        clear exponential;
        clear PI;


    end

    % stationary distribution and "transition matrix" have been computed before
else
    if PARAMS.directed == true
        solution = diag(PARAMS.pi)*expm(time* (PARAMS.P - eye(size(PARAMS.P))) );
        solution = (solution +solution')/2; % symetrization needed for directed case
    else
        solution = diag(PARAMS.pi)*expm(time* (PARAMS.P - eye(size(PARAMS.P))) );
    end
end

stabilities = zeros(1,size(C,2));
for i=1:size(C,2)
    H = transformPartitionVectorToHMatrix(C(:,i));
    stabilities(i) = trace(H' *(solution - PARAMS.pi*PARAMS.pi'   ) *H);
end

index = find(stabilities==max(stabilities),1);

S = stabilities(index);
Cbest = C(:,index);
N = max(Cbest);

end


function [S, N, Cbest, VAROUT] = stability_FCL(Graph, time, C, PARAMS)
% Computes the full combinatorial stability

VAROUT =[]; % init varying outputs


% "transition matrix" and pi unknown
if PARAMS.precomputed == false
    % directed case: M_ij >> from i to j
    if PARAMS.directed == true
        % Directed part so far not implemented, different variants are possible
        error(['This parameter combination is not allowed! If you want to use directed stability, use the normalised Laplacian']);
        % undirected case
    else
        % standard Laplacian and average degree
        av_degree = sum(sum(Graph))/PARAMS.NbNodes;
        Lap=  -(Graph-diag(sum(Graph)));
        % store results for future use
        VAROUT.precomputed = true;
        VAROUT.P = Lap/av_degree;


        % Generate the matrix exponential
        exponential=expm(-time.*Lap/av_degree);
        clear Lap;

        PI=sparse(eye(PARAMS.NbNodes)/PARAMS.NbNodes);  %diag matrix with stat distr
        VAROUT.pi = diag(PI);   % store results for future use
        PARAMS.pi = diag(PI);

        solution=PI*exponential;
        clear exponential;
        clear PI;


    end

    % stationary distribution and "transition matrix" have been computed before
else
    if PARAMS.directed == true
        solution = diag(PARAMS.pi)*expm(-time* PARAMS.P);
        solution = (solution +solution')/2; % symetrization needed for directed case
    else
        solution = diag(PARAMS.pi)*expm(-time* PARAMS.P);
    end
end


stabilities = zeros(1,size(C,2));
for i=1:size(C,2)
    H = transformPartitionVectorToHMatrix(C(:,i));
    stabilities(i) = trace(H' *(solution - PARAMS.pi*PARAMS.pi'   ) *H);
end

index = find(stabilities==max(stabilities),1);

S = stabilities(index);
Cbest = C(:,index);
N = max(Cbest);

end

function H = transformPartitionVectorToHMatrix(pvector)
%Transforms a given partition vector into a valid H matrix (Partition incidence matrix)
% assumes that the ordering is contigous and starts from either zero or
% one.

if min(pvector)==0
    pvector = pvector+1;
end
nr_nodes = length(pvector);

% create H matrix
H = sparse(1:nr_nodes,pvector,1);

end