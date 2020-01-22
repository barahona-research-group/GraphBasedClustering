function z = purity(x, y)
% Compute purity of two clustering x and y.
% Input:
%   x, y: two integer vector of the same length. y is the given class(ground truth).
% Ouput:
%   z: purity
assert(numel(x) == numel(y));
n = numel(x);

t = crosstab(x,y);

z = sum(max(t'))/n;

