function p = transformHMatrixtoPartitionVector(H)
%Transforms a given valid H matrix (Partition incidence
%matrix) into a partition vector 

[n,c] = size(H);

p = zeros(n,1);
    for i = 1:c
        p(H(:,i)==1) = i;
    end


end
