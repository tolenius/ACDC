function [lsame]=compare_clusters(clus1,clus2)
% Input:
%    names of two clusters
% Output:
%    a logical telling if the clusters have the same composition;
%    the molecules don't need to be in the same order (1=yes, 0=no)

clusters={clus1 clus2};

nmols=cell(1,2);
molnames=cell(1,2);

for i=1:2
    [molnames{i},nmols{i}]=parse_cluster(clusters{i});
    % sort the molecules in alphabetical order
    [molnames{i},index]=sort(molnames{i});
    if ~isempty(nmols{i})
        nmols{i}=nmols{i}(index);
    end
end

if isequal(molnames{1},molnames{2}) && isequal(nmols{1},nmols{2})
    lsame=1;
else
    lsame=0;
end

end