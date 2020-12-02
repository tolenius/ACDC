function [totmols]=calc_mols(namelabel,ch_info)
% Input:
%    cluster name, info on charged species' names
% Output:
%    total number of molecules in the cluster;
%    a proton/removed proton is NOT considered to be a separate molecule

excluded=ch_info(end-1,:);

[molnames,nmols]=parse_cluster(namelabel);
totmols=sum(nmols);

for i=1:length(excluded)
    ind=find(strcmp(molnames,excluded{i}));
    if ~isempty(ind)
        totmols=totmols-nmols(ind);
    end
end

end