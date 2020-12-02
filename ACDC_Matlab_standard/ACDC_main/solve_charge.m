function [clus_ch,varargout] = solve_charge(clus,ch_info)
% Input:
%    cluster name, info on charged species' names
% Output:
%    charging state: 1=neutral, 2=neg., 3=pos.
% Optional output:
%    logical telling if clust is a cluster (1) or a generic ion (0)

clus_ch=nan;
clus_id=nan;

% First check if this is a generic ion
gen_ind=find(strcmp(clus,ch_info(end,:)));
if ~isempty(gen_ind)
    clus_id=0;
    clus_ch=gen_ind;
% If this is a molecular cluster, check for charged molecules
elseif check_cluster(clus)
    clus_id=1;
    clus_ch=1;
    [molnames,~]=parse_cluster(clus);
    for ch_ind=2:3
        mol_ch=ismember(molnames,ch_info(:,ch_ind));
        if ismember(1,mol_ch)
            clus_ch=ch_ind;
            break
        end
    end
end

if nargout > 1
    varargout={clus_id};
end

end