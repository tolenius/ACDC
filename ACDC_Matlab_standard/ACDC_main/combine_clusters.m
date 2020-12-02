function [combined]=combine_clusters(clus1,clus2,ch_info)
% Input:
%    cluster names, info on charged species' names
% Output:
%    name of the cluster resulting from combining the two clusters

% NB! The order of the molecules in the combined cluster might be different than in
% your original ACDC name list - if needed, check it by e.g. using the function compare_clusters


combined='undef';
combined_ok=0;

clusters={clus1 clus2};

% Charging states and a logical for being a cluster (1) or a generic ion (0)
clus_ch=zeros(1,2);
clus_id=zeros(1,2);

% Vectors containing the names of different molecules, and
% the corresponding numbers of molecules for both clusters
nmols=cell(1,2);
molnames=cell(1,2);

for nclus=1:2
    [clus_ch(nclus),clus_id(nclus)]=solve_charge(clusters{nclus},ch_info);
    if isnan(clus_ch(nclus)) || isnan(clus_id(nclus))
        return
    end
    if clus_id(nclus)==1
        [molnames{nclus},nmols{nclus}]=parse_cluster(clusters{nclus});
    end
end

% Check all possible collisions, ionizations and recombinations
nmols_temp=cell(1,2);
molnames_temp=cell(1,2);

if all(clus_ch==1) || (all(clus_id==1) && ~all(clus_ch>1))
    % This is two neutral clusters or a neutral + a charged cluster aggregating
    [molsum,sumnames]=add_mols(nmols,molnames);
    combined_ok=1;
else
    if clus_ch(1)==clus_ch(2) || all(clus_id==0)
        % This doesn't produce a cluster
        return
    elseif any(clus_id==0)
        nclus=find(clus_id==1);
        ngen=find(clus_id==0);
        nmols_temp{1}=nmols{nclus};
        molnames_temp{1}=molnames{nclus};
        if clus_ch(nclus)>1
            % This is a recombination by a generic ion
            % Find the ionic species
            for nmol=1:length(molnames{nclus})
                nmol_info=find(strcmp(ch_info(1:end-2,clus_ch(nclus)),molnames{nclus}{nmol}));
                if ~isempty(nmol_info)
                    % Remove the ion
                    nmols_temp{1}(nmol)=nmols_temp{1}(nmol)-1;
                    if ~strcmp(molnames{nclus}{nmol},ch_info{end-1,clus_ch(nclus)})
                        % Also add a neutral molecule if this is a separate ionic
                        % molecule (i.e. not the proton / missing proton)
                        nmols_temp{2}=[1];
                        molnames_temp{2}={ch_info(nmol_info,1)};
                    end
                    [molsum,sumnames]=add_mols(nmols_temp,molnames_temp);
                    combined_ok=1;
                    break
                end
            end
        else
            % This might be an ionization by a generic ion
            % Find out if there is a molecule that can be ionized
            for nmol=1:length(molnames{nclus})
                nmol_info=find(strcmp(ch_info(1:end-2,clus_ch(nclus)),molnames{nclus}{nmol}));
                chmol=ch_info{nmol_info,clus_ch(ngen)};
                if ~strcmp(chmol,'')
                    % Add the ion
                    nmols_temp{2}=[1];
                    molnames_temp{2}={chmol};
                    if ~strcmp(chmol,ch_info{end-1,clus_ch(ngen)})
                        % Also remove a neutral molecule if this is a separate ionic
                        % molecule (i.e. not the proton / missing proton)
                        nmols_temp{1}(nmol)=nmols_temp{1}(nmol)-1;
                    end
                    [molsum,sumnames]=add_mols(nmols_temp,molnames_temp);
                    combined_ok=1;
                    break
                end
            end
        end
    else
        % This is a recombination of clusters
        nmols_temp=nmols;
        molnames_temp=molnames;
        for nclus=1:2
            % Find the ionic species and remove them
            for nmol=1:length(molnames{nclus})
                nmol_info=find(strcmp(ch_info(1:end-2,clus_ch(nclus)),molnames{nclus}{nmol}));
                if ~isempty(nmol_info)
                    nmols_temp{nclus}(nmol)=nmols_temp{nclus}(nmol)-1;
                    if ~strcmp(molnames{nclus}{nmol},ch_info{end-1,clus_ch(nclus)})
                        % Also add a neutral molecule if this is a separate ionic
                        % molecule (i.e. not the proton / missing proton)
                        nmols_temp{length(nmols_temp)+1}=[1];
                        molnames_temp{length(molnames_temp)+1}={ch_info(nmol_info,1)};
                    end
                    break
                end
            end
        end
        [molsum,sumnames]=add_mols(nmols_temp,molnames_temp);
        combined_ok=1;
    end
end

if combined_ok==0
    return
else
    combined=create_label(molsum,sumnames);
end

end

function [label]=create_label(nums,names)
% Takes in 2 vectors containing the numbers and names of the molecules and creates a name label
% e.g. [2 1 1], {'A' 'B' 'N'} -> 2A1B1N

if length(nums)~=length(names)
    error(['The amount of molecule numbers is not equal to the amount of molecule names: ',nums,', ',names])
end

label='';

for i=1:length(nums)
    if nums(i)>0
        label=[label,num2str(nums(i)),names{i}];
    end
end

end

function [molsum,sumnames]=add_mols(nmols,molnames)
% Takes in 2 cells containg vectors for the numbers and names of the molecules of different
% molecular clusters (2 or more), and combines the clusters by adding the molecules together

if ~iscell(nmols) || ~iscell(molnames)
    error(['Input should be cells: ',nmols,', ',molnames])
end

for i=1:length(nmols)
    if length(nmols{i})~=length(molnames{i})
        error(['The amount of molecule numbers is not equal to the amount of molecule names: ',nmols{i},', ',molnames{i}])
    end
end

% Create the vector for the names e.g. from the names in the 1st cluster
sumnames=molnames{1};
% Check if the other clusters have more molecule types
for i=2:length(molnames)
    for j=1:length(molnames{i})
        % Check if this is already in the name vector, if not, add it
        lfound=find(strcmp(sumnames,molnames{i}{j}));
        if isempty(lfound)
            sumnames=[sumnames,molnames{i}{j}];
        end
    end
end

molsum=zeros(1,length(sumnames));

% Sum together all the molecules that are of the same type
for i=1:length(nmols)
    for j=1:length(nmols{i})
        if nmols{i}(j)>0
            nmol=find(strcmp(sumnames,molnames{i}{j}));
            molsum(nmol)=molsum(nmol)+nmols{i}(j);
        end
    end
end

end


