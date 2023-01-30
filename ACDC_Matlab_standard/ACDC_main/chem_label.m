function [label] = chem_label(clus,ch_info,varargin)
% Input:
%    cluster name, info on charged species' names
% Optional input:
%    logical telling if ions (e.g. B = negative A) are converted to
%    their neutral form marked with a plus/minus (e.g. A-)
% Output:
%    cluster name in a "chemical" format, e.g. 2A3N1P -> (A_2N_3)+

l_ion_conv=1;
if ~isempty(varargin)
    l_ion_conv=varargin{1};
end

[molnames,nmols]=parse_cluster(clus);

label='';
ch_str='';

mol_charged='';
mol_neutral='';

ch=solve_charge(clus,ch_info);
if ch > 1
    if ~isempty(ch_info{end-1,ch})
        mol_charged=ch_info{end-1,ch};
    elseif l_ion_conv
        [inters,ind]=intersect(ch_info(1:end-2,ch),molnames);
        mol_charged=char(inters);
        mol_neutral=ch_info{ind,1};
        
        if ~ismember(mol_neutral,molnames)
            molnames=[molnames,mol_neutral];
            nmols=[nmols,0];
        end
    end
    
    if ch == 2
        ch_str='-';
    elseif ch == 3
        ch_str='+';
    end
end

for nmol=1:length(molnames)
    
    nmols_tmp=nmols(nmol);
    
    if strcmp(molnames{nmol},mol_charged)
        nmols_tmp=nmols_tmp-1;
    elseif strcmp(molnames{nmol},mol_neutral)
        nmols_tmp=nmols_tmp+1;
    end
    
    if nmols_tmp > 0
        label=[label,molnames{nmol}];
        if nmols_tmp > 1
            label=[label,'_{',num2str(nmols_tmp),'}'];
        end
    end
    
end

if ch_str ~= ""
    if sum(nmols) > 1
        label=['(',label,')'];
    end
    label=[label,'^',ch_str];
end

end