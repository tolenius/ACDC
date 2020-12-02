function [roundnames,roundvalues] = get_significant(names,values,varargin)
% Input:
%    vectors containing names and corresponding values of entries contributing to an arbitrary quantity
%    (e.g. cluster flux names and values)
% Optional input:
%    a criterion for including individual entries of the input vector to the output
%    (relative to the total sum of the values vector, see the output description below);
%    given as a word-and-number pair 'crit',value
% Output:
%    vectors containing the names and values for all entries that are more than the given limit value relative to
%    the sum of the input number vector; rest of the entries are summed up under the label 'Others'

% default value for the limit (1% of the sum)
crit=0.01;

if ~isempty(varargin)
    for i=1:length(varargin)
        if strcmpi(varargin{i},'crit')
            crit=varargin{i+1};
        end
    end
end

if size(names,1)==1
    names=names';
end

roundvalues=[];
roundnames=cell(1,size(names,2));
k=0;
for i=1:length(values)
    % take just the significant values
    if (values(i)/sum(values,'omitnan') >= crit)
        k=k+1;
        roundvalues(k)=values(i);
        roundnames(k,:)=names(i,:);
    end
end
value_temp=sum(values,'omitnan')-sum(roundvalues,'omitnan');
if value_temp~=0
    k=k+1;
    roundvalues(k)=value_temp;
    roundnames{k,1}='Others';
end

end