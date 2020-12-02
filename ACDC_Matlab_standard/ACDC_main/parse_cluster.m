function [molnames,nmols]=parse_cluster(namelabel)
% Input:
%    cluster name
% Output:
%    a cellstring and a vector containing the molecule names and
%    the corresponding molecule numbers, respectively

molnames=regexp(namelabel,'[A-Za-z]+','match');
nmols=regexp(namelabel,'[0-9]+','match');
nmols=cellfun(@str2num,nmols);

end