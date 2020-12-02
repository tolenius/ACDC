function [lcluster]=check_cluster(namelabel)
% Input:
%    cluster name
% Output:
%    a logical telling if the input is a cluster, i.e. of the form
%    [numbers][letters][numbers][letters]... (1=yes, 0=no)

lcluster=0;

matchstr=regexp(namelabel,'([0-9]+[A-Za-z]+)+','match');
if strcmp(matchstr,namelabel)
    lcluster=1;
end

end