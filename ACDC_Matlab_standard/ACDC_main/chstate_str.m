function [ch_name] = chstate_str(ch,varargin)
% Input:
%    number of the charging state (0=all, 1=neutral, 2=negative, 3=positive)
% Optional input:
%    logical telling if the first letter is capitalized
% Output:
%    name of the charging state as a string

l_upper=0;
if ~isempty(varargin)
    l_upper=varargin{1};
end

if ch == 0
    ch_name='all';
elseif ch == 1
    ch_name='neutral';
elseif ch == 2
    ch_name='negative';
elseif ch == 3
    ch_name='positive';
else
    disp(['Cannot recognize charging state ',num2str(ch)])
    ch_name='';
end

if l_upper == 1
    ch_name(1)=upper(ch_name(1));
end

end