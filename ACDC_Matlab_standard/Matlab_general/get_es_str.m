function [es_str] = get_es_str(num,prec,varargin)
% Function for creating a neat output string for scientific notation
% Example: get_es_str(3.22571e4,1,[-1 1]) gives '3.2\cdot10^{4}'
%          get_es_str(2.3e-1,1,[-1 1])    gives '0.23'
% Input:
%    num     number for which the string is created
%    prec    wanted precision in terms of number of decimals in the coefficient part
% Optional input:
%    range in which fixed-point notation is used instead of exponential representation
%    (given as vector of the min and max values)
% Output:
%    the number as a text string

    if ~isempty(varargin)
        min_max=varargin{1};
        if num >= min_max(1) && num <= min_max(2)
            es_str=num2str(num);
            return
        end
    end

    es_exp=floor(log10(num));
    es_pref=num/10^es_exp;
    if strcmp(sprintf(['%0.',num2str(prec),'f'],es_pref),'10')
        es_str=['10^{',num2str(es_exp+1),'}'];
    elseif strcmp(sprintf(['%0.',num2str(prec),'f'],es_pref),'1')
        es_str=['10^{',num2str(es_exp),'}'];
    else
        es_str=[sprintf(['%0.',num2str(prec),'f'],es_pref),'\cdot','10^{',num2str(es_exp),'}'];
    end
    
end