function [varargout]=get_rate_coefs(clus,ratetype,varargin)
% Get the collision or evaporation rates from the ACDC Matlab files
% Input:
%    clus = the name of the cluster for which the rate is wanted;
%        if clus is a cell string, it contains the colliding/evaporating parties,
%        if clus is a string, it is the evaporating cluster
%    ratetype = either 'coll' or 'evap'
% Optional input:
%    keyword 'ldisp' = print the results to the Matlab workspace
%    keyword 'lrun' = run the ACDC Perl script to create the rate files (see the settings below)
%    keyword 'dp' and the full directory path = directory for the rate files, if it's not the current one
%    keyword 'suffix' and the suffix string = suffix in the name of the rate files,
%        if  they are not simply 'get_coll.m' and 'get_evap.m'
% Optional output:
%    rate = the rate(s) in SI units (m^3 s^-1 for collisions, s^-1 for evaporations)
%    cluslabels = names of the colliding/evaporating parties


ldisp=0;
lrun_perl=0;
dp='';
suffix='';

if ~isempty(varargin)
    for i=1:length(varargin)
        if strcmpi(varargin{i},'ldisp')
            ldisp=1;
        elseif strcmpi(varargin{i},'lrun_perl') | strcmpi(varargin{i},'lrun')
            lrun_perl=1;
        elseif strcmpi(varargin{i},'dp')
            dp=varargin{i+1};
        elseif strcmpi(varargin{i},'suffix')
            suffix=varargin{i+1};
        end
    end
end

if lrun_perl
    disp('Running the Perl script: check that the Perl call is correct!')
    
    addpath('../Matlab_general','-end')
    
    % DeltaG and dip_pol files
    Gfile='./Energy_and_rate_files/HS298.15K_426clusters2016Apr25.txt';
    DPfile='./Energy_and_rate_files/dip_pol_298.15K_426clusters2016Apr25.txt';
    inputfile='./Cluster_set_files/input_AN_neutral.inp';
    rh=0;
    
    ACDCfile=find_newest_file('acdc.pl');
    commandstr=['perl ',ACDCfile,' --i ',inputfile,' --e ',Gfile,' --dip ',DPfile,' --rh ',num2str(rh)];
    lrun=system(commandstr);
    if lrun ~= 0
        error('Running Perl failed!')
    end
    rehash pathreset
end

% "Type" of the process
lcoll=0;
levap=0;
if strcmp(ratetype,'coll')
    lcoll=1;
    ratefile=[dp,'get_coll',suffix,'.m'];
elseif strcmp(ratetype,'evap')
    levap=1;
    ratefile=[dp,'get_evap',suffix,'.m'];
end

% Is this the rate for specific daughter clusters?
% If not, then it's all possible evaporations of the cluster
lspecific=0;
if levap && ischar(clus)
    cluslabels={};
    % The comments in the files should be in this format
    str0=['% ',clus,' -> '];
elseif iscellstr(clus) && length(clus)==2
    lspecific=1;
    cluslabels={[clus{1},' + ',clus{2}]};
    % The comments in the files should be in this format
    str1=[' ',clus{1},' + ',clus{2}];
    str2=[' ',clus{2},' + ',clus{1}];
end

if lspecific
    value=nan;
else
    j=0;
    cluslabels=cell(1,1);
    rate=nan(1,1);
end
    
fid=fopen(ratefile,'r');
line=fgetl(fid);
while ischar(line)
    if lspecific
        found1=strfind(line,str1);
        found2=strfind(line,str2);
        if (~isempty(found1) && (length(str1)==length(line(found1:end))...
                || isempty(regexp(line(found1+length(str1)),'[A-Za-z0-9]','match'))))...
                || (~isempty(found2) && (length(str2)==length(line(found2:end))...
                || isempty(regexp(line(found2+length(str2)),'[A-Za-z0-9]','match'))))
            %disp(line)
            split1=regexp(line,';','split');
            split2=regexp(split1{1},'=','split');
            split3=regexp(split1{end},'%','split');
            cluslabels_info=split3{end};
            eval(['value=',split2{end},';']);
            break
        end
    else
        if ~isempty(strfind(line,str0))
            %disp(line)
            split1=regexp(line,';','split');
            split2=regexp(split1{1},'=','split');
            split3=regexp(split1{end},'->','split');
            j=j+1;
            cluslabels{j}=split3{end};
            eval(['rate(j)=',split2{end},';']);
        end
    end
    line = fgetl(fid);
end
fclose(fid);

if lspecific
    if isnan(value) && ldisp
        disp(['Couldn''t find the rate from ',ratefile,' for ',clus])
    end
    rate=value;
else
    if isequaln(rate,nan(1,1)) && ldisp
        disp(['Couldn''t find any rates for ',clus])
    end
end

if lspecific && ldisp && exist('cluslabels_info','var')
    if lcoll
        disp([cluslabels_info,'    ',sprintf('%.2e',rate),' m^3 s^-1'])
    elseif levap
        disp([cluslabels_info,'    ',sprintf('%.2e',rate),' s^-1'])
    end
elseif ~lspecific
    [rate,ind]=sort(rate,'descend');
    cluslabels=cluslabels(ind);
    for j=1:length(cluslabels)
        if ldisp
            disp([cluslabels{j},'    ',sprintf('%.2e',rate(j)),' s^-1'])
        end
        if ~isempty(strfind(cluslabels{j},','))
            split1=regexp(cluslabels{j},',','split');
            cluslabels{j}=split1{1};
        end
    end
    if ldisp
        disp(['Total rate    ',sprintf('%.2e',sum(rate,'omitnan')),' s^-1'])
    end
end

if nargout>0
    output={rate,cluslabels};
    varargout=output(1:nargout);
end

end