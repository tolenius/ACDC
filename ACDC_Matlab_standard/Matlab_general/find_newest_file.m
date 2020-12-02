function [fn_newest] = find_newest_file(fn_nodate,varargin)
% Finds the newest version of a file based on the date given in the file suffix,
% or optionally based on the timestamp
%
% Input:
%    fn_nodate: filename without the date suffix, e.g. 'file.txt' so that the files to be compared are
%        file_suffix.txt, where 'suffix' is a date of the form file_YYYYMMDD.txt or file_YYYY_MM_DD.txt
%        (no other non-numeric characters than '_' are allowed; YYYY, MM, DD etc. must be in the same order for all filenames)
% Optional input:
%    Keyword 'timestamp': Compare the timestamps of the files instead of filenames
% Output:
%    Full name or path of the newest file, i.e. fn_nodate with the date suffix

fn_newest='';

lstamp=0;

if ~isempty(varargin)
    for i=1:length(varargin)
        if strcmpi(varargin{i},'timestamp')
            lstamp=1;
        else
            error(['Cannot recognize input variable ',varargin{i}])
        end
    end
end

[fp,fn,fext]=fileparts(fn_nodate);
fpre=replace(fn_nodate,fext,'');

files=dir([fpre,'*',fext]);
nf_tot=numel(files);

time_number=nan(1,nf_tot);

for nf=1:nf_tot
    suffix=files(nf).name;
    suffix=replace(suffix,{fn,fext,'_'},'');
    
    if isempty(regexp(suffix,'^\d{8}$','once'))
        %disp(['Suffix ''',suffix,''' excluded'])
        continue
    end
    
    if lstamp == 0
        time_number(nf)=str2double(suffix);
    else
        time_number(nf)=files(nf).datenum;
    end
end

if ~all(isnan(time_number))
    [~,ind]=max(time_number,[],'omitnan');
    fn_newest=files(ind).name;
    
    if ~strcmp(fp,pwd)
        fn_newest=fullfile(fp,fn_newest);
    end
else
    if isfile(fn_nodate)
        fn_newest=fn_nodate;
        disp(['Returning ',fn_nodate,' as no files with date suffices were found'])
    else
        error(['Cannot find any files called ',fp,'/',fn,'[_DATE]',fext])
    end
end

end