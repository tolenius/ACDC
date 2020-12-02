function [lchannel] = check_channel(clus,ch_info,molnames_channel,varargin)
% Input:
%    cluster name, info on charged species' names (can be an empty array if there is no charge info),
%    cell array containing the (electrically neutral) molecule names forming the wanted channel
% Optional input:
%    logical telling if the cluster must contain all molnames_channel;
%    default: 0, i.e. all molecules don't need to be there as long as there is nothing extra
% Output:
%    a logical telling if the cluster is on the wanted clustering channel,
%    i.e. does not contain other molecules than molnames_channel

l_all=0;
if ~isempty(varargin)
    l_all=varargin{1};
end

[molnames,~]=parse_cluster(clus);

% If no charge info is given, assume all the given molecule names
msg='';
if isempty(ch_info)
    st=dbstack;
    msg=['#### NOTE: ',st(1).name,': no charge info given, result may be wrong for ions ####'];
    molnames_avail=unique([molnames,molnames_channel]);
    ch_info=cell(length(molnames_avail)+2,3);
    ch_info(:,:)={char.empty};
    ch_info(1:length(molnames_avail),1)=molnames_avail;
end

% Rows in ch_info for molnames_channel
[~,~,rows]=intersect(molnames_channel,ch_info(:,1));

lchannel=1;

for nmol=1:length(molnames)
    if ~ismember(molnames{nmol},ch_info(rows,:))
        % Generic ions that are able to charge molnames_channel are considered valid
        ind_gen=strcmp(molnames{nmol},ch_info(end,:));
        if ismember(1,ind_gen)
            if any(ch_info(rows,ind_gen) ~= "")
                continue
            end
        end
        disp(msg)
        lchannel=0;
        return
    end
end

if l_all
    for nr=1:length(rows)
        if ~any(ismember(ch_info(nr,1:2),molnames))
            lchannel=0;
            return
        end
    end
end

end