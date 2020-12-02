function [alloutname,allout,roundoutname,roundout] = plotflux_out(clust,ch_info,outflux,varargin)
% Finds the fluxes leading out from a simulation system
% Input:
%    clust: names of the simulated clusters
%    ch_info: info on charged species' names
%    outflux: matrix for the cluster pairs that collide to form the outgoing fluxes
% Optional input:
%    options for output, the criterion for if a flux is taken into account individually or
%    put under the label 'others', ... (see below)
% Output:
%    alloutname: cell string containing the names of the cluster pairs contributing to the outgoing flux, e.g. '5A5N+1A'
%    allout: vector containing the absolute values for all the outgoing fluxes
%    roundoutname: cell string containing the names of the outgrowing cluster pairs for all the fluxes that form at least
%                  a certain fraction of the total flux, determined by the variable crit
%    roundout: vector containing the absolute values for the fluxes that form a certain fraction of the total flux


% Default values for the input variables

crit=0.01;          % Criteria for the most significant fluxes (fraction of the total flux)
ch=0;               % Charging state (0=all, 1=neutral, 2=negative, 3=positive)
                    % (at least) one of the colliders must have the charge when ch > 0
out_channel={};     % Consider only specific outgrowth channels, either in terms of flux out from
                    % specific boundary clusters (out_channel is a cell array of cluster names), or flux out from
                    % clusters with a given composition (out_channel is a cell array of molecule names)

lshowall=0;         % Print all the outgoing fluxes to the screen (instead of only the main fluxes)
lclusters_out=0;    % List also the clusters that got out

% See what is given as the function input

if ~isempty(varargin)
    for i=1:length(varargin)
        if ischar(varargin{i})
            if strcmpi(varargin{i},'ch')
                ch=varargin{i+1};
            elseif strcmpi(varargin{i},'out_channel')
                out_channel=varargin{i+1};
            elseif strcmpi(varargin{i},'crit')
                crit=varargin{i+1};
            elseif strcmpi(varargin{i},'lshowall')
                lshowall=1;
            elseif strcmpi(varargin{i},'lclusters_out')
                lclusters_out=1;
            elseif strcmpi(varargin{i},'fn')
                fn=varargin{i+1};
            end
        end
    end
end

lclust_channel=0;
lmol_channel=0;
if ~isempty(out_channel)
    if check_cluster(out_channel{1})
        lclust_channel=1;
        for nco=1:length(out_channel)
            lfound=0;
            for nc=1:length(clust)
                if compare_clusters(clust{nc},out_channel{nco})
                    out_channel{nco}=clust{nc};
                    lfound=1;
                    break
                end
            end
            if lfound == 0
                disp([out_channel{nco},' is not in the set of simulated clusters'])
                out_channel{nco}='';
            end
        end
        out_channel(out_channel=="")=[];
    else
        lmol_channel=1;
    end
end

% Read the fluxes from a file instead of using the input variables,
% if a filename is given (the old way of obtaining the fluxes)
if exist('fn','var')
    if ~isempty(outflux)
        error('If you give a filename for the fluxes, you can''t also give a flux matrix')
    end
    % Read the header containing the cluster names
    fid=fopen(fn,'r');
    s=textscan(fid,'%s',1,'delimiter','\n');
    fclose(fid);
    s=strtrim(s);
    s2=regexp(s{1,1},'\s+','split');
    % Ignore the percent mark in the beginning of the header
    if strcmp(s2{1,1}{1},'%')
        clust=s2{1,1}(2:end);
    else
        error([s2{1,1}{1},'?'])
    end
    outflux=load(fn);
    if size(outflux,1) ~= length(clust) || size(outflux,2) ~= length(clust)
        error(['Number of rows or columns in ',fn,' doesn''t match with the header of the data set!'])
    end
end

% Find the contributions to the outgoing fluxes

outflux_clean=removedoubles(outflux,clust,ch_info);
outfluxsum=sum(sum(outflux_clean));

allout=[];
alloutname={};
alloutclust={};
l=0;
for i=1:size(outflux_clean,1)
    for j=1:size(outflux_clean,2)
        if outflux_clean(i,j)~=0
            clust_tmp={clust{i},clust{j}};
            if ch > 0
                for nc=1:2
                    ch_tmp(nc)=solve_charge(clust_tmp{nc},ch_info);
                end
                if ~ismember(ch,ch_tmp)
                    %disp(['Not of correct charge: ',clust{i},'+',clust{j}])
                    continue
                elseif ch==1 && ~all(ch_tmp==ch)
                    if calc_mols(clust_tmp{ch_tmp~=ch},ch_info) > calc_mols(clust_tmp{ch_tmp==ch},ch_info)
                        %disp(['Neutral cluster too small: ',clust{i},'+',clust{j}])
                        continue
                    end
                end
            end
            if lclust_channel
                if ~any(ismember({clust{i},clust{j}},out_channel))
                    %disp(['Not on cluster channel: ',clust{i},'+',clust{j}])
                    continue
                end
            elseif lmol_channel
                if ~check_channel(clust{i},ch_info,out_channel) || ~check_channel(clust{j},ch_info,out_channel)
                    %disp(['Not on composition channel: ',clust{i},'+',clust{j}])
                    continue
                end
            end
            l=l+1;
            allout(l)=outflux_clean(i,j);
            alloutname{l}=[clust{i},'+',clust{j}];
            if lclusters_out
                alloutclust{l}=combine_clusters(clust{i},clust{j},ch_info);
            end
        end
    end
end

% Sort them in descending order
[allout,ind]=sort(allout,'descend');
alloutname=alloutname(ind);
if lclusters_out
    alloutclust=alloutclust(ind);
end

% To make the plot clearer, just take the largest fluxes and put the rest under the label 'Others'
if isnan(crit)
    [roundoutname,roundout] = get_significant(alloutname,allout);
else
    [roundoutname,roundout] = get_significant(alloutname,allout,'crit',crit);
end

str_total=['Total flux out ',sprintf('%0.3g',outfluxsum),' cm^{-3} s^{-1}\n'];
str_rel='(relative to the total flux out';
str_spec='';
if ch > 0 || ~isempty(out_channel)
    
    outfluxsum_spec=sum(allout);
    plotsum=outfluxsum_spec;
    
    str_total=[str_total,'Flux out'];
    str_tmp='';
    
    if ch > 0
        str_tmp=[' through ',chstate_str(ch),' charge'];
    end
    if ~isempty(out_channel)
        str_tmp=[str_tmp,' along the given channel(s)'];
    end
    str_total=[str_total,str_tmp];
    str_rel=[str_rel,str_tmp,')'];
    
    if lclust_channel
        str_spec=[str_tmp,' ',strrep(strjoin(out_channel),' ',', ')];
    elseif lmol_channel
        str_spec=[str_tmp,' ',strrep(strjoin(out_channel),' ',' + ')];
    end
    
    str_total=[str_total,' ',sprintf('%0.3g',outfluxsum_spec),' cm^{-3} s^{-1} (',...
        sprintf('%0.0f',outfluxsum_spec/outfluxsum*100),' %%)'];
    
else
    plotsum=outfluxsum;
    str_rel=[str_rel,')'];
end
fprintf(['\n',str_total,'\n\n'])

disp(['Main collisions leading out',str_spec,': ']);
disp(str_rel)
for i=1:length(roundout)
    disp([roundoutname{i},' ',num2str(round(roundout(i)/plotsum*100)),'% '])
end
if lshowall
    disp('All collisions out: ');
    for i=1:length(allout)
        disp([alloutname{i},' ',sprintf('%0.3g',allout(i)/plotsum*100),'% '])
    end
end

figure(11)
set(gca,'LineWidth',2,'FontWeight','bold','FontSize',12);
set(gcf,'Color','white')
pie(roundout/plotsum)
hold on
legend(roundoutname,'Location','BestOutside')
title(sprintf(str_total))
hold off

% List the clusters getting out?
if lclusters_out
    % First see if the same outgoing cluster forms from different collisions
    allout_forclust=allout;
    for i=1:length(alloutclust)
        if ~strcmpi(alloutclust{i},'')
            for j=(i+1):length(alloutclust)
                lsame=compare_clusters(alloutclust{i},alloutclust{j});
                if lsame
                    % add the value to the current element and reset the value of the other element to zero
                    allout_forclust(i)=allout_forclust(i)+allout_forclust(j);
                    allout_forclust(j)=0;
                    alloutclust{j}='';
                end
            end
        end
    end
    % Remove empty elements
    ind=find(allout_forclust==0);
    allout_forclust(ind)=[];
    alloutclust(ind)=[];
    
    % Sort them in descending order
    [allout_forclust,ind]=sort(allout_forclust,'descend');
    alloutclust=alloutclust(ind);
    
    [roundoutclust,roundout_forclust] = get_significant(alloutclust,allout_forclust,'crit',crit);

    disp('Clusters getting out of the system: ');
    for i=1:length(roundout_forclust)
        disp([roundoutclust{i},' ',num2str(round(roundout_forclust(i)/plotsum*100)),'% '])
    end
    
    figure(2)
    set(gca,'LineWidth',2,'FontWeight','bold','FontSize',12);
    set(gcf,'Color','white')
    pie(roundout_forclust/plotsum)
    hold on
    legend(roundoutclust)
    title(sprintf(str_total))
    hold off
end

end

function [outflux_clean] = removedoubles(outflux,clust,ch_info)
% OLDER FORMAT: element (j,i) should be the same as (i,j); set the one corresponding to the smaller cluster
% to zero; in this case also fluxes due to collisions of two similar clusters need to be divided by two
% NEWER FORMAT: just see that the flux is in the element corresponding to the larger cluster
% (and the element corresponding to the smaller cluster is zero); change this, if needed
outflux_clean=zeros(size(outflux));
lnew_format=0;
for i=1:length(clust)
    for j=i:length(clust)
        if outflux(i,j)~=0
            if i~=j
                if outflux(j,i)==0
                    lnew_format=1;
                elseif outflux(i,j)~=outflux(j,i)
                    error(['Fluxes don''t match: ',clust{i},' and ',clust{j},' fluxes out ',...
                        sprintf('%0.3g',outflux(i,j)),' and ',sprintf('%0.3g',outflux(j,i))])
                end
                
                largest=max_mols(clust{i},clust{j},ch_info);
                if strcmpi(largest,clust{i})
                    outflux_clean(i,j)=outflux(i,j);
                    outflux_clean(j,i)=0;
                else
                    outflux_clean(j,i)=outflux(i,j);
                    outflux_clean(i,j)=0;
                end
            else
                outflux_clean(i,j)=outflux(i,j);
            end
        end
    end
end
if lnew_format == 0
    disp('Using the older outflux_matrix format! i.e. dividing fluxes from similar clusters by two')
    for i=1:length(clust)
        % Division by 2 since the fluxes are always from the point of view of the smaller cluster
        outflux_clean(i,i)=outflux_clean(i,i)/2;
    end
end

end
