function [fluxname,fluxval,roundfluxname,roundfluxval] = plotflux(tar,direction,clust_flux,clust,ch_info,flux,varargin)
% Finds the fluxes into or out from a given cluster
% Input:
%    tar: cluster name
%    direction: 'to' or 'from'
%    clust_flux: cluster (or other) names corresponding to the fluxes
%    clust: only the actual simulated clusters, i.e. no wall terms etc. or clusters outside the system boundary
%    ch_info: info on charged species' names
%    flux: matrix of all fluxes from ACDC
% Optional input:
%    options for plotting, the criterion for if a flux is taken into account individually or
%    put under the label 'others', ... (see below)
% Output:
%    all the clusters contributing to the fluxes and the corresponding fluxes, and the names and fluxes where
%    the fluxes that don't meet the criterion are put together under the label 'others'

% Names of the generic ions (they need some specific considerations; see below)
genions=ch_info(end,2:3);

% Logicals that tell the direction of the flux
lfrom=0;
lto=0;
if strcmpi(direction,'to')
    lto=1;
elseif strcmpi(direction,'from')
    lfrom=1;
    % Flip the flux matrix
    flux=flux.';
else
    error('Input argument ''direction'' not well defined; should be ''to'' or ''from''.')
end

% Possible other input arguments

% Default values
crit=0.01;  % criteria for the most significant fluxes (fraction of the total flux)

lsame_ch=0; % track the fluxes based on the charge and not on the absolute number of molecules
lpairs=0;   % return both parties of a collision/evaporation (only the larger one is returned by default)
lgrowth=0;  % when the direction is 'from', growth to sizes that are larger than
                % 2 times the number of molecules in the growing ("target") cluster, i.e. if A+B->C,
                % C_mols<=2*A_mols, are classified as self-coagulation
ldisp=0;    % print the main fluxes to the screen and draw a pie chart
lshowall=0; % print all the fluxes to the screen (instead of just the main fluxes)

if ~isempty(varargin)
    for i=1:length(varargin)
        if strcmpi(varargin{i},'crit')
            crit=varargin{i+1};
        elseif strcmpi(varargin{i},'lsame_ch')
            if length(varargin)>i && isnumeric(varargin{i+1})
                lsame_ch=varargin{i+1};
            else
                lsame_ch=1;
            end
            if lsame_ch
                tar_ch=solve_charge(tar,ch_info);
            end
        elseif strcmpi(varargin{i},'lpairs')
            if length(varargin)>i && isnumeric(varargin{i+1})
                lpairs=varargin{i+1};
            else
                lpairs=1;
            end
        elseif strcmpi(varargin{i},'lgrowth')
            if length(varargin)>i && isnumeric(varargin{i+1})
                lgrowth=varargin{i+1};
            else
                lgrowth=1;
            end
        elseif strcmpi(varargin{i},'ldisp')
            if length(varargin)>i && isnumeric(varargin{i+1})
                ldisp=varargin{i+1};
            else
                ldisp=1;
            end
        elseif strcmpi(varargin{i},'lshowall')
            if length(varargin)>i && isnumeric(varargin{i+1})
                lshowall=varargin{i+1};
            else
                lshowall=1;
            end
        elseif strcmpi(varargin{i},'fn')
            fn=varargin{i+1};
        end
    end
end

% Read the fluxes from a file instead of using the input variables,
% if a filename is given (the old way of obtaining the fluxes)
if exist('fn','var')
    if ~isempty(flux)
        error('If you give a filename for the fluxes, you can''t also give a flux matrix')
    end
    % Read the header containing the cluster names + other flux names
    fid=fopen(fn,'r');
    s=textscan(fid,'%s',2,'delimiter','\n');
    fclose(fid);
    s=strtrim(s);
    s2=regexp(s{1,1}{1},'\s+','split');
    % Ignore the percent mark in the beginning of the header
    if strcmp(s2{1},'%')
        clust_flux=s2(2:end);
    else
        error([s2{1},'?'])
    end
    % If there is another header line, assume it's the names of the explicitly simulated clusters
    s3=regexp(s{1,1}{2},'\s+','split');
    if strcmp(s3{1},'%')
        clust=s3(2:end);
        disp('Found info on the second line of the flux file, assuming it''s the names of the simulated clusters')
    end
    flux=load(fn);
    if size(flux,2) ~= length(clust_flux)
        error(['Number of columns in ',fn,' doesn''t match with the header of the data set!'])
    end
end

% Find the column of the target cluster
col=find(strcmp(clust_flux,tar));
if isempty(col)
    error(['Cluster not found: ',tar])
end

% Size of the target cluster
tar_mols=calc_mols(tar,ch_info);

% See if the target cluster is a boundary cluster outside the true system
lbound_tar=0;
if check_cluster(tar) && ~ismember(tar,clust)
    lbound_tar=1;
end

% ltaken is a logical that tells if the cluster in question has already been taken into account in the flux
% (because fluxes from collisions are otherwise counted twice)
ltaken=zeros(1,length(clust_flux));

if lgrowth
    if ~lfrom
        error('lgrowth can be set to 1 only if you are tracking the fluxes FROM a cluster!')
    end
    % If we're studying cluster growth, save the self-coagulation fluxes to variable self_coag_flux
    self_coag_flux=0.0;
end

% Index in the flux vector
k=0;

for i=1:size(flux,1)
    
    if flux(i,col)~=0 && ltaken(i)==0
        
        % First check if this is a flux from a generic ion; then it should not be taken into account
        % (since the corresponding recombination/ionization flux will be missing, i.e. this is not a net flux)
        lgenion=strcmpi(clust_flux{i},genions);
        if ismember(1,lgenion)
            ltaken(i)=1;
            %disp(['Discarded flux from ',clust_flux{i}])
            continue
        end
        
        % If we're tracking the flux from a boundary cluster outside the system, we are not interested in the evaporating monomers
        % (this has to be excluded here since combining the product clusters doesn't work because there may be more than one monomer)
        if lbound_tar && lfrom
            [~,nmols]=parse_cluster(clust_flux{i});
            if length(nmols)==1
                if nmols(1)==1
                    ltaken(i)=1;
                    disp([tar,': Discarded flux to ',clust_flux{i}])
                    continue
                end
            end
        end
        
        lfound=0;
        
        if check_cluster(clust_flux{i})
            clust_mols=calc_mols(clust_flux{i},ch_info);
            if clust_mols > tar_mols
                % Count this as self-coagulation if the result cluster is "too" large
                if lgrowth && clust_mols >= 2*tar_mols
                    %disp(['Self-coagulation flux to ',clust_flux{i},' from ',tar])
                    if clust_mols > 2*tar_mols
                        self_coag_flux=self_coag_flux+flux(i,col);
                    else
                        % If the cluster collides with another similar cluster, one of them
                        % is considered to "grow", and the other one goes to self-coagulation
                        self_coag_flux=self_coag_flux+flux(i,col)./2;
                        
                        k=k+1;
                        fluxname{k,1}=clust_flux{i};
                        fluxval(k)=flux(i,col)./2;
                        if lpairs
                            fluxname{k,2}='';
                        end
                    end
                else
                    k=k+1;
                    fluxname{k,1}=clust_flux{i};
                    fluxval(k)=flux(i,col);
                    if lpairs
                        fluxname{k,2}='';
                    end
                end
                % No need to try to find the other collider, it just makes the script slow (and can't find it anyway)
                continue
            end
        end
        
        % If the flux is from a collision in the case of flux in, or from an evaporation in the case of flux out
        % (or ionization), find the other party
        for j=(i+1):size(flux,1)
            if ltaken(j)==0
                % Combine the clusters
                combined=combine_clusters(clust_flux{i},clust_flux{j},ch_info);
                % Check if the result is the target cluster
                lsame=compare_clusters(combined,tar);
                if lsame
                    lfound=1;
                    % This cluster has now been taken into account
                    ltaken(j)=1;
                    
                    k=k+1;
                    % Name of the cluster
                    [largest,smallest]=max_mols(clust_flux{i},clust_flux{j},ch_info);
                    if lsame_ch==1
                        % Choose the one with the same charge as the target cluster
                        clust_ch=solve_charge(largest,ch_info);
                        if clust_ch==tar_ch
                            fluxname{k,1}=largest;
                            if lpairs
                                fluxname{k,2}=smallest;
                            end
                        else
                            clust_ch=solve_charge(smallest,ch_info);
                            if clust_ch==tar_ch
                                fluxname{k,1}=smallest;
                                if lpairs
                                    fluxname{k,2}=largest;
                                end
                            else
                                % No match, choose the bigger one
                                fluxname{k,1}=largest;
                                if lpairs
                                    fluxname{k,2}=smallest;
                                end
                            end
                        end
                    else
                        % Choose the bigger one
                        fluxname{k,1}=largest;
                        if lpairs
                            fluxname{k,2}=smallest;
                        end
                    end
					
                    % Numerical value of the flux
                    fluxval(k)=flux(i,col);
                    
                    break
                end
            end
        end
        if lfound==0
            % This flux is not from a collision (for flux in)/evaporation (for flux out)...
            k=k+1;
            fluxname{k,1}=clust_flux{i};
            
            % ... except if it's two identical clusters colliding;
            % then the flux needs to be divided by 2, since it's from the
            % point of view of the colliding parties
            combined=combine_clusters(clust_flux{i},clust_flux{i},ch_info);
            lsame=compare_clusters(combined,tar);
            if lsame
                fluxval(k)=flux(i,col)./2;
                if lpairs
                    fluxname{k,2}=clust_flux{i};
                end
                %disp(['Divided: ',clust_flux{i}])
            else
                fluxval(k)=flux(i,col);
                if lpairs
                    fluxname{k,2}='';
                end
            end
        end
    end

end

if ~exist('fluxval','var')
    % No fluxes found?
    disp(['No fluxes ',direction,' ',tar,'!'])
    fluxval=[];
    fluxname={};
    roundfluxval=[];
    roundfluxname={};
else
    if lgrowth && self_coag_flux > 0
        k=k+1;
        fluxname{k,1}='self_coag';
        fluxval(k)=self_coag_flux;
        if lpairs
            fluxname{k,2}='';
        end
    end
    
    % Sort the fluxes in descending order
    [fluxval,ind]=sort(fluxval,'descend');
    fluxname=fluxname(ind,:);

    % To make the plot clearer, just take the largest fluxes and put
    % the rest under the label 'others'
    if isnan(crit)
        [roundfluxname,roundfluxval] = get_significant(fluxname,fluxval);
    else
        [roundfluxname,roundfluxval] = get_significant(fluxname,fluxval,'crit',crit);
    end
    
    if ldisp
        
        % Print the result to the workspace
        str=['Main fluxes ',direction,' ',tar,': '];        
        for i=1:length(roundfluxval)
            str=[str,'  ',roundfluxname{i,1},' ',num2str(round(roundfluxval(i)/sum(fluxval)*100)),'% '];
        end
        str=[str,'  Total flux ',sprintf('%0.3g',sum(fluxval)),' cm^-3 s^-1'];
        disp(str)
        if lshowall
            disp(['All fluxes ',direction,' ',tar,': ']);        
            for i=1:length(fluxval)
                disp([fluxname{i,1},' ',sprintf('%0.3g',fluxval(i)/sum(fluxval)*100),'% '])
            end
        end

        % Do a pie plot
        figure(12)
        set(gca,'LineWidth',2,'FontWeight','bold','FontSize',12);
        set(gcf,'Color','white')
        pie(roundfluxval/sum(fluxval))
        hold on
        legend(roundfluxname{:,1},'Location','BestOutside')
        title(['Total flux ',direction,' ',tar,' ',sprintf('%0.3g',sum(fluxval)),' cm^{-3}s^{-1}'])
        drawnow
        hold off

    end

end

end