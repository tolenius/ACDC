function [fluxes,main_route,clust_bound] = track_fluxes(clust_flux,clust,ch_info,flux,outflux,varargin)
% Finds the main flux chains through the system; fluxes are also visualized for 2- and 3-component cases
% Input:
%    clust_flux:    cluster (or other) names corresponding to the fluxes
%    clust:         only the actual simulated clusters, i.e. no wall terms etc. or clusters outside the system boundary
%    ch_info:       info on charged species' names
%    flux:          matrix of all fluxes from ACDC
%    outflux:       matrix for the cluster pairs that collide to form the outgoing fluxes
% Optional input:
%    options for the charging state, the criterion for if a flux is taken into account individually or
%    put under the label 'others', ... (see below)
% Output:
%    fluxes:        a structure containing the fluxes (cm^-3 s^-1), and the names of the starting and ending clusters
%    main_route:    cluster names along the main growth route out of the simulation system
%    clust_bound:   clusters from which the outgoing fluxes originate

addpath(genpath('../Matlab_general'),'-end')

% Names of the (electrically neutral) vapors
name_mon=ch_info(1:end-2,1)';

%%%%%%%%% Default settings %%%%%%%%%

% Default criterion for plotting the flux to a certain cluster or out
% (fraction of the total flux into the cluster/out)
crit_out=0.01;
crit_clust=0.01;

% Charge: 1=neutral, 2=neg., 3=pos.
ch=1;

% Growth out of the system along specific channels (clusters or compositions),
% given as a cell array of cluster or (electrically neutral) molecule names
out_channel={};

% Plot the fluxes and write them out in the workspace
ldisp=1;

% See what is given as the function input
if ~isempty(varargin)
    for i=1:length(varargin)
        if ischar(varargin{i})
            if strcmpi(varargin{i},'crit_out')
                crit_out=varargin{i+1};
            elseif strcmpi(varargin{i},'crit_clust')
                crit_clust=varargin{i+1};
            elseif strcmpi(varargin{i},'ch')
                ch=varargin{i+1};
            elseif strcmpi(varargin{i},'out_channel')
                out_channel=varargin{i+1};
            elseif strcmpi(varargin{i},'ldisp')
                ldisp=varargin{i+1};
            elseif strcmpi(varargin{i},'fn')
                fn=varargin{i+1};
            end
        end
    end
end

% Override the setting for plotting the fluxes if there are more than 3 vapors, or only one vapor
% (can't visualize >3D data; 1-component case could be done but hasn't been implemented)
if length(name_mon) == 2 || length(name_mon) == 3
    lplot=ldisp;
else
    lplot=0;
end

%%%%%%%%% End of the settings section %%%%%%%%%

fluxes=struct;

% First see what's coming out
if ~exist('fn','var')
    [~,~,roundoutname,roundout] = plotflux_out(clust,ch_info,outflux,'ch',ch,'out_channel',out_channel,'crit',crit_out,'ldisp',ldisp);
else
    [~,~,roundoutname,roundout] = plotflux_out(clust,ch_info,outflux,'ch',ch,'out_channel',out_channel,'crit',crit_out,'ldisp',ldisp,'fn',fn);
end

k=0;
for i=1:length(roundoutname)
    lvalid=0;
    if ~isempty(strfind(roundoutname{i},'+'))
        clusters=regexp(roundoutname{i},'\+','split');
        for j=1:2
            clusters_ch(j)=solve_charge(clusters{j},ch_info);
        end
        % Include this flux if one of the colliders is of the wanted charge,
        % unless neutral pathways are wanted and the neutral collider is the smaller one
        if ismember(ch,clusters_ch)
            if ch==1
                if all(clusters_ch==ch)
                    lvalid=1;
                    [validclust,~,~,~]=max_mols(clusters{1},clusters{2},ch_info);
                else
                    if calc_mols(clusters{clusters_ch==ch},ch_info) >= calc_mols(clusters{clusters_ch~=ch},ch_info)
                        lvalid=1;
                        validclust=clusters{clusters_ch==ch};
                    end
                end
            else
                lvalid=1;
                validclust=clusters{clusters_ch==ch};
            end
        end
        if lvalid==1
            combined=combine_clusters(clusters{1},clusters{2},ch_info);
            k=k+1;
            fluxes(k).start=validclust;
            fluxes(k).end=combined;
            fluxes(k).val=roundout(i);
        end
    end
end

if k==0
    disp(['No contribution from ch = ',num2str(ch),' to the flux out for crit_out = ',num2str(crit_out)])
    return
end

% Total rate out through the plotted clusters (not including 'Others')
tot_out_plot=sum([fluxes.val],'omitnan');
titlestr=['Outgoing flux ',sprintf('%0.3g',tot_out_plot),' cm^{-3} s^{-1}\nthrough the depicted ',chstate_str(ch),' clusters'];
if ~isempty(out_channel)
    titlestr=[titlestr,' along the given channel(s)'];
end
% Another (optional) criterion for plotting the flux: threshold (compared to
% the total flux out along the major routes along the wanted channels)
%thflux=0.1*tot_out_plot;
thflux=0; % Set this to zero if you don't want to use it

% Save the names of the boundary clusters
clust_bound={fluxes.start};

nflux_out=k;
tracked={};
track_clus={};
for i1=1:nflux_out
    % Solve the origin of this, if it's a cluster
    if check_cluster(fluxes(i1).start) && ~ismember(fluxes(i1).start,tracked)
        start_ch=solve_charge(fluxes(i1).start,ch_info);
        if start_ch==ch
            if ~exist('fn','var')
                [~,~,mainfluxnames,mainfluxes]=plotflux(fluxes(i1).start,'to',clust_flux,clust,ch_info,flux,...
                    'crit',crit_clust,'lsame_ch','ldisp',ldisp);
            else
                [~,~,mainfluxnames,mainfluxes]=plotflux(fluxes(i1).start,'to',clust_flux,clust,ch_info,flux,...
                    'crit',crit_clust,'lsame_ch','ldisp',ldisp,'fn',fn);
            end
            for i2=1:length(mainfluxes)
                if mainfluxes(i2)>=thflux && ~strcmpi(mainfluxnames{i2},'others')
                    k=k+1;
                    fluxes(k).end=fluxes(i1).start;
                    fluxes(k).start=mainfluxnames{i2};
                    fluxes(k).val=mainfluxes(i2);
                    if ~ismember(mainfluxnames{i2},track_clus) && ~ismember(mainfluxnames{i2},tracked)
                        track_clus=[track_clus mainfluxnames{i2}];
                    end
                end
            end
        end
    end
    tracked=[tracked fluxes(i1).start];
end

while ~isempty(track_clus)
    track_clus_new={};
    for i1=1:length(track_clus)
        if check_cluster(track_clus{i1}) && ~ismember(track_clus{i1},tracked)
            start_ch=solve_charge(track_clus{i1},ch_info);
            if start_ch==ch
                if ~exist('fn','var')
                    [~,~,mainfluxnames,mainfluxes]=plotflux(track_clus{i1},'to',clust_flux,clust,ch_info,flux,...
                        'crit',crit_clust,'lsame_ch','ldisp',ldisp);
                else
                    [~,~,mainfluxnames,mainfluxes]=plotflux(track_clus{i1},'to',clust_flux,clust,ch_info,flux,...
                        'crit',crit_clust,'lsame_ch','ldisp',ldisp,'fn',fn);
                end
                for i2=1:length(mainfluxes)
                    if mainfluxes(i2)>=thflux && ~strcmpi(mainfluxnames{i2},'others')
                        k=k+1;
                        fluxes(k).end=track_clus{i1};
                        fluxes(k).start=mainfluxnames{i2};
                        fluxes(k).val=mainfluxes(i2);
                        if ~ismember(mainfluxnames{i2},track_clus_new) && ~ismember(mainfluxnames{i2},tracked)
                            track_clus_new=[track_clus_new mainfluxnames{i2}];
                        end
                    end
                end
            end
        end
        tracked=[tracked track_clus{i1}];
    end
    track_clus=track_clus_new;
end

% Find the main growth route
main_route={};
end_tmp=fluxes(1).start;
while check_cluster(end_tmp)
    main_route=[end_tmp, main_route];
    ind_all=strcmp(end_tmp,{fluxes.end});
    start_all={fluxes(ind_all).start};
    [~,ind_max]=max([fluxes(ind_all).val],[],'omitnan');
    end_tmp=start_all{ind_max};
end

if lplot == 1
    
    % Figure settings
    
    % Number of colors in the arrow plot
    nlevels=5;

    % Font size
    fsize=12;
    
    % Start plotting
    
    figure(100)
    clf(100)
    set(gcf,'Position',[50 150 500 420]*1.5);
    set(gca,'Position',[0.1 0.12 0.8 0.8]);
    set(gca,'FontSize',fsize,'FontWeight','normal');
    set(gcf,'Color','white')
    box on
    hold on
    title(sprintf(titlestr),'FontWeight','normal')
    
    % Find the largest numbers of molecules of each type, in order to
    % mark the system boundaries in the figure
    nmols_max=zeros(1,length(name_mon));
    for nc=1:length(clust)
        if check_cluster(clust{nc})
            nmols_mon=countnumbers(clust{nc},ch_info);
            for nmon=1:length(name_mon)
                if nmols_mon(nmon)>nmols_max(nmon)
                    nmols_max(nmon)=nmols_mon(nmon);
                end
            end
        end
    end
    
    % Set preliminary axes limits to prevent warnings about "ARROW changed the axis limits"
    maxlim=max(nmols_max);
    
    xlabel(['Number of molecules ',name_mon{1}])
    ylabel(['Number of molecules ',name_mon{2}])

    xlim([0 2*maxlim]);
    ylim([0 2*maxlim]);

    if length(name_mon)==3
        %view(3)
        view(150,30)
        zlabel(['Number of molecules ',name_mon{3}])
        zlim([0 2*maxlim]);
    end

    % Set arrow colors and widths
    arrowfluxes=[fluxes.val];
    for i=1:length(arrowfluxes)
        if ~check_cluster(fluxes(i).start) || ~check_cluster(fluxes(i).end)
            arrowfluxes(i)=nan;
        end
    end
    step=(log10(max(arrowfluxes,[],'omitnan'))-log10(min(arrowfluxes,[],'omitnan')))/nlevels;
    limits=log10(min(arrowfluxes,[],'omitnan')):step:log10(max(arrowfluxes,[],'omitnan'));
    Lwidth=1:0.5:0.5*(nlevels+1);
    Alength=14;
    % Arrow colors according to the flux magnitude
    lflux_colors=1;
    if length(name_mon)==3
        Alength=10;
        % 3D plots are easier to read when the coloring is according to
        % the molecular content along the z axis (in this case, the magnitude
        % of fluxes is depicted only by the arrow linewidth)
        lflux_colors=0;
    end
    if lflux_colors==1
        colors=makeColorMap([0 0.45 0.74],[1 0 0],nlevels); % blue -> red
        cmap=makeColorMap([0 0.45 0.74],[1 0 0]);
        %colors=flipud(hot(nlevels));
    else
        colors=[0 0 0; 0 0.45 0.74; 0 0.50 0; ...
            0.87 0.49 0; 0.64 0.08 0.18; 1 0 0];
        % Clusters outside the simulation system are depicted by one color
        cmap=colors(1:nmols_max(3)+2,:);
    end
    
    % Also find largest clusters at arrow endpoints to set axes limits
    nmols_max_end=zeros(1,length(name_mon));
    
    for i=1:size(fluxes,2)

        nmols_mon=countnumbers(fluxes(i).end,ch_info);
        for nmon=1:length(name_mon)
            if nmols_mon(nmon)>nmols_max_end(nmon)
                nmols_max_end(nmon)=nmols_mon(nmon);
            end
        end
        
        % Coordinates of the tracked cluster for arrows and possible text boxes
        endcoord=nmols_mon+0.5;
        text_coord_str=strcat(arrayfun(@(x) num2str(x), endcoord,'UniformOutput',false),',');
        text_coord_str=strcat(text_coord_str{:});
        if check_cluster(fluxes(i).start)
            start_ch=solve_charge(fluxes(i).start,ch_info);
            end_ch=solve_charge(fluxes(i).end,ch_info);
            if start_ch~=ch && end_ch==ch
                % Just mark the source to the figure, if it's not of the same
                % charge (unless it's just a neutral monomer)
                if ~ismember(fluxes(i).start,strcat('1',name_mon))
                    if start_ch==1
                        chcolor='k';
                    elseif start_ch==2
                        chcolor='b';
                    else
                        chcolor='r';
                    end
                    eval(['text(',text_coord_str,'fluxes(i).start,''EdgeColor'',''',chcolor,...
                        ''',''Fontsize'',',num2str(fsize),',''Linewidth'',2);'])
                end
            elseif start_ch==ch && end_ch==ch
                % Draw an arrow
                nmols_mon=countnumbers(fluxes(i).start,ch_info);
                startcoord=nmols_mon+0.5;
                nlev=nan;
                for j=1:length(limits)-1
                    if log10(fluxes(i).val)>=limits(j) && log10(fluxes(i).val)<=limits(j+1)
                        nlev=j;
                        break
                    end
                end
                ncol=nlev;
                if lflux_colors==0
                    ncol=min(nmols_mon(3)+1,nmols_max(3)+2);
                end
                arrow(startcoord,endcoord,'Linewidth',Lwidth(nlev),'Edgecolor',colors(ncol,:),...
                    'Facecolor',colors(ncol,:),'Length',Alength)
            end
        else
            % Mark the source to the figure, if it's not a cluster
            if ~strcmpi(fluxes(i).start,'source')
                eval(['text(',text_coord_str,'fluxes(i).start,''EdgeColor'',[0.8 0.8 0.8],''Fontsize'',',...
                    num2str(fsize),',''Linewidth'',2);'])
            end
        end
        
    end

    % Set colormap and colorbar
    colormap(cmap);
    if lflux_colors==1
        for i=1:length(limits)
            cbtick{i}=get_es_str(10^limits(i),0);
        end
        cb=colorbar('YTick',0:1/nlevels:1,'YTickLabel',cbtick,'FontWeight','normal','FontSize',fsize);
        ylabel(cb,'Flux (cm^{-3}s^{-1})','FontWeight','normal','FontSize',fsize)
    else
        cbtick_val=0:nmols_max(3)+1;
        cbtick=arrayfun(@(x) num2str(x), cbtick_val,'UniformOutput',false);
        cbtick{end}=['>',cbtick{end-1}];
        cb=colorbar('YTick',(0:1/length(cbtick_val):1)+1/length(cbtick_val)/2,'YTickLabel',cbtick,...
            'FontWeight','normal','FontSize',fsize);
        ylabel(cb,['Number of molecules ',name_mon{3}],'FontWeight','normal','FontSize',fsize)
    end
    
    % Limit axes
    maxlim=max(nmols_max_end); % Same limit for all axes in order to get even squares
    xlim([0 maxlim+1]);
    ylim([0 maxlim+1]);
    set(gca,'XTick',(0:maxlim)+0.5,'XTickLabel',0:maxlim);
    set(gca,'YTick',(0:maxlim)+0.5,'YTickLabel',0:maxlim);
    if length(name_mon)==3
        zlim([0 maxlim+1]);
        set(gca,'ZTick',(0:maxlim)+0.5,'ZTickLabel',0:maxlim);
%         % The axes label positioning doesn't work too well; haven't found a very good function for it
%         %axesLabelsAlign3D
%         h=rotate3d;
%         %set(h,'ActionPostCallback',@axesLabelsAlign3D);
%         set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
%         set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')')
    end

    % Points for grid lines
    grid_coord=cell(1,3);
    grid_coord(:)={[0]};
    for nmon=1:length(name_mon)
        if nmon==1
            grid_coord{nmon}=get(gca,'XTick')+0.5;
        elseif nmon==2
            grid_coord{nmon}=get(gca,'YTick')+0.5;
        else
            grid_coord{nmon}=get(gca,'ZTick')+0.5;
        end
        grid_coord{nmon}=[min(grid_coord{nmon}-1),grid_coord{nmon}];
    end
    % Draw the grid
    for nmon=1:length(name_mon)
        grid_coord_ij=cell(1,3);
        grid_coord_ij{nmon}=[grid_coord{nmon}(1) grid_coord{nmon}(end)];
        nother=setdiff(1:3,nmon);
        for i=1:length(grid_coord{nother(1)})
            grid_coord_ij{nother(1)}=[grid_coord{nother(1)}(i) grid_coord{nother(1)}(i)];
            for j=1:length(grid_coord{nother(2)})
                % Grid only at the edges of the box (full 3D grid may be hard to read)
                if (nmon==1 && i>1 && j>1) || ...
                    (nmon==2 && i>1 && j>1) || ...
                    (nmon==3 && i>1 && j>1)
                    continue
                end
                grid_coord_ij{nother(2)}=[grid_coord{nother(2)}(j) grid_coord{nother(2)}(j)];
                % Grid lines along the dimension corresponding to nmon
                h=plot3(grid_coord_ij{1},grid_coord_ij{2},grid_coord_ij{3},...
                    'Color',[0.38 0.38 0.38],'LineWidth',1,'LineStyle','-');
                uistack(h,'bottom')
                if i==length(grid_coord{nother(1)}) || j==length(grid_coord{nother(2)})
                    % Boundaries of the simulated system
                    line_coord=grid_coord_ij;
                    line_coord{nmon}=[0,nmols_max(nmon)+1];
                    [~,n]=max([i,j]);
                    line_coord{nother(n)}=[nmols_max(nother(n))+1,nmols_max(nother(n))+1];
                    line(line_coord{1},line_coord{2},line_coord{3},'Color','k','LineWidth',2.5,'LineStyle','--');
                end
            end
        end
    end

    hold off

end

end

function [nmols_mon]=countnumbers(namelabel,ch_info)
% Find the numbers of molecules for the matrix plot

name_mon=ch_info(1:end-2,1);
nmols_mon=zeros(1,length(name_mon));

[molnames,nmols]=parse_cluster(namelabel);
for nmon=1:length(name_mon)
    if ismember(name_mon{nmon},molnames)
        n=strcmp(molnames,name_mon{nmon});
        nmols_mon(nmon)=nmols(n);
    end
    for nch=2:3
        name_ch=ch_info{nmon,nch};
        if ~strcmp(name_ch,'') && ~strcmp(name_ch,ch_info{end-1,nch})
            % Add the ion to the number of molecules, if this is a charged cluster
            if ismember(name_ch,molnames)
                n_ch=strcmp(molnames,name_ch);
                if nmols(n_ch)~=1
                    error(['Number of ions is ',num2str(nmols(n_ch)),'?'])
                end
                nmols_mon(nmon)=nmols_mon(nmon)+nmols(n_ch);
                break
            end
        end
    end
end

end