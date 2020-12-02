function rates_and_deltags_ABe(fn_input)
% Syntax:
% rates_and_deltags_ABe(fn_input);
%
% Description:
% Create 2-dimensional matrix plots of evaporation rate constants and DeltaGs
%
% The plots are used to characterize the size- and composition-dependent cluster stability,
% and to assess if a cluster system is sufficiently large for ACDC simulations
% The plot axes correspond to the numbers of molecules of 2 components; additional components
% can be included by giving fixed numbers of additional molecules in each cluster
%
% Input is set in the separate input file given as the first input argument


addpath('../Matlab_general','-end')
get_constants;


% Run the input file
[~,~,fext]=fileparts(fn_input);
if ~strcmp(fext,'.m')
    error('The input file must be a .m file')
end
run(fn_input);


%%%%%%%%% Read in the data and do some pre-processing %%%%%%%%%

% Command for running the Perl (will be complemented below)
ACDCfile=find_newest_file('acdc.pl');
commandstr=['perl ',ACDCfile,' --keep_useless_collisions --i ',inputfile,' --e ',Gfile];
if exist('DPfile','var')
    commandstr=[commandstr,' --dip ',DPfile];
end

% Read the free energy file
fid=fopen(Gfile);

% Read the pressure and temperature (1st and 2nd lines)
s=textscan(fid,'%f',2,'Delimiter','\n');
if ~isnumeric(s{1})
    error(['The 1st and 2nd lines in ',Gfile,' should be the pressure and the temperature!'])
else
    pres=s{1}(1);
    if ~exist('temp','var')
        temp=s{1}(2);
    end
end

commandstr=[commandstr,' --temperature ',num2str(temp)];

% Read in the DeltaGs
% Find the number of columns
testline=fgetl(fid);
while isempty(testline) || ismember(1,regexp(testline,'\s*#'))
    testline=fgetl(fid);
end
ndatacols=numel(regexp(testline,'\S+'));
frewind(fid)
if ndatacols == 3
    M=textscan(fid,'%s %f %f','HeaderLines',2,'CommentStyle','#');
    clust=M{1,1};
    deltah=M{1,2};
    deltas=1e-3.*M{1,3}; % cal/molK to kcal/molK
    deltag_ref=deltah-temp*deltas;
elseif ndatacols == 2
    % DeltaGs
    M=textscan(fid, '%s %f','HeaderLines',2,'CommentStyle','#');
    clust=M{1,1};
    deltag_ref=M{1,2};
end

fclose(fid);

% Add the possible RH
if exist('rh','var')
    if rh > 0
        fprintf('\n#### RH > 0 %%: Note that effective DeltaGs can only be calculated for dry clusters ####\n')
        commandstr=[commandstr,' --rh ',num2str(rh)];
    else
        clear rh
    end
end

% Run the Perl script to create the rate constant files
if lrun_perl
    fprintf(['\nGenerating the rate constant files:\n',regexprep(commandstr,'\s+--','\n--'),'\n\n'])
    lrun=system(commandstr);
    if lrun ~= 0
        error('Running the Perl script failed!')
    end
    rehash pathreset
else
    fprintf('\n#### NOTE: Using previously generated ACDC files ####\n')
end

% Get the cluster and ion names from the driver file
fid=fopen('driver_acdc.m','r');
line=0;
while line ~= -1
    line=[fgetl(fid),' '];
    if ~isempty(regexp(line,'^\s*clust\s*=\s*{', 'once'))
        line=strrep(line,'clust','clust_acdc');
        eval(line);
    elseif ~isempty(regexp(line,'^\s*labels_ch\s*=\s*{', 'once'))
        eval(line);
    end
    if exist('clust_acdc','var') && exist('labels_ch','var')
        break
    end
end
fclose(fid);

%%%%%%%%% Find possible ions, charge carriers and extra molecule types %%%%%%%%%

ch_carrier='';
name_ch='';
mon_ch='';
l_ch_nonmol_AB=0;
l_ch_nonmol_extra=0;

% See what the possible extra molecules are
if l_inc_extra
    
    [molnames_extra,nmols_extra]=parse_cluster(clust_extra);
    ch_extra=solve_charge(clust_extra,labels_ch);
    
    if ch_extra > 1
        
        if ch ==1
            % If the extra species is charged, change the ch variable accordingly
            ch=ch_extra;
        elseif ch_extra ~= ch
            error(['Cannot have two opposite charges in the clusters: ch = ',num2str(ch),', clust_extra = ',clust_extra])
        end
        
        if length(molnames_extra) == 1 && ismember(molnames_extra{1},labels_ch(end-1,:))
            % The extra species is actually only a charge (i.e. not an actual molecule)
            ch=find(strcmp(labels_ch(end-1,:),molnames_extra{1}));
            l_inc_extra=0;
            
            fprintf(['\nYou have given ',clust_extra,' as an additional species:\n'])
            disp(['    Assuming that you want to plot charged 2-component ',A,B,'-clusters'])
            disp('    - but you are not very good at communicating it')
        else
            % The extra species may be an ionic molecule, find it here - it may not have a neutral counterpart
            % It may also be charged A or B, but if it's given as an extra species, it's kept as is
            % in order to have the possibility to plot A and B also without including their ion forms
            for nmol=1:length(molnames_extra)
                if solve_charge(['1',molnames_extra{nmol}],labels_ch) > 1
                    ch_carrier=molnames_extra{nmol};
                    name_ch=molnames_extra{nmol};
                    mon_ch=['1',name_ch];
                    break
                end
            end
        end
        
    end
    
end

if ~l_inc_extra
    molnames_extra=nan;
    nmols_extra=nan;
    ch_extra=nan;
end

if ch > 1
    if ch_carrier == "" || ch_carrier == name_ch
        
        % Find the electrically neutral molecule that acts as a charge carrier
        if ch_extra > 1
            molnames_tmp=setdiff(molnames_extra,name_ch);
        else
            molnames_tmp={A, B};
        end
        
        for nmol=1:length(molnames_tmp)
            name_tmp=molnames_tmp{nmol};
            name_ch_tmp=labels_ch{strcmp(labels_ch(1:end-2,1),name_tmp),ch};
            
            if ~isempty(name_ch_tmp)
                ch_carrier=name_tmp;
                name_ch=name_ch_tmp;
                if name_ch == labels_ch{end-1,ch}
                    if name_tmp == A || name_tmp == B
                        l_ch_nonmol_AB=1;
                    else
                        l_ch_nonmol_extra=1;
                    end
                    mon_ch=['1',ch_carrier,'1',name_ch];
                else
                    mon_ch=['1',name_ch];
                end
                break
            end
        end
        
    end
    if name_ch == ""
        error(['Could not find an ion with charge no. ',num2str(ch)])
    end
    
    fprintf('\nConsidering charged clusters:\n')
    disp(['    Charge carrier: ',ch_carrier])
    disp(['    Charged species: ',name_ch])
    disp(['    Charged monomer: ',mon_ch])
end


% Fill in monomer DeltaGs

% Possible non-zero value for a charged monomer (in the case that the DeltaGs are calculated
% with respect to a specific ion, and other ions of the same polarity are scaled accordingly)
scale_ch=nan;

% Reference DeltaGs should be zero for monomers (they are probably not
% listed in the DeltaG file since they are assumed to be zero by default)
for nmol=1:2
    
    if nmol == 1
        name_tmp=A;
    else
        name_tmp=B;
    end
    if ch_carrier == name_tmp
        clustname=mon_ch;
    else
        clustname=['1',name_tmp];
    end
    
    if ismember(clustname,clust)
        % If a charged monomer is found in the DeltaG file with a non-zero DeltaG, set it to zero and
        % subtract the value from the DeltaGs of all charged clusters of the same polarity
        if clustname == mon_ch
            nclust=find(strcmp(clust,clustname));
            if deltag_ref(nclust) ~= 0
                scale_ch=deltag_ref(nclust);
            end
        end
    else
        % Add the monomer
        nclust=length(clust)+1;
        clust{nclust}=clustname;
        deltag_ref(nclust)=0;
    end
end


% Indices of neutral, negative and positive clusters
num=cell(3,1);
for nch=1:3
    num{nch}=[];
end

for nc=1:length(clust)
    clus_ch=solve_charge(clust{nc},labels_ch);
    num{clus_ch}=[num{clus_ch} nc];
end

% Scale the DeltaGs of charged clusters, if needed
if ch > 1 && ~isnan(scale_ch)
    for nc=1:length(num{ch})
        nclust=num{ch}(nc);
        deltag_ref(nclust)=deltag_ref(nclust)-scale_ch;
    end
end

% for nch=1:3
%     figure(nch)
%     hold on
%     plot(1:length(num{nch}),deltag_ref(num{nch}),'-o');
%     set(gca,'XTickLabel',clust(num{nch}),'XTick',1:length(num{nch}));
%     rotateticklabel(gca,45);
% end

% Maximum number of different molecule types in a cluster
if ch == 1
    % Types A and B
    ntypes_max=2;
else
    % Types A and B, and a charge carrier
    ntypes_max=3;
end
% Possible additional molecule type
if l_inc_extra
    ntypes_max=ntypes_max+length(nmols_extra);
    if ismember(name_ch,molnames_extra)
        ntypes_max=ntypes_max-1;
    end
end

% Default maximum numbers of A and B molecules in the clusters, if not given
if strcmpi(B,'N') || strcmpi(B,'W')
    maxA_def=5;
    maxB_def=5;
elseif strcmpi(B,'M')
    maxA_def=3;
    maxB_def=3;
elseif strcmpi(B,'D')
    maxA_def=4;
    maxB_def=4;
elseif strcmpi(B,'T')
    maxA_def=2;
    maxB_def=2;
else
    maxA_def=4;
    maxB_def=4;
end
if ~exist('maxA','var')
    maxA=maxA_def;
end
if ~exist('maxB','var')
    maxB=maxB_def;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Calculate the wanted quantities and create the data matrices %%%%%%%%%

if lmon_A && ismember(2,fignums)
    fprintf('\nAssuming that vapor %s exists as monomers, i.e. (1%s,n%s) cluster formation is negligible\n',A,A,B)
    fprintf('#### DeltaG figures ####: This may NOT be a good assumption for e.g. H2SO4-amine systems;\ninstead, you should run ACDC to solve how large fraction of H2SO4 is clustered\n')
end

for nCb=1:length(Cb_vector)

    Cb=Cb_vector(nCb);
    % Partial pressure of B [Pa]
    if strcmpi(B,'W')
        Pb=Cb/100*exp(-2991.2729*temp^(-2)-6017.0128/temp+18.87643854-0.028354721*temp...
            +0.17838301e-4*temp^2-0.84150417e-9*temp^3+0.44412543e-12*temp^4+2.858487*log(temp));
        disp('Water: Assuming that Cb is given as RH (%)')
        %Pb=Cb*1e6*kB*temp; % cm^-3 -> Pa
        %disp('Assuming that Cb is given in cm^-3')
    else
        if l_ppt_B
            Pb=Cb*1e-12*pres_atm; % ppt -> Pa
        else
            Pb=Cb*1e6*kB*temp; % cm^-3 -> Pa
        end
    end

    for nCa=1:length(Ca_vector)

        if lmon_A
            Ca_CIMS=nan;
            Ca=Ca_vector(nCa);
        else
            Ca_CIMS=Ca_vector(nCa);
            % Read the concentration file from ACDC to find out the true acid monomer concentration
            file_neutral=[dirpath,'/c_neutral_A',sprintf('%0.2e',Ca_vector(nCa)),...
                '_',B,num2str(Cb_vector(nCb)),'_I',num2str(IPR),suffix,'.txt'];
            if ~exist(file_neutral,'file')
                % This might be some sprintf conversion problem
                file_neutral=regexprep(file_neutral,'e\+0','e\+00');
            end
            c_neutral=load(file_neutral);
            % Acid monomer should be on the 2nd row, 1st column (cm^-3)
            Ca=c_neutral(2,1);
        end

        % Partial pressure of A (Pa)
        if l_ppt_A
            Pa=Ca*1e-12*pres_atm;
        else
            Pa=Ca*1e6*kB*temp;
        end

        % DeltaG surface at the given concentrations
        deltag_act=nan(size(deltag_ref));

        deltag_ref_mat=nan(maxA+1,maxB+1); % Reference and actual DeltaG matrices
        deltag_act_mat=nan(maxA+1,maxB+1); % (+1 because of the elements with 0 molec.)
        
        evap_tot_mat=nan(maxA+1,maxB+1); % Total evaporation frequency (s^-1)
        monA_coll_mat=nan(maxA+1,maxB+1); % Monomer collision frequencies (s^-1)
        monB_coll_mat=nan(maxA+1,maxB+1);
        
        clust_mat=cell(maxA+1,maxB+1); % Names of the clusters


        for nc=1:length(num{ch})

            lvalid=0;
            nclust=num{ch}(nc);
            [molnames,nmols]=parse_cluster(clust{nclust});
            
            % First see if the cluster contains the possible extra molecule types
            if l_inc_extra
                for nmol=1:length(molnames_extra)
                    if ~ismember(molnames_extra{nmol},molnames)
                        lvalid0=0;
                        break
                    else
                        if nmols(strcmp(molnames,molnames_extra{nmol})) ~= nmols_extra(nmol)
                            lvalid0=0;
                            break
                        else
                            lvalid0=1;
                        end
                    end
                end
                if ~lvalid0
                    continue
                end
            end
            
            % Then see if A and/or B and a possible charge are present
            if (length(molnames)==ntypes_max && ismember(A,molnames) && ismember(B,molnames)) || ...
                    (length(molnames)==ntypes_max-1 && (ismember(A,molnames) || ismember(B,molnames))) || ...
                    (length(molnames)==ntypes_max-2 && ismember(name_ch,molnames) && ~l_ch_nonmol_AB && ~(ch_extra > 1))
                lvalid=1;
            end

            if ch > 1  && ~ismember(name_ch,molnames)
                lvalid=0;
            end
            
            if ~lvalid
                continue
            end
            
            % nmolsA and nmolsB are for the DeltaG conversion,
            % matrow and matcol are for the axes of the plots

            nmolsA=nmols(strcmp(molnames,A));
            nmolsB=nmols(strcmp(molnames,B));

            if isempty(nmolsA), nmolsA=0; end
            if isempty(nmolsB), nmolsB=0; end
            
            % Note that for the DeltaG conversion of charged clusters, only
            % neutral molecules are considered (i.e. "the ion is not a molecule")
            matrow=nmolsA+1;
            matcol=nmolsB+1;
            if ch > 1
                if strcmp(ch_carrier,A)
                    if l_ch_nonmol_AB
                        nmolsA=nmolsA-1;
                    else
                        matrow=matrow+1;
                    end
                end
                if strcmp(ch_carrier,B)
                    if l_ch_nonmol_AB
                        nmolsB=nmolsB-1;
                    else
                        matcol=matcol+1;
                    end
                end
            end
            
            if matrow <= maxA+1 && matcol <= maxB+1
                
                % Convert to the actual partial pressures
                deltag_act(nclust)=deltag_ref(nclust)-kB*temp*(nmolsA*log(Pa/pres)+nmolsB*log(Pb/pres))*J2kcal;
                
                if l_inc_extra
                    for nmol=1:length(molnames_extra)
                        nmols_tmp=nmols_extra(nmol);
                        if molnames_extra{nmol} == name_ch
                            continue
                        elseif l_ch_nonmol_extra && molnames_extra{nmol} == ch_carrier
                            nmols_tmp=nmols_tmp-1;
                        end
                        if l_ppt_extra
                            Pe=Ce(nmol)*1e-12*pres_atm;
                        else
                            Pe=Ce(nmol)*1e6*kB*temp;
                        end
                        deltag_act(nclust)=deltag_act(nclust)-kB*temp*nmols_tmp*log(Pe/pres)*J2kcal;
                    end
                end
                
                % Fill the matrices
                clust_mat{matrow,matcol}=clust{nclust};
                deltag_ref_mat(matrow,matcol)=deltag_ref(nclust);
                deltag_act_mat(matrow,matcol)=deltag_act(nclust);
                
                clustname=clust{nclust};
                if calc_mols(clustname,labels_ch)>1
                    if isnan(get_rate_coefs(clustname,'evap'))
                        % The molecules might be in different orders in the energy file and
                        % in the ACDC system input file
                        for nc0=1:length(clust_acdc)
                            if compare_clusters(clustname,clust_acdc{nc0})
                                clustname=clust_acdc{nc0};
                                break
                            end
                        end
                    end
                    evap_tot_mat(matrow,matcol)=sum(get_rate_coefs(clustname,'evap'),'omitnan');
                    if evap_tot_mat(matrow,matcol)==0
                        evap_tot_mat(matrow,matcol)=nan;
                    end
                    monA_coll_mat(matrow,matcol)=Pa/kB/temp*get_rate_coefs({clustname,['1',A]},'coll');
                    monB_coll_mat(matrow,matcol)=Pb/kB/temp*get_rate_coefs({clustname,['1',B]},'coll');
                end

            end

        end
             
%%%%%%%%% Plot the figures %%%%%%%%%
        
        str_unit_A=' cm^{-3}'; str_unit_B=str_unit_A;
        if l_ppt_A, str_unit_A=' ppt'; end
        if l_ppt_B, str_unit_B=' ppt'; end
        
        if lmon_A
            str_A=['[',A,']_{mon} = ',get_es_str(Ca,0,[1e-3, 1e3]),str_unit_A,', '];
        else
            str_A=sprintf(['[',A,']_{sum} = ',get_es_str(Ca_CIMS,0,[1e-3, 1e3]),str_unit_A,...
                ', [',A,']_{mon} = ',get_es_str(Ca,0,[1e-3, 1e3]),str_unit_A,',\n']);
        end
        
        info_str=[str_A,'[',B,'] = ',get_es_str(Cb,0,[1e-3, 1e3]),str_unit_B,', {\itT} = ',num2str(temp),' K'];
        
        Lwidth=3;
        Fsize=12;
        
        hadd=0; % Numbers for the figure handles (windows): fignums+hadd
        if l_add_figs
            hnums=get(findobj('Type','figure'),'Number');
            if ~isempty(hnums)
                if ~isnumeric(hnums)
                    hnums=cell2mat(hnums);
                end
                hadd=max(hnums);
            end
        end

        for nfig=fignums
            
            if exist('rh','var') && nfig < 3
                continue
            end
            
            if nfig == 1
                data=deltag_ref_mat';
                data_label='\Delta{\itG}_{ref} (kcal mol^{-1})';
                title_str=['Standard \Delta{\itG}_{ref} at {\itT} = ',num2str(temp),' K'];
            elseif nfig == 2
                data=deltag_act_mat';
                data_label='Actual \Delta{\itG} (kcal mol^{-1})';
                title_str=['Actual \Delta{\itG} at ',info_str];
            elseif nfig == 3
                data=log10(evap_tot_mat');
                data_label='Total evaporation rate \Sigma{\it\gamma} (s^{-1})';
                title_str=['Overall evaporation at {\itT} = ',num2str(temp),' K'];
            elseif nfig == 4
                data=log10(monA_coll_mat'./evap_tot_mat');
                data_label=['{\it\beta}_{',A,'}{\itC}_{',A,'} / \Sigma{\it\gamma}'];
                title_str=['Coll. ',A,' / total evap. at ',info_str];
            elseif nfig == 5
                data=log10(monB_coll_mat'./evap_tot_mat');
                data_label=['{\it\beta}_{',B,'}{\itC}_{',B,'} / \Sigma{\it\gamma}'];
                title_str=['Coll. ',B,' / total evap. at ',info_str];
            end
            pcolor_mat=[data nan(size(data,1),1); nan(1,size(data,2)+1)];
            
            figure(nfig+hadd)
            if ~l_add_figs
                clf(nfig+hadd)
            end
            figpos=[50+520*(1-mod(nfig,2)) 500-450*floor(min(nfig,3)/3) 500 450];
            if nfig > 4, figpos(1)=figpos(1)+520*(nfig-3); end
            set(gcf,'Position',figpos);
            %set(gca,'Position',[0.1 0.12 0.8 0.8]);
            set(gca,'Position',[0.1 0.18 0.8 0.75]);
            set(gca,'FontSize',Fsize,'FontWeight','demi');
            set(gcf,'Color','white');
            h=pcolor(pcolor_mat);
            %set(h,'LineWidth',2);
            if nfig <=2
                colormap('cool')
            else
                colormap('parula')
                brighten(0.7)
            end
            cb=colorbar('FontSize',Fsize,'FontWeight','demi');
            ylabel(cb,data_label,'LineWidth',Lwidth,'FontSize',Fsize,'FontWeight','demi')
            if nfig == 3 || nfig == 4 || nfig == 5
                ytickdata=get(cb,'YTick');
                set(cb,'YTickLabel',cellfun(@(x) get_es_str(10^x,0),num2cell(ytickdata),'UniformOutput',false))
            end
            set(gca,'XTick',(1:maxA+1)+0.5,'XTickLabel',0:maxA,'FontSize',Fsize,'FontWeight','demi');
            set(gca,'YTick',(1:maxB+1)+0.5,'YTickLabel',0:maxB,'FontSize',Fsize,'FontWeight','demi');
            if ch == 2
                txt_str='^-';
            elseif ch == 3
                txt_str='^+';
            end
            if ch_carrier == A
                xlabel(['Molecules ',A,' including ',A,txt_str],'FontSize',Fsize,'FontWeight','demi')
            else
                xlabel(['Molecules ',A],'FontSize',Fsize,'FontWeight','demi')
            end
            if ch_carrier == B
                ylabel(['Molecules ',B,' including ',B,txt_str],'FontSize',Fsize,'FontWeight','demi')
            else
                ylabel(['Molecules ',B],'FontSize',Fsize,'FontWeight','demi')
            end
            t=title(title_str,'FontSize',Fsize,'FontWeight','demi');
            set(t,'Position',get(t,'Position')+[0.5 0 0])
            if l_inc_extra
                txt_str=['All clusters contain also ',chem_label(clust_extra,labels_ch,0)];
                annotation('textbox',[0.22 0.015 0.8 0.06],'String',txt_str,'FontSize',Fsize,'FontWeight','demi',...
                    'Color','r','LineStyle','none')
            end
            
            % Print the data values inside the cells
            pos=get(gca,'Position');
            [rows,cols]=size(data);
            width=pos(3)/cols;
            height =pos(4)/rows;
            for i=1:cols
                for j=rows:-1:1
                    if ~isnan(data(j,i))
                        if nfig == 3 || nfig == 4 || nfig == 5
                            str=get_es_str(10^data(j,i),0);
                        else
                            str=sprintf('%0.2f',data(j,i));
                        end
                        annotation('textbox',[pos(1)+width*(i-1),pos(2)+height*(j-1),width,height], ...
                            'string',str,'LineStyle','none','HorizontalAlignment', ...
                            'center','VerticalAlignment','middle','FontWeight','demi','FontSize',Fsize);
                    end
                end
            end
            
        end

        if length(Ca_vector)>1 || length(Cb_vector)>1
            pause
        end

    end

end

end