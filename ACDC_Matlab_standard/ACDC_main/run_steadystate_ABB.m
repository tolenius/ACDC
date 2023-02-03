function [J,clust,conc,conc_An_tot,Ca_vector,Cb_vector,main_routes] = run_steadystate_ABB(fn_input,varargin)
% Syntax:
% run_steadystate_ABB(fn_input);
% [J,clust,conc,conc_An_tot,Ca_vector,Cb_vector,main_routes] = run_steadystate_ABB(fn_input,'var1',value1,'var2',value2...);
% For options for variables 'var' and the subsequent values 'value', see the cell array var_in in this function
%
% Description:
% Get the simulated steady-state cluster concentrations, formation rate J, and growth pathways
% for a 1-, 2- or 3-component system at given conditions
% Note that here J corresponds to flux of clusters growing out of the given simulation system
% (i.e. the definition of J is not unambiguous, and J depends on the system size - keep this in mind when
% comparing to experiments or other models)
%
% Input is set in the separate input file given as the first input argument
%
% Some variables can also be given as additional input arguments (see the syntax above);
% in this case the corresponding variables in the input file are not used
%
% Output variables can also be returned (instead of only plotting/saving variables in a .mat file);
% see the current output in the function definition above (or add/remove variables according to your needs)


addpath('../Matlab_general','-end')
get_constants;


% Run the input file
[~,~,fext]=fileparts(fn_input);
if ~strcmp(fext,'.m')
    error('The input file must be a .m file')
end
run(fn_input);


% See if some variables are given as input, and fill in/override them
% Same variables are expected in fn_input, except for add_perl_opt_append which is appended to the Perl call
% (while add_perl_opt overrides the add_perl_opt of fn_input)
var_in={'inputfile', 'Gfile', 'DPfile', 'temp', 'rh', 'pw', 'IPR', 'cs_ref', 'dirpath', 'suffix', 'add_perl_opt', 'add_perl_opt_append'};

if ~isempty(varargin)
    for i=1:2:length(varargin)
        if ismember(varargin{i},var_in)
            % eval is not a good practice, but it's easiest for now
            eval([varargin{i},'=varargin{i+1};'])
        else
            if isnumeric(varargin{i}), str=num2str(varargin{i}); else, str=varargin{i}; end
            help run_steadystate_ABB
            error(['Cannot recognize input variable ''',str,''' - see the function information above'])
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    Generate the ACDC Matlab files    %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create the command line string to run the Perl script
ACDCfile=find_newest_file('acdc.pl');
commandstr=['perl ',ACDCfile,' --i ',inputfile,' --e ',Gfile];

% Possible dipole moment file
if exist('DPfile','var')
    commandstr=[commandstr,' --dip ',DPfile];
end

% Temperature (in case H and S are used instead of G)
commandstr=[commandstr,' --temperature ',num2str(temp)];

% RH or partial pressure of water
if exist('rh','var')
    commandstr=[commandstr,' --rh ',num2str(rh)];
end
if exist('pw','var')
    commandstr=[commandstr,' --pw ',sprintf('%.14e',pw)];
end

% Possible loss terms
loss_commandstr='';
loss_suffixstr='';
loss_infostr='';

if lexp_loss == 1
    loss_commandstr=[loss_commandstr,' --use_cs --cs exp_loss --exp_loss_coefficient ',num2str(cs_ref),' --exp_loss_exponent ',num2str(cs_exp)];
    if exist('cs_ref_str','var')
        loss_commandstr=[loss_commandstr,regexprep(cs_ref_str,'^--',' --')];
    elseif exist('B','var')
        % Use approximate diameter of H2SO4, if no specific monomer is given as the reference size
        loss_commandstr=[loss_commandstr,' --exp_loss_ref_size 0.55'];
    end
    loss_suffixstr=[loss_suffixstr,'_cs_exp',num2str(cs_exp),'_ref',sprintf('%.6e',cs_ref)];
    loss_infostr=[loss_infostr,', CS_{ref} = ',get_es_str(cs_ref,0),' s^{-1}'];
end
if lCLOUD_loss == 1
    loss_commandstr=[loss_commandstr,' --use_wl --wl CLOUD4_simple --use_dilution'];
    loss_suffixstr=[loss_suffixstr,'_wl_CLOUD'];
    loss_infostr=[loss_infostr,', CLOUD wall losses'];
end

if loss_commandstr ~= ""
    
    commandstr=[commandstr,loss_commandstr];
    
    loss_suffixstr=loss_suffixstr(2:end);
    loss_infostr=loss_infostr(3:end);
else
    loss_suffixstr='0';
    loss_infostr='no losses';
end

% Additional options
if exist('add_perl_opt','var')
    commandstr=[commandstr,regexprep(add_perl_opt,'^--',' --')];
end
if exist('add_perl_opt_append','var')
    commandstr=[commandstr,regexprep(add_perl_opt_append,'^--',' --')];
end

% Run the Perl script
fprintf(['\nGenerating the equations:\n',regexprep(commandstr,'\s+--','\n--'),'\n\n'])
lrun=system(commandstr);
if lrun ~= 0
    error('Running the Perl script failed!')
end
rehash pathreset


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    Run the ACDC Matlab simulation    %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the simulated cluster and ion names
fid=fopen('driver_acdc.m', 'r');
line=0;
while line ~= -1
    line=[fgetl(fid),' '];
    if ~isempty(regexp(line,'^\s*clust\s*=\s*{', 'once'))
        eval(line)
    elseif ~isempty(regexp(line,'^\s*labels_ch\s*=\s*{', 'once'))
        eval(line);
    end
    if exist('clust','var') && exist('labels_ch','var')
        break
    end
end
fclose(fid);

% Sort the clusters according to the charge
ndistr=1;
n_ch=cell(1,3);
for nc=1:length(clust)
    ch_temp=solve_charge(clust{nc},labels_ch);
    n_ch{ch_temp}=[n_ch{ch_temp} nc];
end

if IPR > 0
    if isempty(n_ch{2}) && isempty(n_ch{3})
        fprintf('\nNo charged clusters found - not including ions\n')
        IPR=0;
    else
        ndistr=3;
    end
end


% Find the indices of clusters containing 2, 3, 4, ... molecules of type A
% Previously used variable which only considered dimers
if exist('ldimers','var')
    if ldimers
        lnmers=1;
    end
else
    ldimers=0;
end
if lnmers
    n_An=cell(1);
    for nc=1:length(clust)
        % Consider only neutral clusters
        if solve_charge(clust{nc},labels_ch) == 1
            [molnames,nmols]=parse_cluster(clust{nc});
            if ismember(A,molnames)
                nmols_A=nmols(strcmp(A,molnames));
                if length(n_An) < nmols_A
                    n_An{nmols_A}=[];
                end
                n_An{nmols_A}=[n_An{nmols_A} nc];
            end
        end
    end
    maxn_An=length(n_An);
    if ldimers
        plot_An=2;
    else
        plot_An=1:maxn_An;
    end
end

% Plot or save growth routes
% Previously only plotting of fluxes was possible
if ~exist('lroute','var')
    if lfluxes
        lroute=1;
    else
        lroute=0;
    end
end

% Draw the given figures (or only calculate the data)
if ~exist('lshow','var')
    lshow=1;
end

%%%%%%%%% Create the vectors and matrices etc. needed  %%%%%%%%%
%%%%%%%%% in the simulation and for saving the results %%%%%%%%%

if ~exist('B','var')
    B={};
end

if isempty(B)
    Cb_vector{1}=[0];
    Cb_vector_ppt{1}=[0];
    Cb_limits{1}=[0 0];
    Cb_points{1}=1;
    
    % If there is only compound A, its input concentration must correspond to monomers
    incl_with_mon={'', '', ''};
end

if length(B) < 2
    Cb_vector{2}=[0];
    Cb_vector_ppt{2}=[0];
    Cb_limits{2}=[0 0];
    Cb_points{2}=1;
end

if all(cellfun(@isempty,incl_with_mon)) && ~isempty(B)
    fprintf('\nUsing true monomer concentration for all vapors - do you really want this?\n')
end

% Create the vapor concentration vectors, if vectors of separate values are not used
if lCa_fixed == 0
    Ca_vector=10.^(log10(Ca_limits(1)):(log10(Ca_limits(end))-log10(Ca_limits(1)))/(Ca_points-1):log10(Ca_limits(end)));
    Ca_vector(1)=Ca_limits(1); % Change back the 1st and last value, otherwise they can be something slightly different
    Ca_vector(end)=Ca_limits(end);
else
    Ca_limits=[Ca_vector(1) Ca_vector(end)];
end

for nmolB=1:length(B)
    if lCb_fixed{nmolB} == 0
        Cb_vector{nmolB}=10.^(log10(Cb_limits{nmolB}(1)):(log10(Cb_limits{nmolB}(end))-log10(Cb_limits{nmolB}(1)))/...
            (Cb_points{nmolB}-1):log10(Cb_limits{nmolB}(end)));
        Cb_vector{nmolB}(1)=Cb_limits{nmolB}(1);
        Cb_vector{nmolB}(end)=Cb_limits{nmolB}(end);
    else
        Cb_limits{nmolB}=[Cb_vector{nmolB}(1) Cb_vector{nmolB}(end)];
    end
    % Convert the B concentration from ppt to cm^-3, if needed
    if lB_ppt{nmolB} == 1
        Cb_vector_ppt{nmolB}=Cb_vector{nmolB};
        Cb_limits_ppt{nmolB}=Cb_limits{nmolB};
        Cb_vector{nmolB}=Cb_vector_ppt{nmolB}*1e-18*pres_atm/kB/temp; % cm^-3
        Cb_limits{nmolB}=Cb_limits_ppt{nmolB}*1e-18*pres_atm/kB/temp;
    else
        Cb_vector_ppt{nmolB}=Cb_vector{nmolB}/(1e-18*pres_atm/kB/temp);
    end
end

% Create matrices/cells for the distributions and formation rate
conc = cell(length(Ca_vector),length(Cb_vector{1}),length(Cb_vector{2}));
J = nan(length(Ca_vector),length(Cb_vector{1}),length(Cb_vector{2}));
if lnmers
    conc_An_tot=cell(1,maxn_An);
    for nmols_A=1:maxn_An
        conc_An_tot{nmols_A} = nan(length(Ca_vector),length(Cb_vector{1}),length(Cb_vector{2}));
    end
else
    conc_An_tot=cell(1,1);
end
main_routes = cell(length(Ca_vector),length(Cb_vector{1}),length(Cb_vector{2}));

if save_data==1
    % .mat file in which the results will be saved
    %matfile=[dirpath,'C_J_',A,B{:},suffix,'.mat'];
    matfile=[dirpath,'J',suffix,'.mat'];
    % Create a new folder for the results, if needed
    if ~exist(dirpath,'dir')
        mkdir(dirpath)
    end
end

% Legend
Cb_leg=cell(1,length(Cb_vector{1})*length(Cb_vector{2}));
nleg=0;

% Settings related to effective vapor concentrations
vapor_label={'', '', ''};
for nmol=1:length(B)+1

    if nmol == 1
        vapor_label{nmol}=['[',A,']'];
    else
        vapor_label{nmol}=['[',B{nmol-1},']'];
    end

    % Include possible other complexes contributing to effective/measurable vapor concentration
    if ~isempty(incl_with_mon{nmol})
        
        % See if the given clusters that contribute are actually included
        clust_tmp=strsplit(incl_with_mon{nmol},'-');
        clust_tmp=clust_tmp(~cellfun('isempty',clust_tmp));

        % Rewrite the string for contributing clusters depending on
        % what is included and on the cluster name label convention
        str_tmp='';
        for nct=1:length(clust_tmp)
            lfound=0;
            for nc=1:length(clust)
                if compare_clusters(clust{nc},clust_tmp{nct})
                    str_tmp=[str_tmp,'-',clust{nc}];
                    lfound=1;
                    break
                end
            end
            if lfound==0
                disp([clust_tmp{nct},' is not in the cluster set, excluding it in the measurable vapor'])
            end
        end
        incl_with_mon{nmol}=str_tmp;
        
        % Add note to legends etc.
        if ~isempty(incl_with_mon{nmol})
            vapor_label{nmol}=[vapor_label{nmol},'_{eff}'];
        end
        
    end
    
end

%%%%%%%%% Loops over all the conditions %%%%%%%%%

%fprintf(['\nBeginning simulations for input files\n',Gfile,'\n',inputfile,'\n'])

% Looping through concentrations of B
for nCb1=1:length(Cb_vector{1})

    nCb{1} = nCb1;
    Cb{1} = Cb_vector{1}(nCb1);

    for nCb2=1:length(Cb_vector{2})

        nCb{2} = nCb2;
        Cb{2} = Cb_vector{2}(nCb2);
        
        % Legend for concentrations of B
        str='';
        for nmolB=1:length(B)
            if lB_ppt{nmolB} == 1
                str=[str,vapor_label{nmolB+1},' = ',num2str(Cb_vector_ppt{nmolB}(nCb{nmolB})),' ppt, '];
            else
                str=[str,vapor_label{nmolB+1},' = ',get_es_str(Cb_vector{nmolB}(nCb{nmolB}),0),' cm^{-3}, '];
            end
        end
        str=str(1:length(str)-2);
        nleg=nleg+1;
        Cb_leg{nleg}=str;
        
        % Looping through concentrations of A
        for nCa=1:length(Ca_vector)

            Ca = Ca_vector(nCa);
            
            %fprintf('\n')
            %disp([vapor_label{1},' = ',sprintf('%0.2e',Ca),' cm^{-3}, ',str]);

            % Setting the input concentrations/sources
            fid=fopen('sources.txt','w');
            
            for nmol=1:length(B)+1
                if nmol == 1
                    fprintf(fid,['constant 1',A,' %e'],Ca);
                else
                    fprintf(fid,['constant 1',B{nmol-1},' %e'],Cb{nmol-1});
                end
                % Include possible other complexes contributing to effective/measurable vapor concentration
                if ~isempty(incl_with_mon{nmol})
                    fprintf(fid,[' ',incl_with_mon{nmol}]);
                end
                fprintf(fid,'\n');
            end
            
            if IPR > 0
                fprintf(fid,'source neg %e\n',IPR);
                fprintf(fid,'source pos %e\n',IPR);
                fprintf(fid,'wall enhancement %e\n',Iwlf);
            end
            
            fclose(fid);

            % Running ACDC
            % Initializing
            [C, T, ok, clust, Cf0]=driver_acdc(1e-8,'Sources_in','sources.txt');
            [C, T, ok, clust, Cf0]=driver_acdc(1e-1,Cf0,'Sources_in','sources.txt');
            %[C, T, ok, clust, Cf0]=driver_acdc(1e-1,C,T,'Sources_in','sources.txt'); % For saving the full time series
            Cf = Cf0;
            % Checking for negative concentrations
            if ok == -1
                fprintf('Negative concentrations for Ca = %.2e, Cb = %.2e, %.2e cm^-3\n',Ca,Cb{:})
                continue
            end
            % Running the solver until the concentrations have converged
            % Driver input and output in string format (for an easier inclusion of different simulation options, although
            % using eval is not really recommended)
            driver_input_str='10^j,Cf,''Sources_in'',''sources.txt'',''repeat''';
            %driver_input_str='10^j,Cf,''Sources_in'',''sources.txt'',''repeat'',''Options'',''''''AbsTol'''',1e-20,''''RelTol'''',1e-10''';
            %driver_input_str='10^j,C,T,''Sources_in'',''sources.txt'''; % For saving the full time series
            driver_output_str='C, T, ok, clust, Cf, labels_ch, clust_flux, J_out';
            if lfluxes || lroute
                driver_output_str=[driver_output_str,', flux, outflux'];
            end
            %for j=[4*ones(1,5) 3*ones(1,5)]
            for j=[4*ones(1,3) 3*ones(1,20)]
                eval(['[',driver_output_str,']=driver_acdc(',driver_input_str,');'])
                if ok==1, break, end
            end
            % Checking convergence
            if ok ~= 1
                fprintf('Running steady state for Ca = %.2e, Cb = %.2e, %.2e cm^-3 failed!\n',Ca,Cb{:})
                continue
            end

            % Save the distribution and J
            conc{nCa,nCb1,nCb2} = Cf;
            J(nCa,nCb1,nCb2) = J_out(end);

            % Plot the fluxes and/or save the main growth route
            if lfluxes || lroute
                main_routes_tmp={};
                nch=1;
                [~,main_routes_tmp{nch},~]=track_fluxes(clust_flux,clust,labels_ch,flux,outflux,'ch',nch,'ldisp',lfluxes);
                if IPR > 0
                    for nch=2:3
                        if lfluxes, pause, end
                        [~,main_routes_tmp{nch},~]=track_fluxes(clust_flux,clust,labels_ch,flux,outflux,'ch',nch,'ldisp',lfluxes);
                    end
                end
                main_routes{nCa,nCb1,nCb2} = main_routes_tmp;
            end


            %%%%%%%%% Additional tasks %%%%%%%%%

            if ldistribution
                % Do a bar plot for the current vapor concentrations
                for nd=1:ndistr
                    figure(20-nd+1)
                    clf(20-nd+1)
                    set(gca,'YScale','log')
                    set(gcf,'Color','white')
                    box on
                    hold on

                    nc_temp=n_ch{nd};
                    
                    %fprintf('Total concentration for charging state %d: %.2e cm^-3\n',nd,sum(Cf(nc_temp)))

                    bar(Cf(nc_temp))
                    xlim([0,length(nc_temp)+1])
                    set(gca,'XTick',1:length(nc_temp),'XTickLabel',clust(nc_temp));

                    ax=gca;
                    ax.XTickLabelRotation=45;
                    ylabel('{\itC}_{steady-state} (cm^{-3})')
                    str_tmp=[vapor_label{1},' = ',get_es_str(Ca,0),' cm^{-3}, '];
                    if ~isempty(B), str_tmp=[str_tmp,Cb_leg{nleg},', ']; end
                    str_tmp=[str_tmp,'{\itT} = ',num2str(temp),' K, ',loss_infostr];
                    title(str_tmp,'FontWeight','normal')
                    drawnow
                end
            end

            if lnmers
                % Save the n-mer concentration
                for nmols_A=plot_An
                    conc_An_tot{nmols_A}(nCa,nCb1,nCb2)=sum(Cf(n_An{nmols_A}),'omitnan');
                end
            end

            % Update the results
            if save_data==1
%                 save(matfile,'clust','conc','conc_An_tot','J','Ca_vector','Cb_vector',...
%                     'temp','loss_infostr','IPR','Gfile','inputfile')
                save(matfile,'J','Ca_vector','Cb_vector','Cb_vector_ppt',...
                    'temp','IPR','cs_ref','cs_exp','loss_infostr',...
                    'B','Gfile','inputfile')
            end
            
            if lfluxes || ldistribution
                if nCa<length(Ca_vector) || nCb{1}<length(Cb_vector{1}) || nCb{2}<length(Cb_vector{2})
                    pause
                end
            end
        
        end

    end
    
end

%%%%%%%%% Figures %%%%%%%%%

if lshow

    title_str=['{\itT} = ',num2str(temp),' K, ',loss_infostr];
    if IPR > 0
        title_str=[title_str,', IPR = ',num2str(IPR),' cm^{-3}s^{-1}'];
    end
    
    xlabel_str=[vapor_label{1},' (cm^{-3})'];

    if lJ
        % Plot the formation rate out of the system at the different A and B concentrations
        figure(21)
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        set(gcf,'Color','white')
        box on
        hold on; set(gca,'ColorOrderIndex',1)
        for nCb1=1:length(Cb_vector{1})
            for nCb2=1:length(Cb_vector{2})
                plot(Ca_vector,J(:,nCb1,nCb2),LStyle)
            end
        end
        xlim([Ca_limits(1),Ca_limits(end)])
        xlabel(xlabel_str)
        ylabel('{\itJ} (cm^{-3} s^{-1})')
        title(['Steady-state {\itJ}, ',title_str],'FontWeight','normal')
        if ~isempty(B), legend(Cb_leg,'location','best'), end
    end

    if lnmers
        % Plot the sum of An clusters (An, AnB1, AnB2, ...) at the different A and B concentrations
        sum_str=['\Sigma [',A];
        for nmolB=1:length(B)
            sum_str=[sum_str,'\cdot',B{nmolB},'_{\itn}'];
            if length(B)>1
                sum_str=[sum_str,'_{',num2str(nmolB),'}'];
            end
        end
        sum_str=[sum_str,'] (cm^{-3})'];

        for nmols_A=plot_An
            figure(22+nmols_A-1)
            set(gca,'XScale','log')
            set(gca,'YScale','log')
            set(gcf,'Color','white')
            box on
            hold on; set(gca,'ColorOrderIndex',1)
            for nCb1=1:length(Cb_vector{1})
                for nCb2=1:length(Cb_vector{2})
                    plot(Ca_vector,conc_An_tot{nmols_A}(:,nCb1,nCb2),LStyle)
                end
            end
            xlim([Ca_limits(1),Ca_limits(end)])
            xlabel(xlabel_str)
            ylabel(strrep(sum_str,A,[A,'_{',num2str(nmols_A),'}']))
            title(['Steady-state concentrations, ',title_str],'FontWeight','normal')
            if ~isempty(B), legend(Cb_leg,'location','best'), end
        end
    end
    
end

end