%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    Simulation set-up for run_steadystate_ABB.m    %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables beginning with l:
% 1 = condition true, i.e. use / do / plot
% 0 = condition false, i.e. do NOT use / do / plot


%%%%%%%%% Quantities to calculate or plot %%%%%%%%%

lJ=1;               % Plot the steady-state J
ldistribution=0;    % Plot the steady-state cluster distributions
ldimers=0;          % Plot the total concentration of electrically neutral clusters containing two molecules A (typically an acid)
lfluxes=0;          % Track and plot the growth pathways

save_data=0;        % Save the simulation output (as .mat)


%%%%%%%%% Figure formats %%%%%%%%%

% For plotting different simulation runs in the same figure, it's convenient to run run_steadystate_ABB.m,
% leave the figure windows open, change the line style below, and re-run - the runs will appear with different line styles

LStyle = '-';       % Line style


%%%%%%%%% Simulation system %%%%%%%%%

% Names of the simulated compounds

% NOTE: When studying a 2-component system, simply give only one molecule name in cellstring B,
% i.e. no need to comment out any settings for the possible 3rd compound

A='A';              % Compound A: Primarily the main driver of clustering, most often an acid
B={'N'};            % Additional compound(s) B: Most often bases, 1 or 2 compounds

% Cluster set file
inputfile=['./Cluster_set_files/input_',A,B{:},'_neutral_neg_pos.inp'];

% DeltaG / DeltaH, DeltaS file
Gfile='./Energy_and_rate_files/B3LYP_RICC2/HS298.15K_example.txt';

% Dipole moments for a system containing charged clusters (optional; comment out if it doesn't exist)
DPfile='./Energy_and_rate_files/B3LYP_RICC2/dip_pol_298.15K_example.txt';

% Path for the directory where the results are saved (if they are saved)
dirpath=['./Testfolder_',A,B{:},'/'];

% Suffix for the result files, e.g. you can write some info on the used free energy data etc.
% (if you want no suffix, set this to '')
suffix='';


%%%%%%%%% Simulation conditions %%%%%%%%%

% Temperature (K)
temp=280;

% RH (percent)
rh=0;

% Vapor concentrations

% Use EITHER fixed, separate values, OR give the min. and max. values and
% the number of data points along the (logarithmic) vapor concentration axis

% Concentration of A (normally acid) (% cm^-3)
lCa_fixed=0;        % Use fixed values
if lCa_fixed == 1
    Ca_vector=[1e6 1e8];
else
    Ca_limits=[1e6 1e8];
    Ca_points=10;
end

% Concentration of compounds B (normally bases) (cm^-3 or ppt)
lB_ppt{1}=1;        % Assume that the concentration is in ppt instead of cm^-3
lCb_fixed{1}=1;
if lCb_fixed{1} == 1
    Cb_vector{1}=[10 100 1000]; % ppt
    %Cb_vector{1}=[1e-3 1e-1 1 10]; % ppt
else
    Cb_limits{1}=[10 1000]; % ppt
    Cb_points{1}=10;
end
lB_ppt{2}=1;
lCb_fixed{2}=1;
if lCb_fixed{2} == 1
    Cb_vector{2}=[1e-3 1e-1 1 10]; % ppt
else
    Cb_limits{2}=[1e6 1e8]; % cm^-3
    Cb_points{2}=10;
end


% Possible external sinks of clusters and vapors

% Coagulational scavenging (with exponential size-dependence)
lexp_loss=1;
if lexp_loss == 1
    cs_exp=-1.6;    % Exponent; typical boundary layer value: -1.6
    cs_ref=1e-3;    % Reference loss rate for a vapor (H2SO4) monomer (s^-1); typical boundary layer values: 1e-4...1e-2
end

% CLOUD experiment: wall and dilution losses
lCLOUD_loss=0;


% Settings related to ions (optional; if no ions are included, IPR is discarded)

% Ionization rate (cm^-3 s^-1)
IPR=3;

% Ion wall loss enhancement factor for CLOUD
if lCLOUD_loss == 1
    Iwlf=3.3;
else
    Iwlf=1.0;
end


%%%%%%%%% Miscellaneous %%%%%%%%%

% Additional options for the Perl script
%add_perl_opt=' --no_evap';

% Use the true A monomer concentration, i.e. do not count clusters consisting of one A and n B molecules
lmon_A=0;           % It's better to keep this as 0! (It's expected to correspond to the "measurable" acid concentration)
nB_for_1A=4;        % Maximum number of B molecules attached to the A monomer in case lmon_A=0

