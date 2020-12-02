%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    Settings for rates_and_deltags_ABe.m    %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables beginning with l:
% 1 = condition true, i.e. use / do / plot
% 0 = condition false, i.e. do NOT use / do / plot


%%%%%%%%% Quantities to calculate or plot %%%%%%%%%

% Indices of the wanted figures
% 1: Standard (reference) DeltaG surface
% 2: Actual DeltaG surface at the given vapor concentrations
% 3: Total evaporation rates
% 4: Ratio of collision rate of vapor A to the total evaporation rate
% 5: Ratio of collision rate of vapor B to the total evaporation rate
fignums=[1 2 3 4 5];

% Keep the open figures (from a previous rates_and_deltags_ABe.m call) and add new figure windows,
% or overwrite previous figures
l_add_figs=1;


%%%%%%%%% Cluster set and conditions %%%%%%%%%

% Charging state
ch=1;               % 1 = neutral, 2 = neg., 3 = pos.

% Names of the included compounds
A='A';              % Compound A: Primarily the main driver of clustering, most often an acid
B='N';              % Compound B: Most often a base

% Clusters included in the plot
maxA=5;             % Up to maxA molecules of type A
maxB=5;             % Up to maxB molecules of type B

% Vapor concentrations (cm^-3 or ppt)
Ca_vector=[1e7];
Cb_vector=[100];

l_ppt_A=0;          % Unit: 0 = cm^-3, 1 = ppt
l_ppt_B=1;


% Possible additional molecule type(s)
% Same numbers of additional molecules in all clusters, i.e. they don't affect the (A,B)-grid
% NOTE: If you only want to include a charge, simply set ch to 2 or 3 - no need to set clust_extra to be an ion
l_inc_extra=0;
clust_extra='1D';

% Concentrations of the extra types (for a possible DeltaG plot)
% Each element of Ce corresponds to one additional vapor type (in the order given by clust_extra)
Ce=[0.1];
l_ppt_extra=1;


% Temperature (K)
temp=280;

% RH (percent)
rh=0;


%%%%%%%%% Free energy and ACDC input files (the latter for obtaining the rate constants) %%%%%%%%%

% Create the rate constant files now or use pre-existing files
lrun_perl=1;

% Cluster set file for obtaining the rate constants from ACDC
inputfile=['./Cluster_set_files/input_',A,B,'_neutral_neg_pos.inp'];

% DeltaG / DeltaH, DeltaS file
Gfile='./Energy_and_rate_files/B3LYP_RICC2/HS298.15K_example.txt';

% Dipole moments for a system containing charged clusters (optional; comment out if it doesn't exist)
DPfile='./Energy_and_rate_files/B3LYP_RICC2/dip_pol_298.15K_example.txt';


%%%%%%%%% Miscellaneous %%%%%%%%%

% Assume that the A concentration is the true monomer concentration, or assume that it is the "measurable" concentration
% and use a simulation result to obtain the true monomer concentration
% For assessing collision/evaporation ratios, this doesn't matter, but for precise actual DeltaGs,
% lmon_A should be set to 1 for systems where e.g. (1acid,nbase)-clusters are likely to form
lmon_A=1;

if ~lmon_A
    % Directory path for obtaining the A monomer concentration
    dirpath=['./Growth_',A,B,'_',num2str(temp)];
    suffix='';
    IPR=3;
end
