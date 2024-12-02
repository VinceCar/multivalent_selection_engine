%% This is a script containing directives for the simulations
flag_continue=false(1);
if flag_continue
    file_to_load=''; % insert here the path to the file to load
    load(file_to_load)
else
    change_concentration_sweep=Inf;
    prefix_for_name='caroprese_vincenzo_06_09_2024_trial_1_'; % prefix used to generate files
    target_concentration=1e-7; % in M
    solution_prefactor=target_concentration/1e5; % this is a prefactor that needs to take into account that we are speaking of concentrations and not particles

    seed=42; 
    max_sweeps_per_cycle=1e2; % recommended 1e7
    max_sweeps_per_cycle_removed_particles=1e1;
    maximum_cycles=7; % the cycles of selex
    Number_of_particles=1e5; % 1 e 6 Standard numbers in an experiment
    Number_of_receptors=round(Number_of_particles/3); %  1 e 4 Standard numbers in an experiment
    Valency=3; % We are going to try multivalent selex and analyze monovalent and multivalent selex
    Receptor_targettable_domains=1; % the number of different domanis that we can target


    % we are going to model the possible kDs that you can obtain as normal
    % distributions centered around a certain average
    vector_of_average_kDs=[-10]; % same units, it has been normalized by M
    vector_of_standard_deviations_kDs=[0.85];
    vector_of_average_len=[2.53]; % same units, it has been in nm
    vector_of_standard_deviations_len=[0.8820]; 
    vector_of_support_len=[1.47]; % this is actua    
    target_fraction_of_relative_population=[ones(Receptor_targettable_domains,1)/Receptor_targettable_domains]; % fraction of particles subunits that target a specific domain
    Source_kD_file = 'Source_grids/grids_for_interpolations_2_12.mat';

    cycle_change = Inf;
    New_Valency = 3;
end
