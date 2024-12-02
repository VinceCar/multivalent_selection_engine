%% This is a Matlab Script that will be use as engine for Selex simulation
%frequency_of_output=1000;
tic
frequency_of_output=100;
input_simulations_06_09_2024_multivalent
number_of_events = 2;
%% Some checks and setting random generator seed
if ~flag_continue
    flag_ok=true;
    check_inputs_and_set_rand_generator
    if ~flag_ok
        return
    end
    current_sweep=0; % variable keeping track of the sweeps
    current_step=0; % variable keeping track of the singe steps
    current_cycle=0; % variable keeping track of the current cycle

    %% Generate the particles and receptors following a specific model
    % The matrix will contain: Favourite targettable per unit (as many as valency), the actual kD (as many as valency), the
    % id of the receptor that a part of the particle is bound (as many parts as
    % valency)
    particle=zeros(Number_of_particles,Valency,3); % Row, specific particle, column is the specific subunit, z 1 is binder family, z 2 is the actual kD per subunit, z 3 is the actual length per subunit;
    particle_status =zeros(Number_of_particles,3); % Row is a specific particle, column 1 is the normalized kon of the particle, column 2 is the normalized multivalent kD and column 3 is the receptor the particle it is bound to
    Receptor=zeros(Number_of_receptors,1); % Row is specific receptor, the column contains the id of a particle it is bound to. it is redundant but kept as a check
    gaussian_model_ref_Wang_et_al % generate the ligands within the particles
    generate_particle_level_properties
    if ~flag_ok
        return
    end
    %% Execute the actual simulations
end
%% Actual selex procedure
while current_cycle<maximum_cycles
    flag_continue=true(1);
    simulation_cycle_end_name=[prefix_for_name,'cycle_',num2str(current_cycle)];
    single_cycle_multivalent_extend_clean
    toc
    current_cycle=current_cycle+1; % move forward cycle
    current_sweep=0; % restart sweep
    current_step=0; % restart steps

    max_sweeps_per_cycle_backup=max_sweeps_per_cycle;
    removing_not_bound_particles
    max_sweeps_per_cycle=max_sweeps_per_cycle_removed_particles;
    frequency_of_output_backup=frequency_of_output;
    frequency_of_output=10;
    single_cycle_multivalent_extend_clean
    
    current_sweep=0; % restart sweep
    current_step=0; % restart steps
    frequency_of_output=frequency_of_output_backup;
    max_sweeps_per_cycle=max_sweeps_per_cycle_backup;
    amplification_smooth_3D
    save([simulation_cycle_end_name,'.mat'])
    if(min(particle_status(:,2)))<0
        return;
    end
end
