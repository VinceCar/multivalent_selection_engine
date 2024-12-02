%% Some checks and setting random generator seed
if exist('seed','var')
    if isempty(seed)==false(1)
        rng(seed)
    else
        fprintf("Careful, random seed not initiated. The seed is going to be the ones saved in the variable random_settings");
        random_settings=rng();
    end
else
    fprintf("Careful, random seed not initiated. The seed is going to be the ones saved in the variable random_settings");
    random_settings=rng();
end
if ~exist('prefix_for_name','var')
    fprintf("Error! You did not specify a prefix to save the variables.");
    flag_ok=false(1);
    return
end
if ~exist('max_sweeps_per_cycle','var')
    fprintf("Error! You did not specify how many sweeps should each cycle last.");
    flag_ok=false(1);
    return
end
if ~exist('maximum_cycles','var')
    fprintf("Error! You did not specify how many cycles should the simulation do.");
    flag_ok=false(1);
    return
end
if ~exist('Number_of_particles','var')
    fprintf("Error! You did not specify how many particles should be in the simulation.");
    flag_ok=false(1);
    return
end
if ~exist('Valency','var')
    fprintf("Error! You did not specify the valency of the particles in the simulation.");
    flag_ok=false(1);
    return
end
if ~exist('Number_of_receptors','var')
    fprintf("Error! You did not specify how many receptors should be in the simulation.");
    flag_ok=false(1);
    return
end
if ~exist('Source_kD_file','var')
    fprintf("Error! You did not specify the model to use for multivalent kD estimation.");
    flag_ok=false(1);
    return
end
if ~exist('Receptor_targettable_domains','var')
    fprintf("Error! You did not specify how many targettable domains there are in each subunit of the receptor.");
    flag_ok=false(1);
    return
end
if ~exist('target_fraction_of_relative_population','var')
    fprintf("Error! You did not specify the relative fractions of strands that can target a specific subunit.");
    flag_ok=false(1);
    return
end
if (abs(sum(target_fraction_of_relative_population)-1)>1e-6)
    fprintf("Error! Check the fraction of relative populations, they do not add up to 1");
    flag_ok=false(1);
    return
end
if ~exist('cycle_change','var')
    fprintf("Error! You did not specify at which cycle we switch from mono to multivalent selex. In case no change is desired, set the variable to Inf.");
    flag_ok=false(1);
    return
end
if ~exist('New_Valency','var')
    fprintf("Error! You did not specify the new valency when we switch from mono to multivalent selex. In case no change is desired, set the variable to Inf.");
    flag_ok=false(1);
    return
end