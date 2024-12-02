%% This is a clean up version of the single cycle script that contains all its components
% start the swwp
while current_sweep<max_sweeps_per_cycle
    %% Resetting the binding
    if current_step==0
        if current_sweep==0
            % Initialize useful variables to log the status
            data_points=round(max_sweeps_per_cycle/frequency_of_output);
            total_particle_bound=zeros(data_points,1);
            total_receptor_bound=zeros(data_points,1);
            ligand_bound_per_family=zeros(data_points,Receptor_targettable_domains);
            average_kD_bound=zeros(data_points,1);
            average_kD_bound_per_family=zeros(data_points,Receptor_targettable_domains);

            % Initialize useful variables for the simulation,
            % comparmentalize particles for easy access
            k_on_whole_particles_cumsum=cumsum(particle_status(:,1)); % Initial vector containing the cumulative rates of association
            interval_particle=1; % sample the vector every interval particles
            if Number_of_particles>10000 % when it is more convenient
                interval_particle=round(sqrt(Number_of_particles)); % define an interval bigger than 1
            end
            num_slices_particles = ceil(Number_of_particles / interval_particle); % Number of slices needed
            array_slices_k_on_particles = zeros(interval_particle,num_slices_particles);                % Reorganize particle information for easy access (association constants)
            array_slices_k_off_particles = zeros(interval_particle,num_slices_particles);                % Reorganize particle information for easy access (dissociation constants)
            for i = 1:num_slices_particles % assign the right association constants
                % Calculate the start and end indices for each slice
                start_idx = (i - 1) * interval_particle + 1;
                end_idx = min(i * interval_particle, Number_of_particles);
    
                % Slice the array and store it in the cell array
                array_slices_k_on_particles(1:(end_idx-start_idx+1),i) = particle_status(start_idx:end_idx,1);
            end
            k_on_whole_particles_cumsum=[k_on_whole_particles_cumsum(interval_particle:interval_particle:Number_of_particles);k_on_whole_particles_cumsum(end)]; % Resample the vector of cumulative rates of association
            k_off_whole_particles_cumsum=zeros(size(k_on_whole_particles_cumsum)); % Subsampled vector of cumulative rates of dissociation
            compensation_k_off_whole_particles_cumsum=zeros(size(k_on_whole_particles_cumsum));
            % Comparmentalize receptors for easy access
            interval_receptors=1; % sample the vector every interval receptor
            if Number_of_receptors>10000 % when it is more convenient
                interval_receptors=round(sqrt(Number_of_receptors)); % define an interval bigger than 1
            end
            num_slices_receptors = ceil(Number_of_receptors / interval_receptors); % Number of slices needed
            array_slices_unbound_receptors = zeros(interval_receptors,num_slices_receptors); % a flag that tells you if a receptor is bound or not               
            for i = 1:num_slices_particles % assign the right association constants
                % Calculate the start and end indices for each slice
                start_idx = (i - 1) * interval_receptors + 1;
                end_idx = min(i * interval_receptors, Number_of_receptors);
    
                % Slice the array and store it in the cell array
                array_slices_unbound_receptors(1:(end_idx-start_idx+1),i) = 1;
            end
            
            k_on_whole_receptors_cumsum=(1:Number_of_receptors)';
            k_on_whole_receptors_cumsum=[k_on_whole_receptors_cumsum(interval_receptors:interval_receptors:Number_of_receptors);k_on_whole_receptors_cumsum(end)]; % Resample the vector of cumulative rates of association receptors
        

        end
        %% This is a module in case we want to perform binding assays simulations changing the concentration (not tested with the current setup)
        if (current_step==0) && (mod(current_sweep,change_concentration_sweep)==0)
            particle_status(:,3)=0; % make sure all particles are unbound
            Receptor=zeros(Number_of_receptors,1); % make sure all receptors are unbound
            solution_prefactor=(solution_prefactor*10^(0.1*(current_sweep>0))); % effective change of concentration the 10^(0.1*(current_sweep>0)) has empirically found to produce nice coverage on logarithmic scale
            
            %% set up the bindings
%             id_particles=1:Number_of_particles; % valid because we do not distinguish crosslinking and free receptors
%             id_receptors=randperm(Number_of_receptors); % cause all receptors are equal

            % Now this are variables that we use to compute the rates and
            % will get updated according to different scenarios

            free_particles=sum(particle_status(:,1)>0); % total count of free particles
            free_receptors=Number_of_receptors;
            rate_events=zeros(number_of_events,1); % rate_association_free; %rate_association_free_partially_bound; rate_association_multivalent; rate_association_crosslink; rate_dissociation
            counter_events=rate_events;
            counter_new_binding=0;
            counter_new_binding_receptors=0;
            sum_kD=0; % total kD measuring dissociation
        end
    end

    %% move time forwards
    current_step=current_step+1;
    if mod(current_step,Number_of_receptors)==0
        current_step=0;
        current_sweep=current_sweep+1;
    end
    %% support variables update
    % notice we are exploiting that we expect association constants to
    % always be in simlar value (everything is normalized by 1) but
    % otherwise we need to multiply ka_i.*[conc|_i. same goes for kD and
    % ecc

    free_in_solution_ma=free_particles*free_receptors; % effective activitity depends on free particle encountering free receptors on free subunits.
    rate_events(1)=solution_prefactor*(free_in_solution_ma); % rate_association_free vs free;
    rate_events(2)=sum_kD; % rate_dissociation

    total_scaled_rate=cumsum(rate_events); % put the rates together to do random extraction
    total_scaled_rate=total_scaled_rate./total_scaled_rate(end); % normalization
    ran_numb=rand(1); % generate a random number
    event=sum(total_scaled_rate-ran_numb<0)+1; % find which event is happening

    switch(event)
        case 1
            %% Binding
            %% Selection random particle
            ran_numb_2=ran_numb;
            ran_numb=ran_numb/total_scaled_rate(1); % rescale rate
            if isnan(ran_numb)
                ran_numb_2;
            end

            counter_new_binding=counter_new_binding+1;

            left = sum(k_on_whole_particles_cumsum < (ran_numb * k_on_whole_particles_cumsum(end)))+1;
            % first and last index in function of left: 1 0 1; 2 1 interval; 3
            % interval + 1 interva.*2... if you drop the first case, (left-1)+1:(left)
            % at the end you have 
            start_id=(left-1)*interval_particle+1;
%             end_id=min((left)*interval_particle,Number_of_particles);
            if left==1
                previous=0;
            else
                previous=k_on_whole_particles_cumsum(left-1);
            end
            ran_numb=(ran_numb * k_on_whole_particles_cumsum(end)-previous)/(k_on_whole_particles_cumsum(left)-previous);
             if isnan(ran_numb)
                ran_numb_2;
            end
           
            t1=cumsum(array_slices_k_on_particles(:,left));
            idx=sum(t1 < ran_numb*t1(end))+1;
            % Select the corresponding id
            if idx==1
                ran_numb=ran_numb./(t1(idx)/t1(end));
            else
                ran_numb=(ran_numb*t1(end)-t1(idx-1))./(t1(idx)-t1(idx-1));
            end
            if isnan(ran_numb)
                ran_numb_2;
            end
       
            id_to_extract = start_id+idx-1;
            
            to_remove=array_slices_k_on_particles(idx,left);
            array_slices_k_on_particles(idx,left)=0;
            k_on_whole_particles_cumsum(left:end) = k_on_whole_particles_cumsum(left:end) - to_remove;
            
            % valid cause all particles have same probabilities. in case you need to
            % rearrange k_ons...
            
            bounding_particle=id_to_extract;
            particle_status(bounding_particle,1)=0;
            to_add=particle_status(bounding_particle,2);
            array_slices_k_off_particles(idx,left)=to_add;
            k_off_whole_particles_cumsum(left:end) = k_off_whole_particles_cumsum(left:end) + to_add;
            sum_kD=sum_kD+to_add;

            %% extract a free receptor 
           
            counter_new_binding_receptors=counter_new_binding_receptors+1;

            left = sum(k_on_whole_receptors_cumsum < (ran_numb * k_on_whole_receptors_cumsum(end)))+1;

            start_id=(left-1)*interval_receptors+1;
%             end_id=min((left)*interval_receptor,Number_of_receptors);
            if left==1
                previous=0;
            else
                previous=k_on_whole_receptors_cumsum(left-1);
            end
            ran_numb=(ran_numb * k_on_whole_receptors_cumsum(end)-previous)/(k_on_whole_receptors_cumsum(left)-previous);
            if isnan(ran_numb)
                ran_numb_2;
            end
            
            t1=cumsum(array_slices_unbound_receptors(:,left));
            idx=sum(t1 < ran_numb*t1(end))+1;

            id_to_extract = start_id+idx-1;
            bounding_receptor=id_to_extract;
            if Receptor(bounding_receptor,1)~=0
                disp('Selecting an already occupied receptor!')
                save('Eroror_log.mat')
                exit
            end
   
            % Remove receptor to be boundable
            to_remove=array_slices_unbound_receptors(idx,left);
            array_slices_unbound_receptors(idx,left)=0;
            k_on_whole_receptors_cumsum(left:end) = k_on_whole_receptors_cumsum(left:end) - to_remove;

            %% Execute the binding

            
            Receptor(bounding_receptor,1)=bounding_particle;
            particle_status(bounding_particle,3)=bounding_receptor;


            %% Update overall quantities
            free_receptors=free_receptors-1;
            free_particles=free_particles-1;
            previous_move=1;

        case 2
            %% Unbounding
            ran_numb=(ran_numb-total_scaled_rate(1))/(total_scaled_rate(2)-total_scaled_rate(1)); % here is different, multiplicities are important, 
            % extract a free particle and subunit
            counter_new_binding=counter_new_binding-1;

            left = sum(k_off_whole_particles_cumsum < (ran_numb * k_off_whole_particles_cumsum(end)))+1;
            % first and last index in function of left: 1 0 1; 2 1 interval; 3
            % interval + 1 interva.*2... if you drop the first case, (left-1)+1:(left)
            % at the end you have 
            start_id=(left-1)*interval_particle+1;
%             end_id=min((left)*interval_particle,Number_of_particles);
            if left==1
                previous=0;
            else
                previous=k_off_whole_particles_cumsum(left-1);
            end
            ran_numb=(ran_numb * k_off_whole_particles_cumsum(end)-previous)/(k_off_whole_particles_cumsum(left)-previous);
            
            t1=cumsum(array_slices_k_off_particles(:,left));
            idx=sum(t1 < ran_numb*t1(end))+1;
            % Select the corresponding id
            if idx==1
                ran_numb=ran_numb./(t1(idx)/t1(end));
            else
                ran_numb=(ran_numb*t1(end)-t1(idx-1))./(t1(idx)-t1(idx-1));
            end
            if isnan(ran_numb)
                ran_numb_2;
            end
            
            id_to_extract = start_id+idx-1;
            unbounding_particle=id_to_extract;

            to_remove=particle_status(unbounding_particle,2);

            k_off_whole_particles_cumsum(left:end) = k_off_whole_particles_cumsum(left:end) - to_remove;

            array_slices_k_off_particles(idx,left)=0;

            % valid cause all particles have same probabilities. in case you need to
            % rearrange k_ons, I recommend to not modify particle_status
            
            particle_status(unbounding_particle,1)=1;
            to_add=particle_status(unbounding_particle,1);
            array_slices_k_on_particles(idx,left)=to_add;
            k_on_whole_particles_cumsum(left:end) = k_on_whole_particles_cumsum(left:end) + to_add;
            sum_kD=sum_kD-particle_status(unbounding_particle,2);

            %% perform unbounding
            unbounding_receptor=particle_status(unbounding_particle,3);
            particle_status(unbounding_particle,3)=0;
            Receptor(unbounding_receptor)=0;
            free_receptors=free_receptors+1;
            free_particles=free_particles+1;

            %% Restore receptor 
            counter_new_binding_receptors=counter_new_binding_receptors-1;
            idx=mod(unbounding_receptor,interval_receptors);
            if(idx==0)
                idx=interval_receptors;
            end
            left=round((unbounding_receptor-idx)/interval_receptors)+1;
            to_add=1;
            array_slices_unbound_receptors(idx,left)=to_add;
            k_on_whole_receptors_cumsum(left:end) = k_on_whole_receptors_cumsum(left:end) + to_add;
            previous_move=2;
        otherwise
            disp('Error case not recognized')
            flag_ok=false(1);
            return
    end
    counter_events(event)=counter_events(event)+1;
    if mod(sum(counter_events),5e4)==0
            k_on_whole_particles_cumsum=cumsum(particle_status(:,1)); % Initial vector containing the cumulative rates of association
            k_on_whole_particles_cumsum=[k_on_whole_particles_cumsum(interval_particle:interval_particle:Number_of_particles);k_on_whole_particles_cumsum(end)]; % Resample the vector of cumulative rates of association
            k_off_whole_particles_cumsum=[cumsum(sum(array_slices_k_off_particles)');sum_kD];
    end
    if current_step==0
        if mod(current_sweep,frequency_of_output)==0

            id_to_consider=Receptor(Receptor(:)>0);
            surviving_kDs=log10(particle(id_to_consider,:,2));
            surviving_kDs=surviving_kDs(:);
            surviving_len=particle(id_to_consider,:,3);
            surviving_len=surviving_len(:);
            surviving_fam=particle(id_to_consider,:,1);
            surviving_fam=surviving_fam(:);
            raw=round(current_sweep/frequency_of_output);
            total_particle_bound(raw)=counter_new_binding;
            total_receptor_bound(raw)=counter_new_binding_receptors;
            average_kD_bound(raw)=mean(particle(id_to_consider,:,2),"all");
            for pp=1:Receptor_targettable_domains
                mean_kD_per_fam=mean(surviving_kDs(surviving_fam==pp),"all");
                mean_kD_per_fam(isnan(mean_kD_per_fam))=0;
                average_kD_bound_per_family(raw,pp)=mean_kD_per_fam;
                ligand_bound_per_family(raw,pp)=mean(surviving_len(surviving_fam==pp));
                ligand_bound_per_family(isnan(ligand_bound_per_family))=0;
            end
            save([simulation_cycle_end_name,'_provv.mat']); % checkpoints
        end
    end
end