% MATLAB script to process M files, calculate distributions, and save to Excel
main_folder='./';
folder_reservoir=[{'sd_2_12'};{'sd_7_12'};{'sd_10_12'};{'sd_12_12'};{'sd_14_12'};{'sd_17_12'};{'sd_22_12'};{'sd_mono'}]';
repetitions=0:9;
replicas=length(repetitions);
% cycles_to_analyse = -1:6;  % Change this to the actual number of files you want to process
cycles_to_analyse = 6;
M=length(cycles_to_analyse);

% Predefine x-axis for particle(:,:,2) and particle(:,:,3)
deltakD=0.1;
x_axis_kD = -14:deltakD:-3; % For particle(:,:,2)
deltaLinker=0.1;
x_axis_linker = 0.5:deltaLinker:4.5; % For particle(:,:,3)
counter_fig=1;
max_subplot=1:length(folder_reservoir);
size_x=4;
size_y=2;
counter_subplot=0;
counter_fig_2=5;
counter_fig_3=9;
max_subplot_2=7;
size_x_2=4;
size_y_2=2;
for folder_name=folder_reservoir
    %% Initialization block
    counter_subplot=counter_subplot+1;
    replicas_per_cycle=zeros(1,M);
    
    % Initialize reservoir matrices for averages
    reservoir_matrix_kD_avg = zeros(length(x_axis_kD), M); % For particle(:,:,2)
    reservoir_matrix_status_avg = zeros(length(x_axis_kD), M); % For particle_status(:,2)
    reservoir_matrix_linker_avg = zeros(length(x_axis_linker), M); % For particle(:,:,3)
    reservoir_matrix_gm_linker_avg = zeros(length(x_axis_linker), M); % For particle(:,:,3)
    reservoir_bound_avg = zeros(1, M);
    % Initialize storage for individual replica data (2D arrays)
    replica_data_kD = zeros(length(x_axis_kD), M * replicas);
    replica_data_status = zeros(length(x_axis_kD), M * replicas);
    replica_data_linker = zeros(length(x_axis_linker), M * replicas);
    replica_data_gm_linker = zeros(length(x_axis_linker), M * replicas);
    replica_data_bound = zeros(M, replicas);
    replica_entropy_kD=zeros(M, replicas);
    replica_entropy_avidity=zeros(M, replicas);
    replica_entropy_linker=zeros(M, replicas);
    replica_entropy_gm_linker=zeros(M, replicas);

    %% Actual calculation of variables
    for repetition = repetitions
        folder =[main_folder,folder_name{1},'_',num2str(repetition)];
        file_suff='caroprese_vincenzo_06_09_2024_trial_1_cycle_';

        % Initialize reservoir matrices for the current replica
        reservoir_matrix_kD = zeros(length(x_axis_kD), M); % For particle(:,:,2)
        reservoir_matrix_status = zeros(length(x_axis_kD), M); % For particle_status(:,2)
        reservoir_matrix_linker = zeros(length(x_axis_linker), M); % For particle(:,:,3)
        reservoir_matrix_gm_linker = zeros(length(x_axis_linker), M); % For particle(:,:,3)
        reservoir_bound = zeros(1, M);
        
        for pos_index = 1:M
            counter_cycle=cycles_to_analyse(pos_index);
            if counter_cycle==-1
                filename = [folder, '/',file_suff, sprintf('%d_provv.mat', 0)];
            else
                filename = [folder, '/', file_suff, sprintf('%d.mat', counter_cycle)];
            end
            if exist(filename,"file")
                replicas_per_cycle(pos_index)=replicas_per_cycle(pos_index)+1;

                % Load the current .mat file
                load(filename,'particle','particle_status','counter_new_binding','Number_of_particles');

                % Distribution for particle(:,:,2) (kD values)
                 particle_kD = particle(:,:,2);
%                 particle_kD=min(log10(particle(:,:,2)),[],2);
                 particle_kD_flatten = log10(particle_kD(:));
%                 particle_kD_flatten=particle_kD(:);
                [counts_kD,entropy_kD]=pdf_and_entropy(particle_kD_flatten,x_axis_kD,deltakD);
                reservoir_matrix_kD(:, pos_index) = reservoir_matrix_kD(:, pos_index)+counts_kD';


                replica_entropy_kD(pos_index,repetition+1)=entropy_kD;
                % Store individual replica data (flattened to 2D)
                replica_data_kD(:, (pos_index-1)*replicas + repetition + 1) = counts_kD';

                % Distribution for particle_status(:,2)
                avidity = log10(particle_status(:, 2));
                [counts_status,entropy_avidity]=pdf_and_entropy(avidity,x_axis_kD,deltakD);

                reservoir_matrix_status(:, pos_index) = reservoir_matrix_status(:, pos_index) + counts_status';
                replica_entropy_avidity(pos_index,repetition+1)=entropy_avidity;
                
                % Store individual replica data (flattened to 2D)
                replica_data_status(:, (pos_index-1)*replicas + repetition + 1) = counts_status';

                % Distribution for particle(:,:,3) (linker length)
                particle_linker = particle(:,:,3);
                particle_linker_flatten = particle_linker(:);
                [counts_linker,entropy_linker]=pdf_and_entropy(particle_linker,x_axis_linker,deltaLinker);

                reservoir_matrix_linker(:, pos_index) = reservoir_matrix_linker(:, pos_index)+ counts_linker';
                replica_entropy_linker(pos_index,repetition+1)=entropy_linker;

                % Store individual replica data (flattened to 2D)
                replica_data_linker(:, (pos_index-1)*replicas + repetition + 1) = counts_linker';

                particle_gm_linker = prod(particle(:,:,3),2).^(1/size(particle,2));
                [counts_gm_linker,entropy_gm_linker]=pdf_and_entropy(particle_gm_linker,x_axis_linker,deltaLinker);
                reservoir_matrix_gm_linker(:, pos_index) = reservoir_matrix_gm_linker(:, pos_index)+ counts_gm_linker';
                replica_entropy_gm_linker(pos_index,repetition+1)=entropy_gm_linker;
                
                % Store individual replica data (flattened to 2D)
                replica_data_gm_linker(:, (pos_index-1)*replicas + repetition + 1) = counts_gm_linker';

                % Calculate and store binding fraction
                if counter_cycle>-1
                    reservoir_bound(:,pos_index)=counter_new_binding/Number_of_particles;
                    replica_data_bound(pos_index, repetition + 1) = counter_new_binding/Number_of_particles;
                end
            end
        end
        figure(counter_fig_2)
        subplot(size_x_2,size_y_2,counter_subplot)
        hold on
        plot(x_axis_kD,counts_kD);
        xlim([-7,-3])
        title(['kD Ligand ',strrep(erase(folder_name{1},'sd_'),'_','-')])
        xlabel('kD [10^x M]')
        ylabel('Pdf')
        figure(counter_fig_2+1)
        subplot(size_x_2,size_y_2,counter_subplot)
        hold on
        plot(x_axis_kD,counts_status);
    xlim([-9,-4])
    title(['Avidity ',strrep(erase(folder_name{1},'sd_'),'_','-')])
    xlabel('kDp [10^x M]')
    ylabel('Pdf')
                figure(counter_fig_2+2)
        subplot(size_x_2,size_y_2,counter_subplot)
        hold on
        plot(x_axis_linker,counts_linker);

       xlim([1,5])
    title(['Linker Ligand ',strrep(erase(folder_name{1},'sd_'),'_','-')])
    xlabel('Length [nm]')
    ylabel('Pdf')

        figure(counter_fig_2+3)
        subplot(size_x_2,size_y_2,counter_subplot)
        hold on
        plot(x_axis_linker,counts_gm_linker);
    xlim([1,5])
    title(['Geometric Mean Linker ',strrep(erase(folder_name{1},'sd_'),'_','-')])
    xlabel('Length [nm]')
    ylabel('Pdf')
        % Accumulate to calculate average over replicas
        reservoir_matrix_kD_avg = reservoir_matrix_kD_avg + reservoir_matrix_kD;
        reservoir_matrix_status_avg = reservoir_matrix_status_avg + reservoir_matrix_status;
        reservoir_matrix_linker_avg = reservoir_matrix_linker_avg + reservoir_matrix_linker;
        reservoir_matrix_gm_linker_avg = reservoir_matrix_gm_linker_avg + reservoir_matrix_gm_linker;
        reservoir_bound_avg = reservoir_bound_avg + reservoir_bound;
    end
    figure(counter_fig_3)
    subplot(size_x_2,size_y_2,counter_subplot)
    hold on
    plot(cycles_to_analyse+1,replica_entropy_kD);
    title(['Entropy Kd ',strrep(erase(folder_name{1},'sd_'),'_','-')])
    xlabel('Cycle [#]')
    ylabel('Entropy')
    ylim([0,1.1])

    figure(counter_fig_3+1)
    subplot(size_x_2,size_y_2,counter_subplot)
    hold on
    plot(cycles_to_analyse+1,replica_entropy_avidity);
    title(['Entropy Avidity ',strrep(erase(folder_name{1},'sd_'),'_','-')])
    xlabel('Cycle [#]')
    ylabel('Entropy')
    ylim([0,1.1])
    figure(counter_fig_3+2)
    subplot(size_x_2,size_y_2,counter_subplot)
    hold on
    plot(cycles_to_analyse+1,replica_entropy_linker);
    title(['Entropy Linker ',strrep(erase(folder_name{1},'sd_'),'_','-')])
    xlabel('Cycle [#]')
    ylabel('Entropy')
    ylim([0,1.1])
    figure(counter_fig_3+3)
    subplot(size_x_2,size_y_2,counter_subplot)
    hold on
    plot(cycles_to_analyse+1,replica_entropy_gm_linker);
    title(['Entropy gm linker ',strrep(erase(folder_name{1},'sd_'),'_','-')])
    xlabel('Cycle [#]')
    ylabel('Entropy')
    ylim([0,1.1])
    % Calculate averages
    reservoir_matrix_kD_avg = reservoir_matrix_kD_avg./replicas_per_cycle;
    reservoir_matrix_status_avg = reservoir_matrix_status_avg./replicas_per_cycle;
    reservoir_matrix_linker_avg = reservoir_matrix_linker_avg./replicas_per_cycle;
    reservoir_matrix_gm_linker_avg = reservoir_matrix_gm_linker_avg./replicas_per_cycle;
    entropy_kD=zeros(1,M);
    entropy_avidity=zeros(1,M);
    entropy_linker=zeros(1,M);
    entropy_gm_linker=zeros(1,M);
    for pos_index = 1:M
    entropy_kD(pos_index)=entropy_avg(reservoir_matrix_kD_avg,pos_index,deltakD);
    figure(counter_fig)
    subplot(size_x,size_y,counter_subplot)
    hold on
    plot(x_axis_kD,reservoir_matrix_kD_avg(:,pos_index));
    xlim([-7,-3])
    title(['kD Ligand ',strrep(erase(folder_name{1},'sd_'),'_','-')])
    xlabel('kD [10^x M]')
    ylabel('Pdf')
    figure(counter_fig+1)
    subplot(size_x,size_y,counter_subplot)
    hold on
    plot(x_axis_kD,reservoir_matrix_status_avg(:,pos_index));
    xlim([-9,-4])
    title(['Avidity ',strrep(erase(folder_name{1},'sd_'),'_','-')])
    xlabel('kDp [10^x M]')
    ylabel('Pdf')
   
    entropy_avidity(pos_index)=entropy_avg(reservoir_matrix_status_avg,pos_index,deltakD);
    figure(counter_fig+2)
    subplot(size_x,size_y,counter_subplot)
    hold on
    plot(x_axis_linker,reservoir_matrix_linker_avg(:,pos_index));
    xlim([1,5])
    title(['Linker Ligand ',strrep(erase(folder_name{1},'sd_'),'_','-')])
    xlabel('Length [nm]')
    ylabel('Pdf')

    entropy_linker(pos_index)=entropy_avg(reservoir_matrix_linker_avg,pos_index,deltaLinker);
    figure(counter_fig+3)
    subplot(size_x,size_y,counter_subplot)
    hold on
    plot(x_axis_linker,reservoir_matrix_gm_linker_avg(:,pos_index));
    xlim([1,5])
    title(['Geometric Mean Linker ',strrep(erase(folder_name{1},'sd_'),'_','-')])
    xlabel('Length [nm]')
    ylabel('Pdf')

    entropy_gm_linker(pos_index)=entropy_avg(reservoir_matrix_gm_linker_avg,pos_index,deltaLinker);
   
    end
    reservoir_bound_avg = reservoir_bound_avg./replicas_per_cycle;
    
    % Define column titles for Excel, including replicas with a white column
    cycleTitles = arrayfun(@(x) sprintf('Cycle %d', x), 0:M-1, 'UniformOutput', false);
    replicaTitles = arrayfun(@(x) sprintf('Replica %d', x), 0:replicas-1, 'UniformOutput', false);
    
    % Create headers for each block of replicas under a cycle title
    replicaHeaders = [];
    for i = 0:M-1
        replicaHeaders = [replicaHeaders, repmat(replicaTitles, 1, 1)];  %#ok<AGROW>
    end
    
    % Use non-merged version of the top header row, repeating the cycle title for each replica
    topCycleHeaders = [{'x axis'}, cycleTitles, {''}];
    for i = 0:M-1
        for r = 1:replicas
            topCycleHeaders = [topCycleHeaders, {sprintf('Cycle %d', i)}]; %#ok<AGROW>
        end
    end
    for i = 1:M
        topCycleHeaders{i+1}='Average';
    end
    % Prepare bottom headers for replicas
    colTitles_kD = [{'kD [10^x M]'}, cycleTitles, {''}, replicaHeaders];
    colTitles_status = [{'Avidity [10^x M]'}, cycleTitles, {''}, replicaHeaders];
    colTitles_linker = [{'Linker length [nm]'}, cycleTitles, {''}, replicaHeaders];
    colTitles_gm_linker = [{'GM Linker length [nm]'}, cycleTitles, {''}, replicaHeaders];

    % Add the x-axis to the matrix and the empty column
    reservoir_matrix_kD_with_x = [x_axis_kD', reservoir_matrix_kD_avg, nan(size(x_axis_kD', 1), 1), replica_data_kD];
    reservoir_matrix_status_with_x = [x_axis_kD', reservoir_matrix_status_avg, nan(size(x_axis_kD', 1), 1), replica_data_status];
    reservoir_matrix_linker_with_x = [x_axis_linker', reservoir_matrix_linker_avg, nan(size(x_axis_linker', 1), 1), replica_data_linker];
    reservoir_matrix_gm_linker_with_x = [x_axis_linker', reservoir_matrix_gm_linker_avg, nan(size(x_axis_linker', 1), 1), replica_data_gm_linker];
    reservoir_matrix_bound_with_x = [(cycles_to_analyse)'+1, reservoir_bound_avg', nan(size(cycles_to_analyse', 1), 1), replica_data_bound];

    reservoir_entropy_kD=[(cycles_to_analyse)'+1,entropy_kD',nan(size(cycles_to_analyse', 1), 1),replica_entropy_kD];
    reservoir_entropy_avidity=[(cycles_to_analyse)'+1,entropy_avidity',nan(size(cycles_to_analyse', 1), 1),replica_entropy_avidity];
    reservoir_entropy_linker=[(cycles_to_analyse)'+1,entropy_linker',nan(size(cycles_to_analyse', 1), 1),replica_entropy_linker];
    reservoir_entropy_gm_linker=[(cycles_to_analyse)'+1,entropy_gm_linker',nan(size(cycles_to_analyse', 1), 1),replica_entropy_gm_linker];
    
    % Save data to Excel with repeated cycle columns and replicas as columns
    excelFilename = [folder,'/', 'partial_caroprese_vincenzo_24_11_2024_distribution_level_',folder_name{1},'.xlsx'];

    %% Compact writing down
    % move all the reservoirs and titles to generic cells containers
    topCycleHeadersLast =[{'x axis'}, {'Average'},{''}];
    for r=1:replicas
        topCycleHeadersLast=[topCycleHeadersLast,{sprintf('Replica %d', r)}];
    end
    colTitles_fraction_bound = [{'Cycle #'},{'Fraction particle bound'}, {''}];
    for r=1:replicas
        colTitles_fraction_bound=[colTitles_fraction_bound,{'Fraction particle bound'}];
    end
    colTitles_entropy = [{'Cycle #'},{'Entropy'}, {''}];
    for r=1:replicas
        colTitles_entropy=[colTitles_entropy,{'Entropy'}];
    end
    
    topCycleHeadersCollection=[{topCycleHeaders},{topCycleHeaders},{topCycleHeaders},{topCycleHeaders},{topCycleHeaders},{topCycleHeadersLast},{topCycleHeadersLast},{topCycleHeadersLast},{topCycleHeadersLast}];
    colTitlesCollection=[{colTitles_kD},{colTitles_status},{colTitles_linker},{colTitles_gm_linker},{colTitles_fraction_bound},{colTitles_entropy},{colTitles_entropy},{colTitles_entropy},{colTitles_entropy}];
    reserovirMatrixCollection=[{reservoir_matrix_kD_with_x},{reservoir_matrix_status_with_x},{reservoir_matrix_linker_with_x},{reservoir_matrix_gm_linker_with_x},{reservoir_matrix_bound_with_x},{reservoir_entropy_kD},{reservoir_entropy_avidity},{reservoir_entropy_linker},{reservoir_entropy_gm_linker}];
    nameSheetCollection=[{'Ligand kD'},{'Particle kD'},{'Ligand Linker'},{'Geometric mean Linker'},{'Fraction particle bound'},{'Entropy kD'},{'Entropy Particle kD'},{'Entropy Ligand Linker'},{'Entropy Geometric mean Linker'}];
    for counterColl=1:size(nameSheetCollection,2)
        write_distrib_to_excel(topCycleHeadersCollection{counterColl},colTitlesCollection{counterColl},reserovirMatrixCollection{counterColl},excelFilename,nameSheetCollection{counterColl});
    end
end

function [counts,entropy]=pdf_and_entropy(quantity,x_axis,deltax)
    [counts_for_entropy, edges_distribution] = histcounts(quantity, 'Normalization', 'pdf');

    [counts, ~] = histcounts(quantity, [x_axis-deltax/2,x_axis(end)+deltax/2], 'Normalization', 'pdf');
    dx=diff(edges_distribution);
    entropy=-sum(log(counts_for_entropy(counts_for_entropy>0)).*counts_for_entropy(counts_for_entropy>0).*dx(counts_for_entropy>0));
end

function entropy=entropy_avg(reservoir,pos,deltax)
    pos_array=reservoir(reservoir(:,pos)>0,pos);

    entropy=-sum(log(pos_array).*pos_array.*deltax);

end
function write_distrib_to_excel(topCycleHeaders,colTitles,reservoir,excelFilename,name_sheet)
    writecell(topCycleHeaders, excelFilename, 'Sheet', name_sheet, 'Range', 'A1');
    % Write bottom row with replica headers
    writecell(colTitles, excelFilename, 'Sheet', name_sheet, 'Range', 'A2');
    writematrix(reservoir, excelFilename, 'Sheet', name_sheet, 'Range', 'A3');
end