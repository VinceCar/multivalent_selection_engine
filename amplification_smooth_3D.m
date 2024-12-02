id_to_consider=Receptor(:)>0;
surviving_kDs=log(particle(Receptor(id_to_consider),:,2));
surviving_kDs=surviving_kDs(:);
surviving_len=particle(Receptor(id_to_consider),:,3);
surviving_len=surviving_len(:);
surviving_fam=particle(Receptor(id_to_consider),:,1);
surviving_fam=surviving_fam(:);
[hist_2d,kds_axis,len_axis]=histcounts2(surviving_kDs,surviving_len,'Normalization','count','XBinLimits',[min(surviving_kDs),max(surviving_kDs)],'YBinLimits',[min(surviving_len),max(surviving_len)],'BinWidth',[0.1,0.1]);
discrete_distribution=zeros(Receptor_targettable_domains,length(kds_axis)-1,length(len_axis)-1);
for ii=1:Receptor_targettable_domains
    id_to_consider_per_fam=surviving_fam==ii;
    if sum(id_to_consider_per_fam)>0
        discrete_distribution(ii,:,:)=histcounts2(surviving_kDs(id_to_consider_per_fam),surviving_len(id_to_consider_per_fam),kds_axis,len_axis);
    end
end
discrete_distribution=discrete_distribution/length(surviving_fam);

% Meshgrid of bin centers for X, Y, and Z (for trilinear interpolation)
kD_centers = (kds_axis(1:end-1) + kds_axis(2:end)) / 2;
len_centers = (len_axis(1:end-1) + len_axis(2:end)) / 2;
domains_axis = 1:Receptor_targettable_domains;
domain_edges=0.5:(Receptor_targettable_domains+0.5);
[X, Y, Z] = ndgrid(kD_centers, len_centers, domains_axis);

% Step 1: Flatten the 3D histogram and compute the cumulative distribution (CDF)
pdf_1D = discrete_distribution(:);  % Flatten the 3D histogram into a 1D vector
cdf_1D = cumsum(pdf_1D);  % Cumulative sum to create CDF

% Step 2: Generate random uniform samples between 0 and 1
random_uniform = rand(Number_of_particles*Valency, 1);

% Step 3: Find the indices where the random numbers map to in the CDF
[~, ~,bin_indices] = histcounts(random_uniform, [0; cdf_1D]);
% Step 4: Convert the 1D bin indices back to 3D indices
[domains_bin,kD_bin, len_bin] = ind2sub(size(discrete_distribution), bin_indices);

% Step 5: Use trilinear interpolation within the grid to get continuous samples
% Generate continuous (x, y, z) samples using trilinear interpolation
kD_samples = zeros(Valency*Number_of_particles, 1);
len_samples = zeros(Valency*Number_of_particles, 1);
domains_samples = zeros(Valency*Number_of_particles, 1);

for i = 1:Valency*Number_of_particles
    % Find the bin boundaries for the current sample
    kD_left = kds_axis(kD_bin(i));       % Left boundary of the bin in x-direction
    kD_right = kds_axis(kD_bin(i) + 1);  % Right boundary of the bin in x-direction
    len_left = len_axis(len_bin(i));       % Left boundary of the bin in y-direction
    len_right = len_axis(len_bin(i) + 1);  % Right boundary of the bin in y-direction
    
    % Generate random (x, y, z) sample within the bin using trilinear interpolation
    kD_samples(i) = kD_left+ (kD_right-kD_left)*rand();
    len_samples(i) = len_left+ (len_right-len_left)*rand();
    domains_samples(i) = domains_bin(i);

end

% Assign everything to Pariticles
matrix_target=reshape(domains_samples,[Number_of_particles,Valency]);
matrix_kD=reshape(exp(kD_samples),[Number_of_particles,Valency]);
matrix_len=reshape(len_samples,[Number_of_particles,Valency]);

[~, sort_idx] = sort(matrix_kD, 2,'descend');  % Sort along the 2nd dimension (columns)
lin_ids = ((sort_idx-1).*size(matrix_kD,1))+repmat([1:size(matrix_kD,1)]',1,Valency);
matrix_kD(:)=matrix_kD(lin_ids(:));
matrix_target(:)=matrix_target(lin_ids(:));
matrix_len(:)=matrix_len(lin_ids(:));
particle(:,:,1)=matrix_target;
particle(:,:,2)=matrix_kD;
particle(:,:,3)=matrix_len;

% computing particle status
generate_particle_level_properties