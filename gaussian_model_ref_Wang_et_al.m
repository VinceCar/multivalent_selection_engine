%% Generate the particles and receptors

%% first check if everything that is needed for this model is there
if ~exist('vector_of_average_kDs','var')
    fprintf("Error! You did not specify any average kD");
    flag_ok=false(1);
    return
end
if ~exist('vector_of_standard_deviations_kDs','var')
    fprintf("Error! You did not specify any deviation for the distribution of kD");
    flag_ok=false(1);    
    return
end

% Check for coherence
if length(vector_of_average_kDs)<Receptor_targettable_domains
    fprintf("Error! The vector of average kDs is not as long as the expected targettable domains");
    flag_ok=false(1);
    return
end
if length(vector_of_standard_deviations_kDs)<length(vector_of_average_kDs)
    fprintf("Error! the number of standard deviations provided is different than the one for the averages");
    flag_ok=false(1);
    return
end


%% Actual generation of kD variables
% Favourite target
matrix_target=Receptor_targettable_domains-reshape(sum(rand(numel(particle(:,:,1)),1)<cumsum(target_fraction_of_relative_population(:))',2),size(particle(:,:,1)))+1; % In this way we randomly generate the favourite units, taking into account quantization processes.

% kD to that target

% Now the kDs will be generated through random number generations sampling
% from a uniform distribution to do that, we need to compute the inverse
% CDF. In this case I am using the inverse of a Gaussian

invCDF=@(x,mu,sigma) mu+sigma.*sqrt(2).*erfinv(2.*x-1);
vector_of_average_kDs_col=vector_of_average_kDs';
vector_of_standard_deviations_col=vector_of_standard_deviations_kDs';
matrix_kD=exp(invCDF(rand(size(particle(:,:,2))),vector_of_average_kDs_col(matrix_target),vector_of_standard_deviations_col(matrix_target)));
[~, sort_idx] = sort(matrix_kD, 2,'descend');  % Sort along the 2nd dimension (columns)
lin_ids = ((sort_idx-1).*size(matrix_kD,1))+repmat([1:size(matrix_kD,1)]',1,Valency);
matrix_kD(:)=matrix_kD(lin_ids(:));
matrix_target(:)=matrix_target(lin_ids(:));
particle(:,:,1)=matrix_target;
particle(:,:,2)=matrix_kD;
%% first check if everything that is needed for this model is there
if ~exist('vector_of_average_len','var')
    fprintf("Error! You did not specify any average of linker lengths");
    flag_ok=false(1);
    return
end
if ~exist('vector_of_standard_deviations_len','var')
    fprintf("Error! You did not specify any deviation for the distribution of linker lengths");
    flag_ok=false(1);    
    return
end
if ~exist('vector_of_support_len','var')
    fprintf("Error! You did not specify the finite support for the distribution of linker lengths");
    flag_ok=false(1);    
    return
end

% Check for coherence
if length(vector_of_average_len)<Receptor_targettable_domains
    fprintf("Error! The vector of average linker lengths is not as long as the expected targettable domains");
    flag_ok=false(1);
    return
end
if length(vector_of_standard_deviations_len)<length(vector_of_average_len)
    fprintf("Error! the number of standard deviations provided is different than the one for the averages");
    flag_ok=false(1);
    return
end


%% Actual generation of kD variables

% Now the len will be generated through random number generations sampling
% from a uniform distribution to do that, we need to compute the inverse
% CDF. In this case I am using the inverse of a truncated Gaussian (you
% just limit the support)

% Compute the CDF values for the truncation bounds
cdf_lower = 0.5 * (1 + erf((-vector_of_support_len) ./ (vector_of_standard_deviations_len * sqrt(2))));
cdf_upper = 0.5 * (1 + erf((vector_of_support_len) ./ (vector_of_standard_deviations_len * sqrt(2))));
cdf_lower=cdf_lower(:);
cdf_upper=cdf_upper(:);
% Scale uniform random numbers to the truncated CDF range
random_numbers_to_sample=rand(size(particle(:,:,2))).*(cdf_upper(particle(:,:,1)) - cdf_lower(particle(:,:,1)))+cdf_lower(particle(:,:,1));

% Apply the inverse CDF to generate numbers from the truncated Gaussian
vector_of_average_len_col=vector_of_average_len';
vector_of_standard_deviations_len_col=vector_of_standard_deviations_len';

particle(:,:,3)=invCDF(random_numbers_to_sample,vector_of_average_len_col(particle(:,:,1)),vector_of_standard_deviations_len_col(particle(:,:,1)));
