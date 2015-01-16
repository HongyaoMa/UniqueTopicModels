% Analyzing Autism Data I --- Basics
% The Autism NMF Project

% Hongyao Ma
% Created:  10/20/2014
% Modified: 10/20/2014

clear;
clc;
close all;
load('E:\My Files\Documents\data\asd_data.mat')

% load('/Users/hma/asd_data.mat')


%% Basic Info

% Size & Rank
[n_patient, n_code] = size(data_set);
n_elements = n_patient * n_code;
% original_rank = rank(data_set); 
original_rank = 2310;

% Extreme Values
any_neg = any(any(data_set < 0))
max_element = max(max(data_set))
min_element = min(min(data_set))

% Histogram of Entries
hist_min = 1;
hist_max = 10;
ind_range = logical((data_set <= hist_max) .* (data_set >= hist_min));
figure;
hist((data_set(ind_range)));

% Non-zero elements
data_nonzero = data_set(data_set > 0);
n_data_nonzero = length(data_nonzero)
mean_data_nonzero = mean(data_nonzero)
median_data_nonzero = median(data_nonzero)


%% Sparsity of Rows

% Sparsity
data_equal_zero = (data_set == 0);
num_zeros = sum(data_equal_zero(:));
ratio_zeros = num_zeros / n_elements

% Distribution of sparsity for each row
ratio_zeros_row = sum(data_equal_zero, 2) / n_code;
figure;
hist(ratio_zeros_row);

% All zero rows
all_zero_row = all(data_equal_zero, 2);
num_all_zero_row = sum(all_zero_row)

% Number of non-zero elements
num_nonzero = n_code - sum(data_equal_zero, 2);
hist_min = 0;
hist_max = 19;
ind_range = logical((num_nonzero <= hist_max) .* (num_nonzero >= hist_min));
figure;
hist(num_nonzero(ind_range))

% Statistics
mean_num_nonzero = mean(num_nonzero)
median_num_nonzero = median(num_nonzero)

%% Sparsity of Columns

% Distribution of sparsity for each column
ratio_zeros_col = sum(data_equal_zero, 1) / n_patient;
figure;
hist(ratio_zeros_col(ratio_zeros_col<0.7));

% All zero columns
all_zero_col = all(data_equal_zero, 1);
num_all_zero_col = sum(all_zero_col)

% Number of non-zero elements
num_nonzero_in_col = n_patient - sum(data_equal_zero, 1);
hist_min = 0;
hist_max = 19;
ind_range = logical((num_nonzero_in_col <= hist_max) .* (num_nonzero_in_col >= hist_min));
figure;
hist(num_nonzero_in_col(ind_range))

% Statistics
mean_num_nonzero_col = mean(num_nonzero_in_col)
median_num_nonzero_col = median(num_nonzero_in_col)

%% Columns that are not sparse
ind_nonsparsecol = find(ratio_zeros_col < 0.7);
cols = data_set(:, ind_nonsparsecol);
figure;
n_show = n_patient;
n_show = 100;
stem(1:n_show, cols(1:n_show,1));
hold on;
plot(1:n_show, cols(1:n_show,2), 'k*');

either_nz = double(cols(:,1) > 0) + double(cols(:,2) > 0);
perc_covered = sum(either_nz > 0)/n_patient


%% Kick out all zero rows and cols

% Raw data
figure;
imshow(data_set, [0, 1]);

% Kick out all zero rows and columns
data_nz = data_set(logical(1 - all_zero_row), logical(1-all_zero_col));
figure;
imshow(data_nz, [0, 1]);
[n_patient_nz, n_code_nz] = size(data_nz);


% Singular Values
[U, S, V] = svd(data_nz);

% Save the data
% save('/Users/hma/asd_data_nz.mat', 'data_set', 'data_nz', 'U', 'S', 'V');
save('E:\My Files\Documents\data/asd_data_nz.mat', 'data_set', 'data_nz', 'U', 'S', 'V');