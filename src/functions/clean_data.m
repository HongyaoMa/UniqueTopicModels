function data_new = clean_data(data_in, counts)
% clean_data takes away 1) words that are very infrequent and 2) documents
% that are too short or contain too few words

% The Autism NMF project
% Hongyao Ma
% Created:  11/25/2014
% Updated:  11/26/2014

% Sanity check
if length(counts) ~= 4
    error('Wrong count size');
end

% Cut the words first
wpw = sum(data_in, 1);
dpw = sum(data_in > 0, 1);
ind2 = logical( (wpw >= counts(3)) .* (dpw >= counts(4)) );
data_new = data_in(:, ind2);

% Then cut the documents
wpd = sum(data_in, 2);
cpd = sum(data_in > 0, 2);
ind1 = logical( (wpd >= counts(1)) .* (cpd >= counts(2)) );
data_new = data_new(ind1, :);