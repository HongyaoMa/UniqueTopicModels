% Analyzing Autism Data IV --- Truncating short documents and write CSV
% files to be analyzed in the Python package
% The Autism NMF Project

% Hongyao Ma
% Created:  10/31/2014
% Modified: 10/31/2014

clear;
clc;
close all;

% Load the data with no all zero rows/columns
load('/Users/hma/asd_data_nz.mat')
% load('E:\My Files\Documents\data/asd_data_nz.mat');

% Size of the data matrix
[m, n] = size(data_nz);

% Parameters
D = m;                  % # of documents, which is # of patients
W = n;                  % # of words, which is # of codes
N = sum(sum(data_nz));  % Total number of words in the collection TODO


%% Construct the Word Count Data

data_out = [];

for i_patient = 1:m
    for i_code = 1:n
        if data_nz(i_patient, i_code) ~= 0
            data_out = [data_out; i_patient, i_code, data_nz(i_patient, i_code)];
        end
    end
end
        
%% Write the CSV Files

% Data Files
csvwrite('bag_of_codes.csv', data_out)
csvwrite('bag_of_codes.txt', data_out)
        
% Vocab File
csvwrite('codes_vocab.txt', (1:n)')   


%% Truncated

% Codes per Patient
wpd_cutoff = 10;
word_per_patient = sum(data_nz, 2);
code_per_Patient = sum(data_nz > 0, 2);

inds = word_per_patient >= wpd_cutoff;

data_nz = data_nz(inds, :);

% Patients per Code
dpw_cutoff = 10;
patient_per_code = sum(data_nz, 1);
data_nz = data_nz(:, patient_per_code >= dpw_cutoff);


% Size of the data matrix
[m, n] = size(data_nz);

% Parameters
D = m;                  % # of documents, which is # of patients
W = n;                  % # of words, which is # of codes
N = sum(sum(data_nz));  % Total number of words in the collection TODO

data_out = [];

for i_patient = 1:m
    for i_code = 1:n
        if data_nz(i_patient, i_code) ~= 0
            data_out = [data_out; i_patient, i_code, data_nz(i_patient, i_code)];
        end
    end
end
        
% Data Files
csvwrite('bag_of_codes_trunc_wpd10.csv', data_out)
csvwrite('bag_of_codes_trunc_wpd10.txt', data_out)
        
% Vocab File
csvwrite('codes_vocab_trunc_wpd10.txt', (1:n)')  
        

        
        