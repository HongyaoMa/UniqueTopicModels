function x = gen_Docs(topics, A, l_Doc, length_type)
% gen_Doc generates documents with topic distribution A and the probability
% of topics "topics". The length of documents generated is l_Doc
%
% INPUT
%   topics - n_doc by k matrix, distribution over topics
%   A - the topic distribution matrix
%   l_Doc - the expected length of documents
%   length_type - how to generate the lengths
%       'fixed' - every document has the same length
%       'uniform' - documents have lengths uniformly distributed [1, 2*l-1]
%
% OUTPUT
%   x - word count matrix, l_Doc by # of words
%

% The Autism NMF Project
% Hongyao Ma
% Created:  01/10/2015
% Modified: 01/10/2015

% Default length type
if nargin == 3
    length_type = 'fixed';
end

% Parameters
n_docs = size(topics, 1);

% Make sure that the sizes match
if size(topics, 2) ~= size(A, 2)
    error('Numbers of topics do not match!');
end

% Distribution over words
word_distr = A * topics';

% Length Type
switch length_type
    case 'fixed'
        array_l_Doc = l_Doc;
    case 'random'
        array_l_Doc = randi([2, 2*l_Doc-2], n_docs, 1); 
    otherwise
        error('Unknown type of the lengths!');
end

% Word counts in each documenet
x = mnrnd(array_l_Doc, word_distr');


