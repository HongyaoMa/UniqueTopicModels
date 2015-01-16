function x = gen_Docs(topics, A, l_Doc, length_type)
% gen_Doc generates documents with topic distribution A and the probability
% of topics "topics". The length of documents generated is l_Doc

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
end

% Word counts in each documenet
x = mnrnd(array_l_Doc, word_distr');


