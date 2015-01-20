function saveAllFig(foldername, varargin)
% SAVEALLFIG saves all currently showed figures in the folder 'foldername'
% INPUT:
%       foldername: The folder in which the figures are saved
%       varargin:   A list of formats to save to figures, e.g. 'fig'
%                   If no type is specified, save the fig file by default

% The Phase Retrieval Project
% Hongyao Ma
% 03/17/2014

mkdir(foldername);

% The number of types 
N = length(varargin);
if N == 0;
    varargin{1} = 'fig';
    N = 1;
end

% Save all the figures
h = get(0,'children');
for i=1:length(h)
    for i_type = 1:N
        saveas(h(i), [foldername, '/figure' num2str(i)], varargin{i_type});    
    end
end
