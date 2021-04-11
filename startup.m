%% Add Assignment Folders in Path

addpath('Assignment1/');
addpath('Assignment1/plots/');
addpath('Assignment2/');
addpath('Assignment3/');
addpath('Assignment4/');
addpath('Assignment5/');
addpath('util/');
addpath('data/');
addpath('data/wind-dataset/');

%% Set Default Graphics Settings

set(groot,  'DefaultLineLineWidth', 1.2, ...
            'DefaultTextInterpreter', 'LaTeX', ...
            'DefaultAxesTickLabelInterpreter', 'LaTeX', ...
            'DefaultAxesFontName', 'LaTeX', ...
            'DefaultLegendInterpreter', 'LaTeX', ...
            'DefaultAxesLineWidth', 1.5, ...
            'DefaultStemLineWidth', 1.5, ...
            'DefaultStemMarker', 'none', ...
            'DefaultAxesTitleFontWeight','bold',...
            'DefaultAxesFontSize', 15, ...
            'DefaultAxesBox', 'on', ...
            'DefaultAxesColor', [1, 1, 1], ...
            'DefaultFigureColor', [1, 1, 1], ...
            'DefaultFigureColormap', parula(10)) 
        
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [0 44 22 9]);

% set(groot,'defaultfigureposition',[0 250 900 750])

%% Environment State
co = get(gca,'ColorOrder');
getcol = @(idx, alpha) [co(idx, :) alpha];
close all;

% set default seed
rng(0);