%{
This is the main file of OpenRocket - Skyward Matlab version

NOTE: the purpose of this work is to show how much openrocket sucks with
respect to Missile Datcom <3
%}

clear all
close all
clc

%% PATH
filePath = fileparts(mfilename('fullpath'));
currentPath = pwd;
if not(strcmp(filePath, currentPath))
    cd (filePath);
    currentPath = filePath;
end

addpath(genpath(currentPath));

%% SETTINGS
settings.Mach = 0.1;
%%% NOSE
settings.Dpitot = [0.0025 0.006 0.006 0.025];
settings.Lpitot = [0.0048 0.0379 0.0268];
settings.Dnose = [0.025 0.15];
settings.Lnose = 0.2898;
settings.LnoseTrue = 0.31;
settings.Tnose = 'HAACK';
settings.bnose = 0;

%%% CENTERBODY
settings.Dcentr = 0.15;
settings.Lcentr = 2.096;

%%% AFTERBODY
settings.Daft = [0.15 0.122];
settings.Laft = 0.07;

%%% FINS
settings.Nfins = 3;
settings.Xle = 1.7440;
settings.xi = [0 0.21 0.31 0.34];
settings.yi = [0 0.15 0.15 0];

%%
settings.Aref = pi * (settings.Dcentr^2) /4;
settings = synths(settings);

%%

out = aerodynamics(settings, linspace(-13*pi/180, 13*pi/180, 27));

save PyxisOpenRocket out