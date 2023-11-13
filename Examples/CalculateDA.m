%{
Example script of how to use FastTouschekTracking to calculate momentum
aperture time-efficiently.

Computation of dynamic apertures (DA) is here done using binary search
%}
clear all;
clc;

Settings.ring = esrf;

%% Define parameters for DA computation. 
% there parameter is DA-method specific. 
Settings.nTurns = 2*512;
Settings.DAmethod = 'binary';
Settings.nlines = 101;
Settings.xmax = 50e-3;
Settings.resolution_x = 100e-6;

%%

Settings.nTurns = 512;
Settings.DAmethod = 'dart';
Settings.xmax = 30e-3;
Settings.resolution_x = 100e-6;


%% Define parameters for interface method
Settings.dpSteps_DA = linspace(-10e-2,10e-2,21);% dp steps for DA computation
Settings.dpResolution = 1e-3;% resolution of momentum aperture ocmputation


