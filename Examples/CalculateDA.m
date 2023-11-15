%{
Example script of how to use FastTouschekTracking to calculate momentum
aperture time-efficiently.

Computation of dynamic apertures (DA) is here done using binary search
%}
clear all;
clc;
addpath('../Tools/')

% Settings.ring = atradon(atSetCavityPhase(esrf));
% Settings.outputName = 'FTTexample_ESRF';
% Settings.Ib = 200e-3/992;% bunch charge in mA

% Settings.ring = atradon(atSetCavityPhase(soleil));
% Settings.outputName = 'FTTexample_SOLEIL';
% Settings.Ib = 500e-3/312;% bunch charge in mA


% Settings.ring = atradon(atSetCavityPhase(australian_synchrotron));
% Settings.outputName = 'FTTexample_AusSync';
% Settings.Ib = 200e-3/330;% bunch charge in mA


%% Define parameters for DA computation. 
% these parameters are DA-method specific. 

Settings.nTurns = 512;% number of turns to track during DA computation
Settings.DAmethod = 'dart';% DA computation method
Settings.nlines = 101;% number of lines for theta in DA computation method
Settings.xmax = 50e-3;% maximum x-value of DA computation
Settings.resolution_x = 1000e-6;% resolution of DA computation




%% Define parameters for interface method and Touschek lifetime calc.
Settings.dpResolution = 1e-3;% resolution of momentum aperture ocmputation
Settings.dpmax = 15e-2;% maximum value of dp/p considered
Settings.dpSliceStep = 1e-2;% separation between x,x',dp/p slices


% calculate some beam parameters
[~,params] = atx(Settings.ring,0,1);
Settings.params = params;
Settings.params.modemittance(2) = 10e-12;%assume 10 pm vert. emittance

%% Now run computation
TL = FastTouschekTracking(Settings);

%% Do some plotting
load(Settings.outputName);

nDA = numel(x_da);
cmap0 = [215,48,39;244,109,67;253,174,97;254,224,144;255,255,191;224,243,248;171,217,233;116,173,209;69,117,180]/255;
step0 = linspace(0,1,size(cmap0,1));
step1 = linspace(0,1,nDA);

cmap(:,1) = interp1(step0,cmap0(:,1),step1);
cmap(:,2) = interp1(step0,cmap0(:,2),step1);
cmap(:,3) = interp1(step0,cmap0(:,3),step1);

figure('renderer','painters')
hold on;
for a = 1:numel(x_da)
    plot3(Settings.dpoffsets_DA(a)*ones(size(x_da{a}))*1e2,x_da{a}*1e3,xp_da{a}*1e3,'Color',cmap(a,:),'LineWidth',2)
end
view(3)
grid on
set(gca,'FontSize',14)
xlabel('dp/p [%]','FontSize',18)
ylabel('x [mm]','FontSize',18)
zlabel('x'' [mrad]','FontSize',18)



figure('renderer','painters')
spos = cat(1,LinData.SPos);
hold on;
plot(spos,MA_neg,'b')
plot(spos,MA_pos,'b')
grid on
set(gca,'FontSize',14)
xlabel('s [m]','FontSize',18)
ylabel('dp/p [%]','FontSize',18)


