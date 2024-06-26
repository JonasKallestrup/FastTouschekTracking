%{
Example script of how to use FastTouschekTracking to calculate momentum
aperture time-efficiently in the Soleil storage ring. 
Note: time-gain by FTT is most relevant for complex rings, e.g., 4th Gen
storage rings. 

Computation of dynamic apertures (DA) is here done using binary search
%}
clear all;
clc;
addpath('../Tools/')

Settings.Ib = 300e-3/416;% bunch charge in mA
Settings.ring = atSetRingProperties(atenable_6d(soleil));
Settings.outputName = 'FTTexample_soleil';


%% Define parameters for DA computation. 
% these parameters are DA-method specific. 

Settings.nTurns = 1024;% number of turns to track during DA computation
Settings.DAmethod = 'binary';% DA computation method
Settings.nlines = 25;% number of lines for theta in DA computation method
Settings.interpDASteps = 0;% number of interpolation steps for DA-reconstruction. Ensure interpDASteps > nlines
Settings.xmax = 100e-3;% maximum x-value of DA computation
Settings.resolution_x = 100e-6;% resolution of DA computation


%% Define parameters for interface method and Touschek lifetime calc.
Settings.dpResolution = 1e-3;% resolution of momentum aperture ocmputation
Settings.dpmax = 10e-2;% maximum value of dp/p considered
Settings.dpSliceStep = 1e-2;% separation between x,x',dp/p slices


% calculate some beam parameters needed for the Touschek Lifetime estimation
[~,params] = atx(Settings.ring,0,1);
Settings.params = params;
Settings.params.modemittance(2) = 10e-12;%assume 10 pm vertical emittance

%% Now run computation
[TL,MA_pos,MA_neg]= FastTouschekTracking(Settings);

%% Do some plotting
load('FTTexample_soleil')
nDA = numel(x_da);
cmap0 = [215,48,39;244,109,67;253,174,97;254,224,144;255,255,191;224,243,248;171,217,233;116,173,209;69,117,180]/255;
step0 = linspace(0,1,size(cmap0,1));
step1 = linspace(0,1,nDA);

cmap(:,1) = interp1(step0,cmap0(:,1),step1);
cmap(:,2) = interp1(step0,cmap0(:,2),step1);
cmap(:,3) = interp1(step0,cmap0(:,3),step1);

fig_polyhedron = figure('renderer','painters');
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



figure('renderer','painters');
spos = cat(1,LinData.SPos);
hold on;
p_FTT = plot(spos,MA_neg,'b');
plot(spos,MA_pos,'b')
grid on
set(gca,'FontSize',14)
xlabel('s [m]','FontSize',18)
ylabel('dp/p [%]','FontSize',18)

%% Calculate a few random locations in the ring using standard momentum acceptance calculation
clear all;
load('FTTexample_soleil')
ring = Settings.ring;
index = 1:length(ring);
length_elems = atgetfieldvalues(ring,index,'Length');
index(length_elems==0) = [];% don't track where there's no length

% pick a few random points to test to speed up computation
index = sort(randperm(numel(index),11));

nTurns = 1024;
dpResolution = 1e-3;
dpmax = 10e-2;
[ap_pos,ap_neg] = momentumAperture_binary(ring,index,nTurns,dpResolution,dpmax);
save('MAexample_soleil_binary','ap_pos','ap_neg','index')
%% now compare results
clear all;

load('FTTexample_soleil.mat')

figure('renderer','painters');
spos = cat(1,LinData.SPos);
hold on;
p_FTT = plot(spos,MA_neg,'b');
plot(spos,MA_pos,'b')
grid on
set(gca,'FontSize',14)
xlabel('s [m]','FontSize',18)
ylabel('dp/p [%]','FontSize',18)

load('MAexample_soleil_binary.mat')

spos_index = findspos(Settings.ring,index);
p_MA = plot(spos_index,ap_pos,'ro');
plot(spos_index,ap_neg,'ro')
legend([p_FTT,p_MA],'MA from FTT','MA from binary search','Location','east')
