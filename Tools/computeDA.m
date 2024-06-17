function [u_da,up_da,x_da,xp_da,Settings] = computeDA(Settings)
%{
Function used to compute the dynamic apertures for various values of dp/p. 
A binary search finds the values of dp/p that gives the boundary for the
largest and smallest value of dp/p that still gives a non-zero dynamic
aperture. 
In between the two dp/p boundaries, a fixed dp/p resolution is used to
sample the x,x',dp/p 3d-volume.

REQUIRED fields of "Settings" struct: 
Settings.ring           - lattice
Settings.dpResolution   - boundary resolution of x,x',dp/p
Settings.dpmax          - maximum dp/p used for dynamic aperture search
Settings.dpSliceStep    - dp/p separation of two adjecent x,x' DA slices
Settings.G              - matrix to normalize x,x' phase space 
                            (inherited automatically from FastTouschekTracking)
%}


if strcmp(Settings.DAmethod,'line')
    DAmethod = @(dpoffset,Settings) DAmethod_line(dpoffset,Settings);
elseif strcmp(Settings.DAmethod,'sieve')
    DAmethod = @(dpoffset,Settings) DAmethod_sieve(dpoffset,Settings);
elseif strcmp(Settings.DAmethod,'binary')
    DAmethod = @(dpoffset,Settings) DAmethod_binary(dpoffset,Settings);
elseif strcmp(Settings.DAmethod,'reverse')
    DAmethod = @(dpoffset,Settings) DAmethod_reversesearch(dpoffset,Settings);
else
    error('DA method not recognized');
end


fprintf('      Computing DA for dp = 0\n')
[x_tmp,xp_tmp] = DAmethod(0,Settings);% Compute DA for dp=0
if all(isnan(x_tmp))% step failed
    error('Lattice has no dynamic aperture for dp = 0')
end
dpSteps_good(1) = 0;
x_da{1} = x_tmp;
xp_da{1} = xp_tmp;


% Do a binary search of x-x'-dp slices to find maximum/minimum relevant
% dp/p to within the dpResolution.
[dp_positive,x_da_positive,xp_da_positive] = binarySearch_DAdp(DAmethod,Settings.dpResolution, Settings.dpmax,Settings.dpResolution,Settings);
fprintf(['          Positive DA-DP boundary found: ',num2str(max(dp_positive)),'\n'])

[dp_negative,x_da_negative,xp_da_negative] = binarySearch_DAdp(DAmethod,-Settings.dpResolution, -Settings.dpmax,-Settings.dpResolution,Settings);
fprintf(['          Negative DA-DP boundary found: ',num2str(min(dp_negative)),'\n'])

dpSteps_good = [dpSteps_good,dp_positive,dp_negative];
x_da = [x_da,x_da_positive,x_da_negative];
xp_da = [xp_da,xp_da_positive,xp_da_negative];


%{
After boundaries have been found, we sample the DA-DP 4D space by computing
DAs for discrete values of dp. The values of dp to compute is decided by
the dpSteps_DA variable specified by the user. 
%}
dp_toCheck = min(dpSteps_good):Settings.dpSliceStep:max(dpSteps_good);
dp_toCheck = dp_toCheck(~ismember(dp_toCheck,dpSteps_good));
for a = 1:numel(dp_toCheck)
    fprintf(['      Computing DA for intermediate step. dp = ',num2str(dp_toCheck(a)),'\n'])
    [x_da_tmp,xp_da_tmp] = DAmethod(dp_toCheck(a),Settings);
    if ~any([isempty(x_da_tmp),any(isnan(x_da_tmp)),all(x_da_tmp==0)])
        dpSteps_good(end+1) = dp_toCheck(a);
        x_da{end+1} = x_da_tmp;
        xp_da{end+1} = xp_da_tmp; 
    end
end


%{
The dynamic aperture in real phase space is then converted to normalized
phase-space
%}
for a = 1:numel(x_da)
    if ~all(isnan(x_da{a}))
        utmp = Settings.G*[x_da{a};xp_da{a}];
        [u_da{a},up_da{a},theta_out] = discretisePolygon(utmp(1,:),utmp(2,:));
    else
        u_da{a} = nan;
        up_da{a} = nan;
    end
end

%{ 
Lastly, DA computations from both binary search and the equally-spaced dp steps
are combined and sorted according to their dp-value. 
%}
[Settings.dpoffsets_DA,I] = sort(dpSteps_good);
x_da = x_da(I);
xp_da = xp_da(I);
u_da = u_da(I);
up_da = up_da(I);



function [dp_stable,x_da_stable,xp_da_stable] = binarySearch_DAdp(DAmethod,value1, value2,resolution,Settings)
values = value1:resolution:value2;
low = 1;
high = numel(values);

dp_stable= [];
x_da_stable = {};
xp_da_stable = {};

while floor(high-low) > 1
    k = floor((high-low)/2)+low;
    fprintf(['      Computing DA for dp = ',num2str(values(k)),'\n'])
    [xout,xpout] = DAmethod(values(k),Settings);
    if any([isempty(xout),any(isnan(xout)),all(xout==0)])
        high = k;
    else
        low = k;
        dp_stable(end+1) = values(k);
        x_da_stable{end+1} = xout;
        xp_da_stable{end+1} = xpout;
    end
end






