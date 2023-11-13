function [u_da,up_da,x_da,xp_da,Settings] = computeDA(Settings)
% Select Dynamic Aperture method

if strcmp(Settings.DAmethod,'grid')
    DAmethod = @(dpoffset,Settings) DAmethod_grid(dpoffset,Settings)
elseif strcmp(Settings.DAmethod,'polargrid')
    DAmethod = @(dpoffset,Settings) DAmethod_polargrid(dpoffset,Settings)
elseif strcmp(Settings.DAmethod,'sieve')
    DAmethod = @(dpoffset,Settings) DAmethod_sieve(dpoffset,Settings);
elseif strcmp(Settings.DAmethod,'binary')
    DAmethod = @(dpoffset,Settings) DAmethod_binary(dpoffset,Settings);
elseif strcmp(Settings.DAmethod,'floodfill')
    DAmethod = @(dpoffset,Settings) DAmethod_floodfill(dpoffset,Settings);
elseif strcmp(Settings.DAmethod,'polarfloodfill')
    DAmethod = @(dpoffset,Settings) DAmethod_polarfloodfill(dpoffset,Settings);
elseif strcmp(Settings.DAmethod,'dart')
    DAmethod = @(dpoffset,Settings) DAmethod_dart(dpoffset,Settings);
else
    error('DA method not recognized');
end


%{
A Binary Search is used to find the maximum and minimum value of dp that
still has non-zero dynamic aperture. 
%}


k=1;
[x_tmp,xp_tmp] = DAmethod(0,Settings);% Compute DA for dp=0
if all(isnan(x_tmp))% step failed
    error('Lattice has no dynamic aperture for dp = 0')
end
dpSteps_good(k) = 0;
x_da{k} = x_tmp;
xp_da{k} = xp_tmp;
k=k+1;

dpStep_lost = [];
dpStep_cur = max(Settings.dpSteps_DA);

% First compute the positive boundary
boundaryFound= 0;
i = 0;
while ~boundaryFound
    i = i+1;
    [x_tmp,xp_tmp] = DAmethod(dpStep_cur,Settings);
    if any(isnan(x_tmp)) || any(x_tmp==0)% step failed
        dpStep_lost = [dpStep_lost,dpStep_cur];
        dpStep_cur = max(dpSteps_good)+(min(dpStep_lost)-max(dpSteps_good))/2;
    else
        dpSteps_good(k) = dpStep_cur;
        x_da{k} = x_tmp;
        xp_da{k} = xp_tmp;
        k=k+1;
        if i == 1
           break; %the maximum dp value was still OK.
        end
        dpStep_cur = max(dpSteps_good)+ (min(dpStep_lost)-max(dpSteps_good))/2;

    end
    boundaryFound = min(abs(dpStep_lost-dpStep_cur))<Settings.dpResolution;
end
disp(['Positive DA-DP boundary found: ',num2str(max(dpSteps_good))])

% Next, compute the negative boundary
boundaryFound= 0;
dpStep_lost = [];
i = 0;
dpStep_cur = min(Settings.dpSteps_DA); 
while ~boundaryFound
    i = i+1;
    [x_tmp,xp_tmp] = DAmethod(dpStep_cur,Settings);
    if any(isnan(x_tmp)) || any(x_tmp==0)% step failed
        dpStep_lost = [dpStep_lost,dpStep_cur];
        dpStep_cur = min(dpSteps_good)+(max(dpStep_lost)-min(dpSteps_good))/2;
    else
        dpSteps_good(k) = dpStep_cur;
        x_da{k} = x_tmp;
        xp_da{k} = xp_tmp;
        k=k+1;
        if i == 1
           break; % The minimum dp value was still OK.
        end
        dpStep_cur = min(dpSteps_good)+ (min(dpStep_lost)-min(dpSteps_good))/2;
        
    end
    boundaryFound = min(abs(dpStep_lost-dpStep_cur))<Settings.dpResolution;
end
disp(['Negative DA-DP boundary found: ',num2str(min(dpSteps_good))])



%{
After boundaries have been found, we sample the DA-DP 4D space by computing
DAs for discrete values of dp. The values of dp to compute is decided by
the dpSteps_DA variable specified by the user. 
%}

dp_toCheck = Settings.dpSteps_DA(all([Settings.dpSteps_DA<max(dpSteps_good); Settings.dpSteps_DA>min(dpSteps_good)]));
dp_toCheck = dp_toCheck(~ismember(dp_toCheck,dpSteps_good));
for a = 1:numel(dp_toCheck)
    disp(['Computing dpStep #',num2str(a),' out of ',num2str(numel(dp_toCheck)),'. dp = ',num2str(dp_toCheck(a))])
    [x_da{k},xp_da{k}] = DAmethod(dp_toCheck(a),Settings);
    if all(x_da{k}==0)
        x_da{k} = NaN;
        xp_da{k} = NaN;
    end
    k=k+1;
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
dpoffsets_DA = [dpSteps_good,dp_toCheck];
[Settings.dpoffsets_DA,I] = sort(dpoffsets_DA);
x_da = x_da(I);
xp_da = xp_da(I);
u_da = u_da(I);
up_da = up_da(I);





