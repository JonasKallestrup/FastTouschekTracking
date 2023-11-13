
function [x_da,xp_da] = DAmethod_sieve(dpoffset,Settings)
%{
Method to compute dynamic aperture using 'sieve'-search. 

REQUIRED fields of "Settings" struct:
Settings.ring - lattice
Settings.nTurns - number of turns to track
Settings.nlines - number of lines to track along
Settings.Jxmax - maximum action for search
Settings.resolution_Jx - resolution of action
Settings.coarseness - coarseness of first search (typically 10)

OPTIONAL fields of "Settings" struct: 
Settings.orbit - orbit to do grid search around
Settings.Fx - function that transforms J,theta to x
Settings.Fxp - function that transforms J,theta to x'

%}

ring = Settings.ring;
nTurns = Settings.nTurns;
nlines = Settings.nlines;
Jxmax = Settings.Jxmax;
resolution_Jx = Settings.resolution_Jx;
coarseness = Settings.coarseness;

theta = linspace(0,2*pi,nlines);

if isfield(Settings,'orbit')
    orbit = Settings.orbit;
else
    if atGetRingProperties(ring,'is_6d')
        orbit = findorbit6(ring,1);
    else
        orbit = [findorbit4(ring,0,1);0;0];
    end
end
if isfield(Settings,'Fx') && isfield(Settings,'Fxp')
    Fx = Settings.Fx;
    Fxp = Settings.Fxp;
else
    [~,LinData] = atlinopt2(ring,1);
    Fx = @(J,theta) sqrt(2*J*LinData.beta(1)).*cos(theta);
    Fxp = @(J,theta) sqrt(2*J/LinData.beta(1)).*(LinData.alpha(1)*cos(theta)+sin(theta));
end

dr = [1e-9;0;1e-9;0;0;0];
     
% perform very coarse grid search
for nn=1:numel(theta)
    Jx_coarse = 0:resolution_Jx*coarseness:Jxmax;
    x_coarse = Fx(Jx_coarse,theta(nn));
    xp_coarse = Fxp(Jx_coarse,theta(nn));
    npoints_coarse = numel(Jx_coarse);

    r0_coarse = zeros(6,npoints_coarse);
    r0_coarse(1,:) = x_coarse;
    r0_coarse(2,:) = xp_coarse;
    r0_coarse(5,:) = dpoffset;
    r0_coarse = r0_coarse+orbit+dr;

    r1_coarse=ringpass(ring,r0_coarse,nTurns);

    loss = find(isnan(r1_coarse(1,end-npoints_coarse+1:end)));
    if isempty(loss)
        x_da(nn) = x_coarse(end);
        xp_da(nn) = xp_coarse(end);
    elseif loss(1) == 1% all points are lost
        x_da(nn) = nan;
        xp_da(nn) = nan;
    else
        survived = loss(1)-1;% last point that survives
        if survived == 1
            Jx_fine= 0:resolution_Jx:Jx_coarse(2);
        else
            Jx_fine = Jx_coarse(survived-1):resolution_Jx:Jx_coarse(survived+1);
        end

        x_fine = Fx(Jx_fine,theta(nn));
        xp_fine = Fxp(Jx_fine,theta(nn));
        npoints_fine = numel(Jx_fine);
        r0_fine = zeros(6,npoints_fine);
        r0_fine(1,:) = x_fine;
        r0_fine(2,:) = xp_fine;
        r0_fine(5,:) = dpoffset;
        r0_fine = r0_fine+orbit+dr;
        r1_fine=ringpass(ring,r0_fine,nTurns);
        loss = find(isnan(r1_fine(1,end-npoints_fine+1:end)));
        if isempty(loss)
            x_da(nn) = x_fine(end);
            xp_da(nn) = xp_fine(end);
        elseif loss(1) == 1
            x_da(nn) = nan;
            xp_da(nn) = nan;
        else
            x_da(nn) = x_fine(loss(1)-1);
            xp_da(nn) = xp_fine(loss(1)-1);
        end
    end
end
end
