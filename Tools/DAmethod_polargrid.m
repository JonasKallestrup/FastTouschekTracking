function [x_da,xp_da] = DAmethod_polargrid(dpoffset,Settings)
%{
Method to compute dynamic aperture using a polar-grid search. 

REQUIRED fields of "Settings" struct: 
Settings.ring - lattice
Settings.nTurns - number of turns to track
Settings.nlines - number of lines to separate grid into
Settings.Jxmax - maximum horizontal action of grid
Settings.resolution_Jx - resolution of horizontal action of grid


OPTIONAL fields of "Settings" struct:
Settings.orbit - orbit to do polar grid search around. 
Settings.Fx - function that transforms J,theta to x
Settings.Fxp - function that transforms J,theta to x'

%}

ring = Settings.ring;
nTurns = Settings.nTurns;
nlines = Settings.nlines;
Jxmax = Settings.Jxmax;
resolution_Jx = Settings.resolution_Jx;

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

Jx = 0:resolution_Jx:Jxmax;
theta = linspace(0,2*pi,nlines);

[Jx,theta] = meshgrid(Jx,theta);
Jx = reshape(Jx,1,[]);
theta = reshape(theta,1,[]);

x_particle = Fx(Jx,theta);
xp_particle = Fxp(Jx,theta);

r0 = repmat(orbit,1,numel(x_particle));
r0(1,:) = r0(1,:) + x_particle;
r0(2,:) = r0(2,:) + xp_particle;
r0(5,:) = r0(5,:) + dpoffset;
r0 = r0 + [1e-9;0;1e-9;0;0;0];
[~,lost] = ringpass(ring,r0,nTurns);

if all(lost)
    x_da = nan;
    xp_da = nan;
else
    x_survived = x_particle(~lost);
    xp_survived = xp_particle(~lost);

    k = boundary(x_survived',xp_survived');
    x_da = x_survived(k);
    xp_da = xp_survived(k);
end


