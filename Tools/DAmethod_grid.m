function [x_da,xp_da,x_stable,xp_stable,x_unstable,xp_unstable] = DAmethod_grid(dpoffset,Settings)
%{
Method to compute dynamic aperture using a grid search. 

REQUIRED fields of "Settings" struct:
Settings.ring - lattice
Settings.nTurns - number of turns to track
Settings.xmax - maximum value of x (search is from -x to + x)
Settings.xpmax - maximum value of x' (serach is from -x' to +x')
Settings.resolution_x - resolution of x of grid
Settings.resolution_xp - resolution of x' of grid

OPTIONAL fields of "Settings" struct: 
Settings.orbit - orbit to do grid search around

The output of x_stable, xp_stable, x_unstable and xp_unstable is optional

%}

ring = Settings.ring;
nTurns = Settings.nTurns;
xmax = Settings.xmax;
xpmax = Settings.xpmax;
resolution_x = Settings.resolution_x;
resolution_xp = Settings.resolution_xp;

if isfield(Settings,'orbit')
    orbit = Settings.orbit;
else
    if atGetRingProperties(ring,'is_6d')
        orbit = findorbit6(ring,1);
    else
        orbit = [findorbit4(ring,0,1);0;0];
    end
end


xrange = -xmax:resolution_x:xmax;
xprange = -xpmax:resolution_xp:xpmax;

[X,XP] = meshgrid(xrange,xprange);
X = reshape(X,1,[]);
XP = reshape(XP,1,[]);
Rin = repmat(orbit,1,numel(X));
Rin(1,:) = Rin(1,:) + X;
Rin(2,:) = Rin(2,:) + XP;
Rin(5,:) = Rin(5,:)+dpoffset;
Rin = Rin + [1e-9;0;1e-9;0;0;0];

% [~,lost] = ringpass(ring,Rin,nTurns);
[~,lost] = ringpass_parallel(ring,Rin,nTurns);
lost = logical(lost);
x_stable = X(~lost);
xp_stable = XP(~lost);
x_unstable = X(lost);
xp_unstable = XP(lost);

bndry = boundary(x_stable',xp_stable');
x_da = x_stable(bndry);
xp_da = xp_stable(bndry);
