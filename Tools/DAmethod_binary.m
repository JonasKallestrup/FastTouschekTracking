function [x_da,xp_da] = DAmethod_binary(dpoffset,Settings)
%{
Method to compute dynamic aperture using a binary search 

REQUIRED fields of "Settings" struct: 
Settings.ring - lattice
Settings.nTurns - number of turns to track
Settings.nlines - number of radial lines
Settings.xmax - maximum value of x
Settings.resolution_x - resolution of x of binary search. 

OPTIONAL fields of "Settings" struct:
Settings.orbit - orbit to do binary search around
Settings.G - matrix to normalize x,x' phase space 
Settings.interpDASteps - number or steps to interpolate the found DA
%}


ring = Settings.ring;
nTurns = Settings.nTurns;
nlines = Settings.nlines;
xmax = Settings.xmax;
resolution_x = Settings.resolution_x;

if isfield(Settings,'orbit')
    orbit_onenergy = Settings.orbit;
else
    if atGetRingProperties(ring,'is_6d')
        orbit_onenergy = findorbit6(ring,1);
    else
        orbit_onenergy = [findorbit4(ring,0,1);0;0];
    end
end
orbit_offenergy = findorbit4(ring,dpoffset,1);% estimate 4D closed orbit for off-momentum particle
orbit = orbit_onenergy+[orbit_offenergy;dpoffset;0];% 

theta = linspace(0,2*pi,nlines);
theta0 = theta;
r00 = repmat(orbit,1,nlines)+[1e-9;0;1e-9;0;0;0];

% convert xmax and x-resolution into an maximum ampltidue and resolution
if isfield(Settings,'G')
    G = Settings.G;
else
    [~,LinData] = atlinopt2(ring,1);
    G = [1/sqrt(LinData.beta(1)),0;LinData.alpha(1)/sqrt(LinData.beta(1)),sqrt(LinData.beta(1))];
end
resolution_A = G*[resolution_x;0];
resolution_A = sqrt(sum(resolution_A.^2));
Amax = G*[xmax;0];
Amax = sqrt(sum(Amax.^2));
Ginv = inv(G);

A_check = Amax*ones(1,nlines);
A_lim_surv = zeros(1,nlines);
A_lim_lost = Amax*ones(1,nlines);
A_da = nan(1,nlines);

n=0;
while ~isempty(theta)
    n=n+1;
    xvec_check = Ginv*[A_check.*cos(theta);A_check.*sin(theta)];
    r0 = r00;
    r0(1,:) = r0(1,:)+xvec_check(1,:);
    r0(2,:) = r0(2,:)+xvec_check(2,:);

    [~,lost] = ringpass(ring,r0,nTurns);

    toDelete = [];
    for a = 1:numel(lost)
        if ~lost(a) % particle survived
            if n == 1
                dA(a) = 0;
            else
                dA(a) = (A_lim_lost(a) - A_check(a))/2;% change in amplitude for next run
                A_lim_surv(a) = A_check(a);
            end
        else % particle is lost
            dA(a) = (A_lim_surv(a) - A_check(a))/2;% change in amplitude for next run
            A_lim_lost(a) = A_check(a);
        end
        if abs(dA(a)) < resolution_A
            if A_lim_surv(a) == 0
                A_da(a) = NaN;
            else
               A_da(a) = A_lim_surv(a); 
            end
            toDelete = [toDelete,a];
        else
            A_check(a) = A_check(a) + dA(a);
        end

    end

    % remove lines where boundary has already been discovered
    theta(toDelete) = [];
    A_lim_surv(toDelete) = [];
    A_lim_lost(toDelete) = [];
    A_check(toDelete) = [];
    r00(:,toDelete) = [];
    
end
if isfield(Settings,'interpDASteps') % do cubic interpolation of DA if requested
    if Settings.interpDASteps > 0 
        theta_interp = linspace(0,2*pi,Settings.interpDASteps);
        A_da = interp1(theta0,A_da,theta_interp,'cubic');
        theta0 = theta_interp;
    end
    
end
xvec_da = Ginv*[A_da.*cos(theta0);A_da.*sin(theta0)];
x_da = xvec_da(1,:)+orbit(1)+1e-9;
xp_da = xvec_da(2,:)+orbit(2);

