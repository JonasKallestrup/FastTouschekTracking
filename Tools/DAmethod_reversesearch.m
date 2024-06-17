function [x_da,xp_da] = DAmethod_reversesearch(dpoffset,Settings)
%{
Method to compute dynamic aperture using a reverse search 
Technically, the search is performed by action, A, but the user specifies
the wished resolution in terms of "x", since horizontal offset is more
commonly used.

The code is optimized to track multiple lines simultaneously to reduce
tracking overhead.

REQUIRED fields of "Settings" struct: 
Settings.ring - lattice
Settings.nTurns - number of turns to track
Settings.nlines - number of radial lines
Settings.xmax - maximum value of x
Settings.resolution_x - resolution of x of binary search. 

OPTIONAL fields of "Settings" struct:
Settings.G - matrix to normalize x,x' phase space 
Settings.interpDASteps - number or steps to interpolate the found DA
%}


ring = Settings.ring;
nTurns = Settings.nTurns;
nlines = Settings.nlines;
xmax = Settings.xmax;
resolution_x = Settings.resolution_x;

orbit0 = findorbit(ring);
orbit_4D = findorbit4(atradoff(ring),'dp',dpoffset,'guess',[orbit0(1:4);dpoffset;0]);
orbit = findorbit(ring,'dp',dpoffset,'guess',[orbit_4D(1:4);orbit0(5)+dpoffset;orbit0(6)]);
if ~isOrbitClosedEnough(orbit,orbit0,dpoffset)
    if isnan(orbit_4D(1))
        x_da = NaN;
        xp_da = NaN;
        return
    else
        orbit = [orbit_4D(1:4);orbit0(5)+dpoffset;orbit0(6)];
    end
end

theta = linspace(0,2*pi,nlines);

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
Avalues = 0:resolution_A:Amax;

linesBeingChecked = 1:nlines;
A_da = zeros(1,nlines);
k = numel(Avalues);
while ~isempty(linesBeingChecked) && k >0 
    xvec_check = Ginv*[Avalues(k).*cos(theta(linesBeingChecked));Avalues(k).*sin(theta(linesBeingChecked))];
    
    r0 = repmat(orbit,1,numel(linesBeingChecked))+[1e-9;0;1e-9;0;0;0];
    r0(1,:) = r0(1,:)+xvec_check(1,:);
    r0(2,:) = r0(2,:)+xvec_check(2,:);
    [~,lost] = ringpass(ring,r0,nTurns);
   
    A_da(linesBeingChecked(~lost)) = Avalues(k);%update DA
    k=k-1;
    % remove lost lines
    linesBeingChecked(~lost) = [];
end

A_da(A_da==0) = NaN;
if isfield(Settings,'interpDASteps') % do cubic interpolation of DA if requested
    if Settings.interpDASteps > 0
        theta_interp = linspace(0,2*pi,Settings.interpDASteps);
        A_da = interp1(theta,A_da,theta_interp,'cubic');
        theta = theta_interp;
    end

end
xvec_da = Ginv*[A_da.*cos(theta);A_da.*sin(theta)];
x_da = xvec_da(1,:)+orbit(1)+1e-9;
xp_da = xvec_da(2,:)+orbit(2);


function yesno = isOrbitClosedEnough(orbit,orbit0,dpoffset)
dpoffset_reached = orbit(5)-orbit0(5);
yesno = ~(or(isnan(orbit(1)),abs(dpoffset) > 0 && abs((dpoffset-dpoffset_reached)/dpoffset) > 1e-2));
