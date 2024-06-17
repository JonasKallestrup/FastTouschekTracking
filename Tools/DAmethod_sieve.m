function [x_da,xp_da] = DAmethod_sieve(dpoffset,Settings)
%{
Method to compute dynamic aperture using a seive search 
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
Settings.resolution_x_coarse -  coarse resolution of x search
Settings.resolution_x_fine - fine resolution of x search

OPTIONAL fields of "Settings" struct:
Settings.G - matrix to normalize x,x' phase space 
Settings.interpDASteps - number or steps to interpolate the found DA
%}


ring = Settings.ring;
nTurns = Settings.nTurns;
nlines = Settings.nlines;
xmax = Settings.xmax;
resolution_x_coarse = Settings.resolution_x_coarse;
resolution_x_fine = Settings.resolution_x_fine;

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
resolution_A_coarse = G*[resolution_x_coarse;0];
resolution_A_fine = G*[resolution_x_fine;0];
resolution_A_coarse = sqrt(sum(resolution_A_coarse.^2));
resolution_A_fine = sqrt(sum(resolution_A_fine.^2));

Amax = G*[xmax;0];
Amax = sqrt(sum(Amax.^2));
Ginv = inv(G);
Avalues_coarse = 0:resolution_A_coarse:Amax;


linesBeingChecked = 1:nlines;
k = 1;
while ~isempty(linesBeingChecked) && k <= numel(Avalues_coarse)
    xvec_check = Ginv*[Avalues_coarse(k).*cos(theta(linesBeingChecked));Avalues_coarse(k).*sin(theta(linesBeingChecked))];
    
    r0 = repmat(orbit,1,numel(linesBeingChecked))+[1e-9;0;1e-9;0;0;0];
    r0(1,:) = r0(1,:)+xvec_check(1,:);
    r0(2,:) = r0(2,:)+xvec_check(2,:);
    [~,lost] = ringpass(ring,r0,nTurns);
   
    k_da(linesBeingChecked(~lost)) = k;%update DA
    k=k+1;
    % remove lost lines
    linesBeingChecked(lost) = [];
end



Avalues_low = Avalues_coarse(k_da);
Avalues_high = Avalues_coarse(k_da+1);
for a =1:nlines
    Avalues_fine(:,a) = Avalues_low(a):resolution_A_fine:Avalues_high(a);
end

A_da = Avalues_low;
linesBeingChecked = 1:nlines;

k = 1;
while ~isempty(linesBeingChecked) && k <= numel(Avalues_fine)
    xvec_check = Ginv*[Avalues_fine(k,:).*cos(theta(linesBeingChecked));Avalues_fine(k,:).*sin(theta(linesBeingChecked))];
    
    r0 = repmat(orbit,1,numel(linesBeingChecked))+[1e-9;0;1e-9;0;0;0];
    r0(1,:) = r0(1,:)+xvec_check(1,:);
    r0(2,:) = r0(2,:)+xvec_check(2,:);
    [~,lost] = ringpass(ring,r0,nTurns);
   
    A_da(linesBeingChecked(~lost)) = Avalues_fine(k,~lost);%update DA
    k=k+1;
    % remove lost lines
    linesBeingChecked(lost) = [];
    Avalues_fine(:,lost) = [];
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
