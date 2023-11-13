function [x_da,xp_da] = DAmethod_dart(dpoffset,Settings)
%{
Method to compute dynamic aperture using dart algorithm. 

REQUIRED fields of "Settings" struct:
Settings.ring - lattice
Settings.nTurns - number of turns to track
Settings.xmax - maximum x coordinate
Settings.resolution_x - resolution of x-coordinate


OPTIONAL fields of "Settings" struct: 
Settings.orbit - orbit to do dart method around
Settings.G - matrix to normalize x,x' phase space 

%}

ring = Settings.ring;
nTurns = Settings.nTurns;
xmax = Settings.xmax;
resolution_x = Settings.resolution_x;

if isfield(Settings,'orbit')
    orbit_onEnergy = Settings.orbit;
else
    if atGetRingProperties(ring,'is_6d')
        orbit_onEnergy = findorbit6(ring,1);
    else
        orbit_onEnergy = [findorbit4(ring,0,1);0;0];
    end
end
% ring_offEnergy = atsetcavity(Settings.ring,'Frequency','nominal','dp',dpoffset);
% orbit_offEnergy = findorbit6(ring_offEnergy,1);
orbit_offEnergy = findorbit4(ring,dpoffset,1);

if isfield(Settings,'G')
    G = Settings.G;
else
    [~,LinData] = atlinopt2(ring,1);
    G = [1/sqrt(LinData.beta(1)),0;LinData.alpha(1)/sqrt(LinData.beta(1)),sqrt(LinData.beta(1))];
end


resolution_A = G*[resolution_x;0];resolution_A = sqrt(sum(resolution_A.^2));
Amax = G*[xmax;0];Amax = sqrt(sum(Amax.^2));
Ginv = inv(G);

A_check = Amax;
A_lim_surv = 0;
A_lim_lost = Amax;

boundaryFound = 0;

n=0;
while ~boundaryFound
    n=n+1;
    xvec = Ginv*[A_check;0];
    Rin = orbit_onEnergy+[orbit_offEnergy(1:4);0;0];
    Rin(1:2) = Rin(1:2)+xvec;
    Rin(5) = Rin(5)+dpoffset;
    Rin = Rin+[1e-9;0;1e-9;0;0;0];
    
    [~,lost] = ringpass(ring,Rin,nTurns);
    
    if ~lost% particle survived
        if n == 0
            dA = 0;
        else
            dA = (A_lim_lost - A_check)/2;% change in amplitude for next run
            A_lim_surv = A_check;
            Rin_surv = Rin;
        end
    else  % particle is lost
        dA = (A_lim_surv - A_check)/2;% change in amplitude for next run
        A_lim_lost = A_check;
    end
    if abs(dA) < resolution_A
        boundaryFound = 1;
        if A_lim_surv == 0
            x_da = NaN;
            xp_da = NaN;
        else
            % maximum amplitude found. Now find equivalent positions for constant dp
            [Rout_da,~] = ringpass(atradoff(ring),Rin_surv,Settings.nTurns);
            xout_surv = Rout_da(1,:);
            xpout_surv = Rout_da(2,:);
            [x_da,xp_da] = getBoundary(xout_surv(~isnan(xout_surv)),xpout_surv(~isnan(xout_surv)));
        end
    else
        A_check = A_check + dA;
    end
end

end
function [x,xp] = getBoundary(xout,xpout)
k = boundary(xout',xpout');
x = xout(k);
xp = xpout(k);
end