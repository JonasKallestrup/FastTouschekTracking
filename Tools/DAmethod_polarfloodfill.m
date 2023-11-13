function [x_da,xp_da] = DAmethod_polarfloodfill(dpoffset,Settings)
%{
Method for computing dynamic aperture using polar floodfill

REQUIRED fields of "Settings" struct:
Settings.ring - lattice
Settings.nTurns - number of turns to track
Settings.Jxmax - maximum horizontal action of grid
Settings.resolution_Jx - resolution of horizontal action of grid
Settings.resolution_theta - resolution of theta

OPTIONAL fields of "Settings" struct:
Settings.orbit - orbit to do binary search around
Settings.Fx - function that transforms J,theta to x
Settings.Fxp - function that transforms J,theta to x'

%}


ring = Settings.ring;
nTurns = Settings.nTurns;
Jxmax = Settings.Jxmax;
resolution_Jx = Settings.resolution_Jx;
resolution_theta = Settings.resolution_theta;

Jx = 0:resolution_Jx:Jxmax;
theta = 0:resolution_theta:2*pi-resolution_theta;
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

% compress test values into Nx2 matrix used to see if computation is done. 
k = 1;
for a = 1:numel(Jx)
    for b = 1:numel(theta)
        testvalues(k,1:2) = [Jx(a),theta(b)];
        k=k+1;
    end
end


% Create queue and start at maximum Jx but at the middle-value of theta
% (i.e. to go from pi towards 2*pi and 0. Makes the calculation easier)
queue = {};
queue{1} = [max(Jx),theta(ceil(numel(theta)/2))];

dr = [1e-9;0;1e-9;0;0;0];
h=1;
nchecks=1;% to calculate how many evaluations are made
x_da = [];
xp_da = [];
while ~isempty(queue) && ~isempty(testvalues)

    nchecks = nchecks+1;

    nqueue = numel(queue);
    Jx_check = [];
    theta_check = [];
    x_particle = [];
    xp_particle = [];

    for b = 1:nqueue
        Jx_check(b) = queue{1}(1);
        if Jx_check(b) < 0
            Jx_check(b)=0;
        end
        theta_check(b) = queue{1}(2);
        x_particle(b) = Fx(Jx_check(b),theta_check(b));
        xp_particle(b) = Fxp(Jx_check(b),theta_check(b));
        queue(1) = [];
    end

    r0 = repmat(orbit,1,nqueue);
    r0(1,:) = r0(1,:)+x_particle;
    r0(2,:) = r0(2,:)+xp_particle;
    r0(5,:) = r0(5,:) + dpoffset;

    r0 = r0 +dr;
    try
    r1 = ringpass(ring,r0,nTurns);
    
    catch
        disp('heh')
    end

    loss = isnan(r1(1,end-nqueue+1:end));
    for b = 1:nqueue
        if ~loss(b)
            x_da(h) = x_particle(b);% add values these values of Jx and theta to the "surviving" points
            xp_da(h) = xp_particle(b);
            h=h+1;
        else
            combination1 = [Jx_check(b),theta_check(b)-resolution_theta];% check same Jx as currently, but smaller value of theta
            combination2 = [Jx_check(b),theta_check(b)+resolution_theta];% check same Jx as currently, but higher value of theta
            combination3 = [Jx_check(b)-resolution_Jx,theta_check(b)];% check smaller Jx, but same value of theta

            [validTest1,I1] = ismembertol(combination1,testvalues,'ByRows',1e-10);% check that the new combinations are still part of the testvalues (or have already been checked and therefore doesn't appear in testvalues)
            [validTest2,I2] = ismembertol(combination2,testvalues,'ByRows',1e-10);
            [validTest3,I3] = ismembertol(combination3,testvalues,'ByRows',1e-10);
            doDelete = [];
            if validTest1
                queue{end+1} = combination1;% add to queue
                doDelete = [doDelete,I1];
            end
            if validTest2
                queue{end+1} = combination2;
                doDelete = [doDelete,I2];
            end
            if validTest3
                queue{end+1} = combination3;
                doDelete = [doDelete,I3];
            end
            testvalues(doDelete,:) = [];% remove new combinations from testvalues; they are now in the queue
        end

    end
end

if isempty(x_da)
    x_da = nan;
    xp_da = nan;
else
    k = boundary(x_da',xp_da');
    x_da = x_da(k);
    xp_da = xp_da(k);
end


