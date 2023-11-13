function [x_da,xp_da] = DAmethod_floodfill(dpoffset,Settings)
%{
Method to compute dynamic aperture using flood-fill algorithm.
Note: this method is not particularly efficient in MATLAB/AT due to
overhead of atpass method

REQUIRED fields of "Settings" struct: 
Settings.ring - lattice
Settings.nTurns - number of turns to track
Settings.xmax - maximum value of x (search is from -x to + x)
Settings.xpmax - maximum value of x' (serach is from -x' to +x')
Settings.resolution_x - resolution of x of floodfill-grid
Settings.resolution_xp - resolution of x' of floodfill-grid

OPTIONAL fields of "Settings" struct: 
Settings.orbit - orbit to do floodfill around


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
nStepsX  = numel(xrange);
nStepsXP  = numel(xprange);
M = zeros(nStepsXP,nStepsX);

queue(1,1:2) = [1,1];
k=1;
while ~isempty(queue)
    p = [];
    q = [];
    pp = [];
    qq = [];
    Rin = [];
    m = 1;
    for a = 1:size(queue,1)
        pq = queue(1,:);
        queue(1,:) = [];
        p(a) = pq(1);
        q(a) = pq(2);
        if and(p(a) >= 1 & p(a) <= nStepsXP, q(a) >= 1 & q(a) <= nStepsX) && M(p(a),q(a)) == 0
            Rin(:,m) = orbit + [X(p(a),q(a));XP(p(a),q(a));0;0;dpoffset;0]+[1e-9;0;1e-9;0;0;0];
            pp(m) = p(a);
            qq(m) = q(a);
            m=m+1;
        end
    end
    if ~isempty(Rin)
      [~,lost] = ringpass(ring,Rin,nTurns);
    else 
        lost = [];
    end
    for a = 1:numel(lost)
        M(pp(a),qq(a)) = lost(a);
        if M(pp(a),qq(a))
            queue(end+1,1:2) = [pp(a)+1,qq(a)];
            queue(end+1,1:2) = [pp(a)-1,qq(a)];
            queue(end+1,1:2) = [pp(a),qq(a)+1];
            queue(end+1,1:2) = [pp(a),qq(a)-1];
        end
    end
    queue = unique(queue,'rows');
    k=k+1;
end

Xgood = X(~M);
XPgood = XP(~M);
bndry = boundary(Xgood,XPgood);
x_da = Xgood(bndry);
xp_da = XPgood(bndry);


