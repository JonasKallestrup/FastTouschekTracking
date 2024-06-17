function [ap_pos,ap_neg] =  momentumAperture_binary(ring,index,nturns,resolution,dpmax)
xy      = 1e-9;

%%%%%%%%%%% negative aperture
% generate a beam to track
orbit    = findorbit6(ring,index);

disp('Tracking positive momentum aperture')
[ap_pos,~] = runMAcomp(ring,orbit,index,nturns,dpmax,resolution,xy);
disp('Tracking negative momentum aperture')
[ap_neg,~]  = runMAcomp(ring,orbit,index,nturns,-dpmax,-resolution,xy);

end

function [ap,nTurnsTracked] = runMAcomp(ring0,orbit,index,nturns,dpmax,resolution,xy)

values = 0:resolution:dpmax;


nTurnsTracked = 0;
ap = [];
for nn = 1:numel(index)
    low = 1;
    high = numel(values);
    
    ring = atrotatelattice(ring0,index(nn));
    r0 = orbit(:,nn);
    r0([1,3]) = r0([1,3])+xy;
    while floor(high-low) > 1
        k = floor((high-low)/2)+low;
        r1 = r0;
        r1(5) = r1(5) + values(k);
        [~,lost,nTurnsTmp]=ringpass(ring,r1,nturns);

        nTurnsTmp(isnan(nTurnsTmp)) = nturns;
        nTurnsTracked = nTurnsTracked+sum(nTurnsTmp);
        if lost
            high = k;
        else
            low = k;
            ap(nn) = values(k);
        end

    end
end

end


