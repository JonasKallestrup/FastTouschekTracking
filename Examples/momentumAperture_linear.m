function [ap_pos,ap_neg] = momentumAperture_linear(ring,index,nturns,dpstep,dpmax)
% Calculate the momentum aperture for specified indices
xy      = 1e-9;

%%%%%%%%%%% negative aperture
% generate a beam to track
orbit    = findorbit6(ring,index);

disp('Tracking positive momentum aperture')
[ap_pos,dp_tested_pos] = runMAcomp(ring,orbit,index,nturns,dpstep,dpmax,xy);
disp('Tracking negative momentum aperture')
[ap_neg,dp_tested_neg] = runMAcomp(ring,orbit,index,nturns,-dpstep,-dpmax,xy);

end

function [ap,dp_tested] = runMAcomp(ring0,orbit,index,nturns,dpstep,dpmax,xy)

ap = zeros(size(index));
for nn = 1:numel(index)
    dp_tested{nn} = [];
    dp_test = dpstep:dpstep:dpmax;
    ring = atrotatelattice(ring0,index(nn));
    
    r0= zeros(6,numel(dp_test));
    r0([1,3],:) = xy;
    r0(5,:) = dp_test;
    r0 = r0+orbit(:,nn);

    
    [~,lost]=ringpass(ring,r0,nturns);
    dp_tested{nn} = [dp_tested{nn},dp_test];

    idx_surv =  find(lost,1,'first')-1;
    if  or(isempty(idx_surv),idx_surv ==0)
        ap(nn) = 0;
    else
        ap(nn) = dp_test(idx_surv);
    end
end
end


