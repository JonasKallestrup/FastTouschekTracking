
function [u_out,up_out,theta] = discretisePolygon(u_in,up_in)
u0 = mean(u_in);
up0 = mean(up_in);

theta = linspace(0,2*pi,101);
for a1 = 1:numel(theta)
    u_line = [u0,2*max(abs(u_in))*cos(theta(a1))];
    up_line = [up0,2*max(abs(up_in))*sin(theta(a1))];
    [aa,bb] = intersections(u_line',up_line',u_in,up_in);
    if numel(aa) == 1
        u_out(a1) = aa;
        up_out(a1) = bb;
    else
        D = sqrt(aa.^2+bb.^2);
        u_tmp  = aa(D==min(D));% if there are two crossings, take the one closest to the centre
        up_tmp = bb(D==min(D));
        if any([isempty(u_tmp),isempty(up_tmp)])
            u_out(a1) = nan;
            up_out(a1) = nan;
        else
            u_out(a1) = u_tmp(1);
            up_out(a1) = up_tmp(1);
        end
    end
end
end