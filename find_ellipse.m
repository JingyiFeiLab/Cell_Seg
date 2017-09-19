function [ellipticity, tilt, a, c] = find_ellipse(E)

% INPUT : Ellipticity Structure from Fit_ellipse.m
% OUTPUT : Ellipticity (oblate Spheroid)

if isfield(E,'angleToX') == 0
    ellipticity = 0;
    tilt = 0;
    a = 0;
    c = 0;
    return
end
    
a = E.long_axis; % Equatorial Radius
%a = structfun(@double , a , 'uniformoutput', 1);

c = E.short_axis; % Polar Radius
%c = structfun(@double , c , 'uniformoutput', 1);

ellipticity = sqrt((a^2 - c^2)/a^2);
tilt = E.angleToX;

end




















