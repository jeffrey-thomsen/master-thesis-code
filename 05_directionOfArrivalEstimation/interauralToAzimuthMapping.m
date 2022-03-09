% Map interaural parameters (IPD or ITD) to azimuth angles of incidence in
% degrees using a corresponding mapping function (polynomial)
%
% Input:
% interauralParameter - vector of numeric values (IPD in radians or ITD in
% seconds)
% lookuptable - struct containing the mapping polynomial coefficients p, S
% and MU, generated with the polyfit function in
% interauralToAzimuthLookuptable.m
%
% Output:
% azimuthDeg - vector of azimuth values in degrees
function azimuthDeg = interauralToAzimuthMapping(interauralParameter, lookuptable)
    
    % by calling the output S and MU, phi is z-scored, thus improving the
    % fitting
    azimuthDeg = polyval(lookuptable.p, interauralParameter, lookuptable.S, ...
        lookuptable.MU);

%     % neglect angles > 95°. Warning => maybe systematic underestimation for
%     %  azi ~ 90°
%     azimuthDeg(abs(azimuthDeg)>95) = NaN;
end

