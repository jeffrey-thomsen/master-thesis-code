function azimuthDeg = interauralToAzimuthMapping(ipdRad, lookuptable)
    
    % by calling the output S and MU, phi is z-scored, thus improving the
    % fitting
    azimuthDeg = polyval(lookuptable.p, ipdRad, lookuptable.S, ...
        lookuptable.MU);
%     % neglect angles > 95°. Warning => maybe systematic underestimation for
%     %  azi ~ 90°
    azimuthDeg(abs(azimuthDeg)>95) = NaN;
end
