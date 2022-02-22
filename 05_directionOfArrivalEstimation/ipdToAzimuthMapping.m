function azimuthRad = ipdToAzimuthMapping(ipdRad, lookuptable)
    
    % by calling the output S and MU, phi is z-scored, thus improving the
    % fitting
    azimuthRad = polyval(lookuptable.p, ipdRad, lookuptable.S, ...
        lookuptable.MU);
%     % neglect angles > 95°. Warning => maybe systematic underestimation for
%     %  azi ~ 90°
%     azimuthRad(abs(azimuthRad)>95) = NaN;
end

