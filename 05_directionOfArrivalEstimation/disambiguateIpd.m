function ipdRad = disambiguateIpd(ipdRad, ildDb)
    % set sign of IPD to sign of ILD when abs(ILD) greater than 2.5 dB
    ildUsedLogical = abs(ildDb)>2.5;
    ildUsedSign = sign(ildDb(ildUsedLogical));

    ipdRad(ildUsedLogical) = abs(ipdRad(ildUsedLogical)).*...
        ildUsedSign;
end