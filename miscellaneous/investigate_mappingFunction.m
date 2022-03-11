    for i=[1,5,10,15,20]
        figure(15); plot(azimuthAngleDeg,ipdRadMedianArray{i}); hold on;
        figure; plot(azimuthAngleDeg,ipdRadDisMedianArray{i})
        hold on; plot(azimuthAngleDeg,unwrap(ipdRadMedianArray{i}));
        plot(azimuthAngleDeg,(ipdRadMedianArray{i}));
        plot(azimuthAngleDeg,unwrap(ipdRadMedianArray{i})-2*pi);
    end

    for iBand=1:30
        figure(16); 
        plot(ipdRadMedianArray{iBand},azimuthAngleDeg); 
        hold on; 
        ylim([-100, 100]);
        xlim([-pi, pi]);
        figure(17); 
        plot(linspace(-pi, pi), polyval(lookuptable{iBand}.p,...
            linspace(-pi, pi), lookuptable{iBand}.S, lookuptable{iBand}.MU)); 
        hold on; 
        ylim([-100, 100]);
        xlim([-pi, pi]);
    end

    for iBand=[1,5,10,15,20]
        figure(18); plot(ildMedianArray{iBand},azimuthAngleDeg); hold on; ylim([-100, 100]); xlim([-10, 10]);
        figure(19); plot(linspace(-10,10),polyval(pILD{iBand}, linspace(-pi,pi), SILD{iBand}, ...
            MUILD{iBand})); hold on; ylim([-100, 100]); xlim([-10, 10]);
    end