function depth = getUnitDepth(meta,obj)
% given a single sessions obj, return array depth
% that contains for each unit in obj.clu the depth along the shank

% here, we just have the assumed depth given the manipulator software
% coords

% sessions missing ex, I am manually hardcoding experimental data

for iprobe = 1:numel(meta.probe)
    
    prb = meta.probe(iprobe);

    % get obj.clu
    clu = obj.clu{prb};

    % get sites / channel
    if isfield(clu,'site')
        sites = [clu(:).site];
    elseif isfield(clu,'channel')
        sites = [clu(:).channel];
    else
        depth = nan;
        return
    end


    if isstruct(obj.ex)
        % which probe
        probe = obj.ex.probe(prb).type;
        % penetration depth
        penetrationDepth = obj.ex.probe(prb).depth;
    else
        if strcmp(meta.anm,'JEB7') && strcmp(meta.date,'2021-04-29')
            probe = 'H2';
            penetrationDepth = 1200;
        elseif strcmp(meta.anm,'JEB7') && strcmp(meta.date,'2021-04-30')
            probe = 'H2';
            penetrationDepth = 1200;
        elseif strcmp(meta.anm,'EKH3') && strcmp(meta.date,'2021-08-11')
            probe = 'H2';
            penetrationDepth = 1200;
        elseif strcmp(meta.anm,'JGR2') && strcmp(meta.date,'2021-11-16')
            probe = 'H2';
            penetrationDepth = 1200;
        elseif strcmp(meta.anm,'JGR2') && strcmp(meta.date,'2021-11-17')
            probe = 'H2';
            penetrationDepth = 1200;
        elseif strcmp(meta.anm,'JGR3') && strcmp(meta.date,'2021-11-18')
            probe = 'H2';
            penetrationDepth = 1100;
        end
    end



    if strcmpi(probe,'h2')
        siteLoc = [0, 0; 0, 25; 0, 50; 0, 75; 0, 100; 0, 125; 0, 150; 0, 175; 0, 200; 0, 225; 0, 250; 0, 275; 0, 300; 0, 325; 0, 350; 0, 375; 0, 400; 0, 425; 0, 450; 0, 475; 0, 500; 0, 525; 0, 550; 0, 575; 0, 600; 0, 625; 0, 650; 0, 675; 0, 700; 0, 725; 0, 750; 0, 775; 250, 0; 250, 25; 250, 50; 250, 75; 250, 100; 250, 125; 250, 150; 250, 175; 250, 200; 250, 225; 250, 250; 250, 275; 250, 300; 250, 325; 250, 350; 250, 375; 250, 400; 250, 425; 250, 450; 250, 475; 250, 500; 250, 525; 250, 550; 250, 575; 250, 600; 250, 625; 250, 650; 250, 675; 250, 700; 250, 725; 250, 750; 250, 775]; % (formerly mrSiteXY) Site locations (in m) (x values in the first column, y values in the second column)
        siteMap = 1:64;
    elseif strcmpi(probe,'neuropixels')
        siteLoc = [43,20;11,20;59,40;27,40;43,60;11,60;59,80;27,80;43,100;11,100;59,120;27,120;43,140;11,140;59,160;27,160;43,180;11,180;59,200;27,200;43,220;11,220;59,240;27,240;43,260;11,260;59,280;27,280;43,300;11,300;59,320;27,320;43,340;11,340;59,360;27,360;43,380;11,380;59,400;27,400;43,420;11,420;59,440;27,440;43,460;11,460;59,480;27,480;43,500;11,500;59,520;27,520;43,540;11,540;59,560;27,560;43,580;11,580;59,600;27,600;43,620;11,620;59,640;27,640;43,660;11,660;59,680;27,680;43,700;11,700;59,720;27,720;43,740;11,740;59,760;27,760;43,780;11,780;59,800;27,800;43,820;11,820;59,840;27,840;43,860;11,860;59,880;27,880;43,900;11,900;59,920;27,920;43,940;11,940;59,960;27,960;43,980;11,980;59,1000;27,1000;43,1020;11,1020;59,1040;27,1040;43,1060;11,1060;59,1080;27,1080;43,1100;11,1100;59,1120;27,1120;43,1140;11,1140;59,1160;27,1160;43,1180;11,1180;59,1200;27,1200;43,1220;11,1220;59,1240;27,1240;43,1260;11,1260;59,1280;27,1280;43,1300;11,1300;59,1320;27,1320;43,1340;11,1340;59,1360;27,1360;43,1380;11,1380;59,1400;27,1400;43,1420;11,1420;59,1440;27,1440;43,1460;11,1460;59,1480;27,1480;43,1500;11,1500;59,1520;27,1520;43,1540;11,1540;59,1560;27,1560;43,1580;11,1580;59,1600;27,1600;43,1620;11,1620;59,1640;27,1640;43,1660;11,1660;59,1680;27,1680;43,1700;11,1700;59,1720;27,1720;43,1740;11,1740;59,1760;27,1760;43,1780;11,1780;59,1800;27,1800;43,1820;11,1820;59,1840;27,1840;43,1860;11,1860;59,1880;27,1880;43,1900;11,1900;59,1920;27,1920;43,1940;11,1940;59,1960;27,1960;43,1980;11,1980;59,2000;27,2000;43,2020;11,2020;59,2040;27,2040;43,2060;11,2060;59,2080;27,2080;43,2100;11,2100;59,2120;27,2120;43,2140;11,2140;59,2160;27,2160;43,2180;11,2180;59,2200;27,2200;43,2220;11,2220;59,2240;27,2240;43,2260;11,2260;59,2280;27,2280;43,2300;11,2300;59,2320;27,2320;43,2340;11,2340;59,2360;27,2360;43,2380;11,2380;59,2400;27,2400;43,2420;11,2420;59,2440;27,2440;43,2460;11,2460;59,2480;27,2480;43,2500;11,2500;59,2520;27,2520;43,2540;11,2540;59,2560;27,2560;43,2580;11,2580;59,2600;27,2600;43,2620;11,2620;59,2640;27,2640;43,2660;11,2660;59,2680;27,2680;43,2700;11,2700;59,2720;27,2720;43,2740;11,2740;59,2760;27,2760;43,2780;11,2780;59,2800;27,2800;43,2820;11,2820;59,2840;27,2840;43,2860;11,2860;59,2880;27,2880;43,2900;11,2900;59,2920;27,2920;43,2940;11,2940;59,2960;27,2960;43,2980;11,2980;59,3000;27,3000;43,3020;11,3020;59,3040;27,3040;43,3060;11,3060;59,3080;27,3080;43,3100;11,3100;59,3120;27,3120;43,3140;11,3140;59,3160;27,3160;43,3180;11,3180;59,3200;27,3200;43,3220;11,3220;59,3240;27,3240;43,3260;11,3260;59,3280;27,3280;43,3300;11,3300;59,3320;27,3320;43,3340;11,3340;59,3360;27,3360;43,3380;11,3380;59,3400;27,3400;43,3420;11,3420;59,3440;27,3440;43,3460;11,3460;59,3480;27,3480;43,3500;11,3500;59,3520;27,3520;43,3540;11,3540;59,3560;27,3560;43,3580;11,3580;59,3600;27,3600;43,3620;11,3620;59,3640;27,3640;43,3660;11,3660;59,3680;27,3680;43,3700;11,3700;59,3720;27,3720;43,3740;11,3740;59,3760;27,3760;43,3780;11,3780;59,3800;27,3800;43,3820;11,3820;59,3840;27,3840]; % (formerly mrSiteXY) Site locations (in m) (x values in the first column, y values in the second column)
        siteMap = 1:384;
        sites(sites==385) = 384;
    end

    yLocOnShank = siteLoc(sites,2);
    depth{iprobe} = penetrationDepth - yLocOnShank;

    % figure; plot(locOnShank); hold on; plot(locOnShank,'.','MarkerSize',15)
    % figure; plot(depth); hold on; plot(depth,'.','MarkerSize',15)


end

if numel(meta.probe)==1
    depth = depth{1};
end

end % getDepth()







