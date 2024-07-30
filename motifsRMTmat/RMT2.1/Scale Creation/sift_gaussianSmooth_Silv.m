function [tempFrames,descriptor_silv,gss,dogss,depd,idm, TIMESCALE, TIMEKEYPOINTS, TIMEDESCRIPTORS]=sift_gaussianSmooth_Silv(I, LocM, Ot, Od, St, Sd, sigmaTime ,sigmaDepd, NBP, gthresh, r,sBoundary, eBoundary)
    %, featureDepdScale, featureTimeScale]
    % [M,N,C] = size(I) ;
    % O = floor(log2(min(M,N)))-2 ; % up to 8x8 images
    % time  = zeros(1, Ot*Od);
    time = zeros(Ot, Od);
    timeDescr = zeros(Ot, Od);
    timee = zeros(1,2);

    featureTimeScale = [];
    featureDepdScale = [];
    % thresh = 0.04 / St / 2 ;
    thresh = 0.04 / 3 / 2 ; %value picked from the vidaldi code
    NBO    = 8;
    NBP_Time = 4;
    NBP_Depd = 4;
    magnif = 3.0;
    tempFrames = [];
    descriptors = [] ;
    descriptor_silv=[];
    ktime = 2^(1/(St-3));
    kdepd = 2^(1/(Sd-3));

    % Try this function
    stmin=0;%-1;
    sdmin=0;%-1;
    otmin=0;
    odmin=0;

    [gss, depd, idm,TIMESCALE] = gaussianss_Silv_fromorg_justUserOctave(I, LocM, Ot,Od,St,Sd,otmin,odmin,stmin,sdmin,St+1,Sd+1, sigmaTime, sigmaDepd,gthresh,-1,-1);

    TIMEKEYPOINTS=zeros(1,4);
    dogss = diffss_asynchronous_justUserOctave(gss,Od,Ot);

    otime = Ot;
    odepd = Od;

    tic;
    % Local maxima of the DOG octave
    scaleDiff = St - 1;

    %% used
    forwardIdx = siftlocalmax_directed_100(dogss.octave{otime, odepd}{3},dogss.octave{otime, odepd}{2},dogss.octave{otime, odepd}{1}, 0.8*thresh, NormalizeF(depd{odepd}), NormalizeB(depd{odepd}'), scaleDiff);

    [i,j, s1] = ind2sub( size( dogss.octave{otime, odepd}{3}), forwardIdx );
    y=i-1;
    x=j-1;
    s1=s1-1+gss.sminT;
    % s1=s1-1-1 ; SICONG
    s2 = s1;
    forwardIdx = [x(:)';y(:)';s1(:)'];

    dogss.octave{otime, odepd}{1} = -dogss.octave{otime, odepd}{1};
    dogss.octave{otime, odepd}{2} = -dogss.octave{otime, odepd}{2};
    dogss.octave{otime, odepd}{3} = -dogss.octave{otime, odepd}{3};

    %% used
    backwardIdx = siftlocalmax_directed_100(dogss.octave{otime, odepd}{3},dogss.octave{otime, odepd}{2},dogss.octave{otime, odepd}{1}, 0.8*thresh, NormalizeF(depd{odepd}), NormalizeB(depd{odepd}'), scaleDiff);

    [i,j, s1] = ind2sub( size( gss.octave{otime, odepd} ), backwardIdx ) ;
    y=i-1;
    x=j-1;
    s1=s1-1+gss.sminT;
    s2 = s1;
    backwardIdx = [x(:)';y(:)';s1(:)'] ;

    oframes = [forwardIdx, backwardIdx];
    [C, ia, ic] = unique(oframes', 'rows');
    oframes = C';
    for hsize = 1:size(dogss.octave{otime, odepd}{3},3)% iterate on the scale in the specific octave
        HY(:,:,hsize) = (NormalizeByRow(depd{odepd})*(-dogss.octave{otime, odepd}{3}(:,:,hsize)'))';
        HY2(:,:,hsize) = (NormalizeByRow(depd{odepd}')*(-dogss.octave{otime, odepd}{3}(:,:,hsize)'))';
    end
    rad = 0;
    sel= ... % select feature that are centered in the timeseries
        oframes(2,:)-rad >= 1  & ...
        oframes(2,:)+rad <= size(gss.octave{otime, odepd},1)      ;
    oframes=oframes(:,sel) ;

    oframes = siftrefinemx_directed(oframes, -dogss.octave{otime, odepd}{3},HY,HY2,gss.sminT,thresh,r,0) ;

    TIMEKEYPOINTS(otime + odepd)= toc;
    clear HY HY2
    TIMEDESCRIPTORS=zeros(1,otime + odepd);
    tic;
    if size(oframes,2) >0
        pricurRatio = oframes(4,:);%zeros(1,size(oframes,2));%
    else
        pricurRatio = zeros(1,0);
    end

    if(size(oframes, 2) ~=0)
        oframes(4,:) = odepd; % Save the Octave Dependency
        oframes(5,:) = otime; % Save the Octave Time
        % Store frames
        x = oframes(1, :); % the original code report to the variate of the  specific scale 2^(o-1+gss.omin) * oframes(1,:) ;
        y  = 2^(gss.ds{otime, odepd}(1)+gss.otmin-1) * oframes(2,:) ; % otmin starts from 0 for timeseries
        %report tievalue to the original scale
        
        tempDepd = oframes(1,:) ;
        tempTime = oframes(2,:) ;
        dependencyScale = oframes(3,:);
        %dependencyScale = oframes(3,:)+1; %SICONG
        % timeScale = oframes(4,:)+1;
        timeScale = oframes(3,:);% SICONG
        
        sigmad =  2^(odepd-1+gss.odmin) * gss.sigmad * 2.^(oframes(3,:)/gss.Sd);
        sigmat =  2^(otime-1+gss.odmin) * gss.sigmat * 2.^(oframes(3,:)/gss.Sd);
        
        % append difference-of-Gausssian values to output
        TimescaleSicongNormalized = timeScale+1;
        [timeDoGs, depdDoGs, bothDoGs] = appendDogs(dogss.octave{otime, odepd}, tempDepd, tempTime, dependencyScale, TimescaleSicongNormalized);

        %% unique goes here
        tempFrames = [tempFrames, [x(:)'+ones(1,size(x,1)) ; y(:)' ; sigmad(:)' ;sigmat(:)' ; oframes(4,:); oframes(5,:); oframes(3,:);pricurRatio; timeDoGs(:)'; depdDoGs(:)'; bothDoGs(:)']];%[x(:)'

    end

    dogss.octave{otime, odepd}{1} = -dogss.octave{otime, odepd}{1};
    dogss.octave{otime, odepd}{2} = -dogss.octave{otime, odepd}{2};
    dogss.octave{otime, odepd}{3} = -dogss.octave{otime, odepd}{3};
    
    % 1 means directed graph
    [fgss_silv,Pseudo_centerVaraite]= computeFeatureMatrix_Silv(gss.octave{otime, odepd},gss.sminT,gss.sminD,gss.sigmad,gss.St,gss.Sd,NormalizeByRow(depd{odepd}),NormalizeByRow(depd{odepd}'),1);

    % Descriptors
    if(size(oframes, 2) > 0)
        for f=1:size(oframes,2)
            oframeSilv=oframes(:,f);
            centerV= Pseudo_centerVaraite(1,oframeSilv(3,1)-gss.sminD+1);
            oframeSilv(1,1)=centerV;
            oframeSilv(4,1)=oframeSilv(3,1);
            oframeSilv(5,1)=0;
            sh_silv=siftdescriptor_Silv(...
                fgss_silv{oframeSilv(3,1)-gss.sminD+1,oframeSilv(4,1)-gss.sminT+1,oframes(1,f)+1},...
                oframeSilv(:,1),...%oframes(:,f), ...
                gss.sigmat, ...
                gss.sigmad	,...
                gss.St, ...
                gss.Sd, ...
                gss.sminT	, ...
                gss.sminD, ...
                magnif, ...
                NBP_Time, ...
                NBP_Depd, ...
                NBO) ;
            descriptor_silv=[descriptor_silv,sh_silv];                
        end          
    end
    TIMEDESCRIPTORS(otime + odepd)= toc;
    clear fOframes bOframes fgss    
end
