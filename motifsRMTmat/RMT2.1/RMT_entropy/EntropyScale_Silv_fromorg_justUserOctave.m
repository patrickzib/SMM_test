function [ SS,depd, IDMall,TIMESCALE] = EntropyScale_Silv_fromorg_justUserOctave(I,LocM,Ot,Od,St,Sd,ominT,ominD,sminT,sminD,smaxD,smaxT,sigmaT0,sigmaD0,DepThreshold,dsigmaT0,dsigmaD0)%
%I : timeseries column are variate rows are timedata,
%LocM: location matrix
% sigmaNT :  Nominal smoothing of the input timeseries across time
% sigmaND :  Nominal smoothing of the input timeseries Dependency
% Ot      :  Numeber of desired octave Time
% Od      :  Numeber of desired octave Dependency
% St      :  Number of maximum scale over Time
% Sd      :  Number of maximum scale over Dependency
% ominT   :  usually setted to 0 is the minimum octave of time
% ominD   :  usually fixed to 0 is the minimum octave for the dependency
% sminT   :  minimum  scale over time (it require to compute an offset is <0)
% sminD   :  minimum  scale over dependency (it require to compute an offset is <0)
% smaxD   :  minimum  scale over time usually >=3
% smaxT   :  minimum  scale over time usually >=3
% sigmaT0 : Smoothing of the level 0 of octave 0 of the scale space. (Note that Lowe's 1.6 value refers to the level -1 of octave 0.)
% sigmaD0 : Smoothing of the level 0 of the dependency
% dsigmaT0: step between the scale time if ==1 then the scale
% dsigmaD0: step between the scale dependency
% defineStepFactor: it is a flag to define a diferent step factor between different scales.

% Scale multiplicative step
TIMESCALE= zeros(1,4);
timestart1=0;
tic
ktime = 2^(1/St) ;
kdepd = 2^(1/Sd) ;

if sigmaT0 <0.5
    sigmaT0 = 1.6 * ktime ;
end

if sigmaD0 <=0
    sigmaD0 = 0.6 ;
end

if dsigmaT0 < 0
    dsigmaT0 = sigmaT0 * sqrt(1 - 1/ktime^2) ; % Scale step factor Time between each scale
end
if dsigmaD0 < 0
    dsigmaD0 = sigmaD0 * sqrt(1 - 1/kdepd^2) ; % Scale step factor Dependency between each scale
    %   dsigmaD0 = sigmaD0;
end

sigmaND =0.5;
sigmaNT =0.5;

% Scale space construction
% Save parameters
SS.Od          = Od;
SS.Ot          = Ot;
SS.St         = St;
SS.Sd         = Sd;

SS.sigmat     = sigmaT0;
SS.sigmad     = sigmaD0;
SS.odmin       = 0;
SS.otmin       = 0;
SS.sminT       = sminT ;
SS.smaxT       = smaxT ;
SS.sminD       = sminD ;
SS.smaxD       = smaxD ;
%% Starting from negative octave
% if otmin < 0
% 	for o=1:-omin
% 		I = doubleSize(I) ;
% 	end
% elseif otmin > 0
% 	for o=1:omin
% 		I = halveSize(I) ;
% 	end
% end

% starting octave
otcur = 1;
odcur = 1;

% Data size
[M, N] = size(I);

% Index offset
soT = -sminT+1 ;
soD = -sminD+1 ;

% For dependency matrix
IDM = (1:N)';
IDMall{odcur} = IDM;
DistM = computeDist(LocM, IDM, odcur);         %DistM: distance matrixLocM;%

%% Silv Normalization
maxim = max(DistM(:));
minim = min(DistM(:));
DistM= (DistM - minim) /abs(maxim-minim);
%silv start to modify
isover=false;
otact=1;
odact=1;
tic;

H = dependencyALL(DistM , DepThreshold, IDM, odcur);
depd{odcur} = H;
SS.ds{otcur, odcur} = [1, 1];

SS.octave{otcur, odcur} = zeros(M, N,smaxD,smaxT) ;
SS.smoothmatrix{otcur, odcur} = zeros(N, N, smaxD);
% From first octave
SDepdsigmafor_OT1_OD1 =sigmaD0; %sqrt((sigmaD0*kdepd^sminD)^2  - (sigmaND/2^ominD)^2);%
Smatrix = ComputeDependencyScale(depd{odcur}, SDepdsigmafor_OT1_OD1); %% can be chaged to
%supposing that the  data are presmoothed with a sigmaNT over time we can
%do a smoothing as:
STimegsigmafor_OT1_OD1 = sqrt((sigmaT0*ktime^sminT)^2  - (sigmaNT/2^ominT)^2);
[SS.octave{otcur,odcur}(:,:,1,1),SS.smoothmatrix{otcur,odcur}(:,:,1,1)] = smooth(I, Smatrix, STimegsigmafor_OT1_OD1);
EntropyInputData=I;
InputEntropyQuaantized= squeeze(SS.octave{otcur,odcur}(:,:,1,1));
SS.Entropyoctave{otcur,odcur}(:,:,1,1) = computeEntropyScale_1(InputEntropyQuaantized,STimegsigmafor_OT1_OD1,Smatrix,DepThreshold);

if (otact-1)~=0
    SS.ds{otact, odact} = [SS.ds{otact-1, odact}(1)+1, SS.ds{otact-1, odact}(2)];
end
LocMTemp = LocM;

doubleSigmaTime = dsigmaT0*St;%    7.800628871;%9.620988313;%sigmaT0*2;
doubleSigmaDepd = dsigmaD0*Sd;%1.378969393;%1.70076652;%sigmaD0*2;
COMEFROM=1;% 1 BOTH; 2=DEPD; 3=TIME;

while(~isover)
    if (otact-1)~=0
        SS.ds{otact, odact} = [SS.ds{otact-1, odact}(1)+1, SS.ds{otact-1, odact}(2)];
    end
    LocMTemp = LocM;
    
    if(otact<Ot & odact<Od)
        % smooth of 2sigmadepd and 2sigmatime over both
        
        Smatrix =ComputeDependencyScale(depd{odcur}, doubleSigmaDepd);
        SS = Smooth_Asyn_justDiagonal(SS, otact, odact, ktime, kdepd, sminT, sminD,smaxT, smaxD ,dsigmaT0,      dsigmaD0,        depd, Smatrix,soT,soD);
        %         [SS.octave{otact,odact}(:,:,Sd,St),SS.smoothmatrix{otcur,odcur}(:,:,Sd)] = smooth(squeeze(SS.octave{otact,odact}(:,:,1,1)), Smatrix, doubleSigmaTime);
        % increase the octave time and depd
        
        otact=otact+1;
        odact=odact+1;
        % reduce  bothsize
        sbest_time = min(sminT + St, smaxT) ;
        sbest_Dep = min(sminD + Sd, smaxD) ;
        TMP=halveSizeTime(squeeze(SS.octave{otact-1,odact-1}(:,:,sbest_Dep+soD,sbest_time+soT)));
        target_sigmaT = sigmaT0 * ktime^sminT ;
        prev_sigmaT = sigmaT0 * ktime^(sbest_time - St) ;
        target_sigmaD = sigmaD0 * kdepd^sminD ;
        prev_sigmaD = sigmaD0 * kdepd^(sbest_Dep - Sd) ;
        if(target_sigmaD > prev_sigmaD)
            Smatrix=ComputeDependencyScale(depd{odcur},sqrt(target_sigmaD^2 - prev_sigmaD^2));
            TMP = smoothJustDependencySilv(TMP, sqrt(target_sigmaT^2 - prev_sigmaT^2),Smatrix ) ;
        end
        if(target_sigmaT > prev_sigmaT)
            TMP = smoothJustTimeSilv(TMP, sqrt(target_sigmaT^2 - prev_sigmaT^2),Smatrix ) ;
        end
        [SS.octave{otact,odact}(:,:,1,1), LocMTemp, IDM] = HalfDependencyMote(LocMTemp, squeeze(TMP));
        DistM = computeDist(LocMTemp, IDM, odcur);
        %% Silv Normalization
        maxim = max(DistM(:));
        minim = min(DistM(:));
        DistM= (DistM - minim) /abs(maxim-minim);
        %% Silv Normalization
        DepThreshold = DepThreshold*1.1;
        H = dependencyALL(DistM , DepThreshold, IDM, odcur);
        
        %         TMP=halveSizeTime(squeeze(SS.octave{otact-1,otact-1}(:,:,Sd,St)));
        %         [SS.octave{otact,odact}(:,:,1,1), LocMTemp, IDM] = HalfDependencyMote(LocMTemp, squeeze(TMP));
        %         DistM = computeDist(LocMTemp, IDM, odcur);
        %         %% Silv Normalization
        %         maxim = max(DistM(:));
        %         minim = min(DistM(:));
        %         DistM= (DistM - minim) /abs(maxim-minim);
        %         %% Silv Normalization
        %         DepThreshold = DepThreshold*1.1;
        %         H = dependencyALL(DistM , DepThreshold, IDM, odcur);
        if(otact==Ot & odact==Od)
            isover=true;
        end
        COMEFROM=1;
    elseif(otact==Ot & odact<Od)
        %         smooth of 2 sigmadepd over depd
        %         Smatrix=ComputeDependencyScale(depd{odcur},doubleSigmaDepd);
        %        [SS.octave{otact,odact}(:,:,Sd,1)] = smoothJustDependencySilv(squeeze(SS.octave{otact,odact}(:,:,1,1)), 0,Smatrix ) ;
        SS = Smooth_Asyn_justDepd(SS, otact, odact, ktime, kdepd, sminT, sminD,smaxT, smaxD ,dsigmaT0,      dsigmaD0,        depd, Smatrix,soT,soD);
        % increase the octave depd
        odact=odact+1;
        % reduce  over dependency size
        sbest_Dep = min(sminD + Sd, smaxD) ;
        %half size dependency
        [SS.octave{otact,odact}(:,:,1,1), LocMTemp, IDM] = HalfDependencyMote(LocMTemp, squeeze(SS.octave{otact,odact-1}(:,:,sbest_Dep+soD,1)));
        %distance matrix
        DistM = computeDist(LocMTemp, IDM, odcur);
        %% Silv Normalization of distance matrix
        maxim = max(DistM(:));
        minim = min(DistM(:));
        DistM= (DistM - minim) /abs(maxim-minim);
        %% Silv Normalization
        DepThreshold = DepThreshold*1.1;
        H = dependencyALL(DistM , DepThreshold, IDM, odcur);
        
        %        [SS.octave{otact,odact}(:,:,1,1), LocMTemp, IDM] = HalfDependencyMote(LocMTemp, squeeze(SS.octave{otact,odact}(:,:,Sd,1)));
        %        DistM = computeDist(LocMTemp, IDM, odcur);
        %        %% Silv Normalization
        %        maxim = max(DistM(:));
        %        minim = min(DistM(:));
        %        DistM= (DistM - minim) /abs(maxim-minim);
        %        %% Silv Normalization
        %        DepThreshold = DepThreshold*1.1;
        %        H = dependencyALL(DistM , DepThreshold, IDM, odcur);
        if(otact==Ot & odact==Od)
            isover=true;
        end
        COMEFROM=2;
    elseif(otact<Ot & odact==Od)
        %         smooth of 2 sigmatime over time
        SS = Smooth_Asyn_justTime(SS, otact, odact, ktime, kdepd, sminT, sminD,smaxT, smaxD ,dsigmaT0,      dsigmaD0,        depd, Smatrix,soT,soD);
        % [SS.octave{otact,odact}(:,:,1,St)] = smoothJustTimeSilv(TMP, doubleSigmaTime,Smatrix ) ;
        % increase the octave time
        otact=otact+1;
        % reduce  over time size
        sbest_time = min(sminT + St, smaxT) ;
        %half size of time
        TMP= halveSizeTime(squeeze(SS.octave{otact-1,odact}(:,:,1,sbest_time+soT)));
        target_sigmaT = sigmaT0 * ktime^sminT ;
        prev_sigmaT = sigmaT0 * ktime^(sbest_time - St) ;
        if(target_sigmaT > prev_sigmaT)
            TMP = smoothJustTimeSilv(TMP, sqrt(target_sigmaT^2 - prev_sigmaT^2),Smatrix ) ;
        end
        SS.octave{otact,odact}(:,:,1,1)=TMP;
        %     SS.octave{otact,odact}(:,:,1,1)= halveSizeTime(squeeze(SS.octave{otact,odact}(:,:,1,St)));
        if(otact==Ot & odact==Od)
            isover=true;
        end
        %     COMEFROM=3;
        
    end
    
    if(odact~=1 && size(depd,2)<odact)
        IDMall{odact} = IDM;
        depd{odact} = H;
    end
    SS.smoothmatrix{otact, odact} = zeros(size(SS.octave{otact, odact},2),size(SS.octave{otact, odact},2),size(SS.octave{otact, odact},3));
    if odact ~= 1
        SS.ds{otact, odact} = [otact, odact];%[SS.ds{otact, odact-1}(1), SS.ds{otact, odact-1}(2)+1];
    end
end
% SS.smoothmatrix{otact, odact} = zeros(size(SS.octave{otact, odact},2),size(SS.octave{otact, odact},2),size(SS.octave{otact, odact},3));
SS = Smooth_Asyn_Entropy(SS, otact, odact, ktime, kdepd, sminT, sminD,smaxT, smaxD ,dsigmaT0, dsigmaD0, depd, Smatrix,soT,soD,DepThreshold,EntropyInputData);

TIMESCALE = toc;

function [SS] = Smooth_Asyn_Entropy(SS, CurrentTimeOct, CurrentDepdOct, ktime, kdepd,stmin,sdmin, stmax, sdmax,sigmaTscaleStep,sigmaDscaleStep,depd, Smatrix,soT,soD,DepThreshold,EntropyInputData)
for sd=sdmin:sdmax
         dsigmaD = kdepd^(sd) * sigmaDscaleStep ;%kdepd^(sd+1) * sigmaDscaleStep ;
         Smatrix = ComputeDependencyScale(depd{CurrentDepdOct}, dsigmaD);
    %     SS.smoothmatrix{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD) = Smatrix;
    if sd==(sdmin)
        for st=stmin+1:stmax
            dsigmaT =  ktime^(st) * sigmaTscaleStep ;%ktime^(st+1) * sigmaTscaleStep ;
            if st== (stmin)
                % this scale is already computed
            else
                [SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD,st+soT)] = smoothJustTimeSilv(squeeze(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD, st+soT-1)),dsigmaT,Smatrix);
                %smoothTime(squeeze(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD, st+soT-1)),dsigmaT,Smatrix);
                InputEntropyQuaantized= SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD,st+soT);%globalQuantization(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD,st+soT));%singlevariateQuantization(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD,st+soT));
                [SS.Entropyoctave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD,st+soT)] = computeEntropyScale_1(InputEntropyQuaantized,dsigmaT,Smatrix,DepThreshold);
                %                ...
                %                                                                 computeEntropyScale_1(squeeze(round(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD,st+soT))),dsigmaT,Smatrix,DepThreshold);%computeEntropyScale_1(EntropyInputData,dsigmat,Smatrix,threshold);%
            end
        end
    else
        dsigmaD = kdepd^(sd) * sigmaDscaleStep ;%kdepd^(sd+1) * sigmaDscaleStep ;
        Smatrix = ComputeDependencyScale(depd{CurrentDepdOct}, dsigmaD);
        SS.smoothmatrix{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD) = Smatrix;
        for st=stmin:stmax
            dsigmaT = ktime^(st) * sigmaTscaleStep ;
            if st== (stmin)
                % compute for the first eelment of the border the smoothing
                % just using the  dependency
                SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD, st+soT) = smoothJustDependencySilv(squeeze(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD-1, st+soT)),dsigmaT,Smatrix);
                %smoothTime(squeeze(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD-1, st+soT)),0.5,Smatrix);
                InputEntropyQuaantized= SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD,st+soT);%globalQuantization(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD,st+soT));%singlevariateQuantization(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD,st+soT));
                [SS.Entropyoctave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD, st+soT)] = computeEntropyScale_1(InputEntropyQuaantized,dsigmaT,Smatrix,DepThreshold);
                %                ...
                %                                                                 computeEntropyScale_1(squeeze(round(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD, st+soT))),dsigmaT,Smatrix,DepThreshold);%computeEntropyScale_1(EntropyInputData,dsigmat,Smatrix,threshold);%
                
            else
                [SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD,st+soT)] = smoothBothSilv(squeeze(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD-1, st+soT-1)),dsigmaT,Smatrix);
                %smoothTime(squeeze(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD-1, st+soT-1)),dsigmaT,Smatrix);
                InputEntropyQuaantized= globalQuantization(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD,st+soT));%singlevariateQuantization(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD,st+soT));
                [SS.Entropyoctave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD, st+soT)] = computeEntropyScale_1(InputEntropyQuaantized,dsigmaT,Smatrix,DepThreshold);
                %                ...
                %                                                                 computeEntropyScale_1(squeeze(round(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD, st+soT))),dsigmaT,Smatrix,DepThreshold);%computeEntropyScale_1(EntropyInputData,dsigmat,Smatrix,threshold);%
            end
        end
    end
end

function [SS] = Smooth_Asyn_justDepd(SS, CurrentTimeOct, CurrentDepdOct, ktime, kdepd,stmin,sdmin, stmax, sdmax,sigmaTscaleStep,sigmaDscaleStep,depd, Smatrix,soT,soD)
for sd=sdmin:sdmax
    st=0;
    if st== (stmin) & sd==(sdmin)
    else
        dsigmaT =  ktime^(st) * sigmaTscaleStep ;
        dsigmaD = kdepd^(sd) * sigmaDscaleStep ;%kdepd^(sd+1) * sigmaDscaleStep ;
        Smatrix = ComputeDependencyScale(depd{CurrentDepdOct}, dsigmaD);
        SS.smoothmatrix{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD) = Smatrix;
        SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD, st+soT) = smoothJustDependencySilv(squeeze(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD-1, st+soT)),dsigmaT,Smatrix);
    end
end

function [SS] = Smooth_Asyn_justTime(SS, CurrentTimeOct, CurrentDepdOct, ktime, kdepd,stmin,sdmin, stmax, sdmax,sigmaTscaleStep,sigmaDscaleStep,depd, Smatrix,soT,soD)
for st=stmin:stmax
    sd=0;
    if st== (stmin) & sd==(sdmin)
    else
        dsigmaT =  ktime^(st) * sigmaTscaleStep ;
        dsigmaD = kdepd^(sd) * sigmaDscaleStep ;%kdepd^(sd+1) * sigmaDscaleStep ;
        Smatrix = ComputeDependencyScale(depd{CurrentDepdOct}, dsigmaD);
        SS.smoothmatrix{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD) = Smatrix;
        SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD, st+soT) = smoothJustTimeSilv(squeeze(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD-1, st+soT)),dsigmaT,Smatrix);
    end
end


function [SS] = Smooth_Asyn_justDiagonal(SS, CurrentTimeOct, CurrentDepdOct, ktime, kdepd,stmin,sdmin, stmax, sdmax,sigmaTscaleStep,sigmaDscaleStep,depd, Smatrix,soT,soD)
for sd=sdmin:sdmax
    st=sd;
    if st== (stmin) & sd==(sdmin)
    else
        dsigmaT =  ktime^(st) * sigmaTscaleStep ;
        dsigmaD = kdepd^(sd) * sigmaDscaleStep ;%kdepd^(sd+1) * sigmaDscaleStep ;
        Smatrix = ComputeDependencyScale(depd{CurrentDepdOct}, dsigmaD);
        SS.smoothmatrix{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD) = Smatrix;
        [SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD,st+soT)] = smoothBothSilv(squeeze(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD-1, st+soT-1)),dsigmaT,Smatrix);
    end
end

function [SS] = Smooth_Asyn(SS, CurrentTimeOct, CurrentDepdOct, ktime, kdepd,stmin,sdmin, stmax, sdmax,sigmaTscaleStep,sigmaDscaleStep,depd, Smatrix,soT,soD)
for sd=sdmin:sdmax
    if sd==(sdmin)
        for st=stmin+1:stmax
            dsigmaT =  ktime^(st) * sigmaTscaleStep ;%ktime^(st+1) * sigmaTscaleStep ;
            if st== (stmin)
                % this scale is already computed
            else
                [SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD,st+soT)] = smoothJustTimeSilv(squeeze(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD, st+soT-1)),dsigmaT,Smatrix);
                %smoothTime(squeeze(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD, st+soT-1)),dsigmaT,Smatrix);
            end
        end
    else
        dsigmaD = kdepd^(sd) * sigmaDscaleStep ;%kdepd^(sd+1) * sigmaDscaleStep ;
        Smatrix = ComputeDependencyScale(depd{CurrentDepdOct}, dsigmaD);
        SS.smoothmatrix{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD) = Smatrix;
        for st=stmin:stmax
            dsigmaT = ktime^(st) * sigmaTscaleStep ;
            if st== (stmin)
                % compute for the first eelment of the border the smoothing
                % just using the  dependency
                SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD, st+soT) = smoothJustDependencySilv(squeeze(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD-1, st+soT)),dsigmaT,Smatrix);
                %smoothTime(squeeze(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD-1, st+soT)),0.5,Smatrix);
            else
                [SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD,st+soT)] = smoothBothSilv(squeeze(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD-1, st+soT-1)),dsigmaT,Smatrix);
                %smoothTime(squeeze(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD-1, st+soT-1)),dsigmaT,Smatrix);
            end
        end
    end
end
%% This code create with sdmax=6 just 8 scale not 8  because the first scale is dropped.
% for sd=sdmin+1:sdmax
%     dsigmaD = kdepd^(sd+1) * sigmaDscaleStep ;
%     Smatrix = ComputeDependencyScale(depd{CurrentDepdOct}, dsigmaD);
%     SS.smoothmatrix{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD) = Smatrix;
%     if sd==(sdmin+1)
%         for st=stmin+1:stmax
%             dsigmaT = ktime^(st+1) * sigmaTscaleStep ;
%             if st== (stmin+1)
%                % this scale is already computed
%             else
%                [SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD-1,st+soT-1)] = smoothTime(squeeze(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD-1, st)),dsigmaT,Smatrix);
%             end
%         end
%     else
%         for st=stmin+1:stmax
%             dsigmaT = ktime^(st+1) * sigmaTscaleStep ;
%             if st== (stmin+1)
%                % compute for the first eelment of the border the smoothing
%                % just using the  dependency
%                SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD-1, st+soT-1) = smoothTime(squeeze(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd, st+soT-1)),0.5,Smatrix);
%             else
%                [SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd+soD-1,st+soT-1)] = smoothTime(squeeze(SS.octave{CurrentTimeOct,CurrentDepdOct}(:,:,sd, st)),dsigmaT,Smatrix);

%             end
%         end
%     end
% end
%% Sicong Function
function H = dependencyALL(Mat, t, IDM, odcur)
% function of generating dependency matrix.
%  -M is the distance matrix of 53 sensors
%  -t is the threshold for the selected scale
%  -H will be a 0-1 sparse matrix containning neighborhood information
Numberof = size(Mat,1);
H = zeros([Numberof Numberof]);
for i = 1:Numberof
    for j = 1:Numberof
        if(Mat(i,j)<=t)
            H(i,j) = 1;
        end
    end
end
H = H - eye([Numberof Numberof]);

function DistM = computeDist(LocM, IDM, odcur)
Num = size(LocM,1);
DistM = zeros(Num, Num);
for i=1:Num
    for j=1:Num
        DistM(i,j) = norm(LocM(i,:)-LocM(j,:),2);
    end
end

function J = doubleSizeTime(I)
[M,N]=size(I) ;
J = zeros(2*M,N) ;
J(1:2:end,:) = I ;
J(2:2:end-1,:) = ...
    0.25*I(1:end-1,:) + ...
    0.25*I(2:end,:) + ...
    0.25*I(1:end-1,:) + ...
    0.25*I(2:end,:) ;
J(2:2:end-1,:) = ...
    0.5*I(1:end-1,:) + ...
    0.5*I(2:end,:) ;
J(1:2:end,:) = ...
    0.5*I(:,1:end-1) + ...
    0.5*I(:,2:end) ;
