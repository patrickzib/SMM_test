function dss = diffss_asynchronous_Xiaolan(ss)
dss.otmin = ss.otmin;
dss.odmin = ss.odmin ;
dss.Od = ss.Od ;
dss.Ot = ss.Ot;
%dss.S = ss.S ;
dss.St = ss.St;
dss.Sd = ss.Sd;
dss.sigmat = ss.sigmat ;
dss.sigmad = ss.sigmad ;


for ODepd = 1 : size(ss.octave,2)
    OTime = ODepd;
    [M, N, Time, Depd] = size(ss.octave{OTime, ODepd});
    Vector = zeros(M, N, Depd-1); % difference
    for i = 1: (Depd-1)
        Vector(:, :, i, i) = ss.octave{OTime, ODepd}( :, :, i+1, i+1)-ss.octave{OTime, ODepd}( : ,:, i, i);
    end
    dss.octave{OTime, ODepd} = Vector; % depd difference
    clear Vector vector
end