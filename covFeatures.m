function feats = covFeatures(x, fs)
 
feats = zeros(1,3);
 
[c, lags] = xcorr(x);
 
minprom = 0.0005;
mindist_xunits = 0.3;
minpkdist = floor(mindist_xunits/(1/fs));
[pks,locs] = findpeaks(c,...
    'minpeakprominence',minprom,...
    'minpeakdistance',minpkdist);
 
tc = (1/fs)*lags;
tcl = tc(locs);
% Feature 1 - peak height at 0
if(~isempty(tcl))   % else f1 already 0
    feats(1) = pks((end+1)/2);
end
% Features 2 and 3 - position and height of first peak 
if(length(tcl) >= 3)   % else f2,f3 already 0
    feats(2) = tcl((end+1)/2+1);
    feats(3) = pks((end+1)/2+1);
end
end