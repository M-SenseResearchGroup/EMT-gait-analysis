function feats = spectralPeaksFeatures(x, fs)
 
mindist_xunits = 0.3;
 
feats = zeros(1,12);
 
N = 4096;
minpkdist = floor(mindist_xunits/(fs/N));
 
[p, f] = pwelch(x,rectwin(length(x)),[],N,fs);
 
[pks,locs] = findpeaks(p,'npeaks',20,'minpeakdistance',minpkdist);
if(~isempty(pks))
    mx = min(6,length(pks));
    [spks, idx] = sort(pks,'descend');
    slocs = locs(idx);
 
    pks = spks(1:mx);
    locs = slocs(1:mx);
 
    [slocs, idx] = sort(locs,'ascend');
    spks = pks(idx);
    pks = spks;
    locs = slocs;
end
fpk = f(locs);
 
% Features 1-6 positions of highest 6 peaks
feats(1:length(pks)) = fpk;
% Features 7-12 power levels of highest 6 peaks
feats(7:7+length(pks)-1) = pks;
end
 