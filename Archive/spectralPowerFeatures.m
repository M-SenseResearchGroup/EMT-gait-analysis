function feats = spectralPowerFeatures(x, fs)
 
feats = zeros(5,1);
 
edges = [0.5, 1.5, 5, 10, 15, 20];
 
[p, f] = periodogram(x,[],4096,fs);
    
for kband = 1:length(edges)-1
    feats(kband) = sum(p( (f>=edges(kband)) & (f<edges(kband+1)) ));
end
end
