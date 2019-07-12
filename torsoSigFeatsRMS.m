function feat = torsoSigFeatsRMS(x)
    % Initialize feature vector
    feat = nan(7,1);
    try
    % Time domain
    feat(1) = mean(x);
    feat(2) = max(x);
    feat(3) = min(x);
    feat(4) = std(x);
    feat(5) = range(x);
    feat(6) = skewness(x);
    feat(7) = kurtosis(x);
    catch
    end
    
end
 