function feat = signal_features(x,fs)
    % Initialize feature vector
    feat = zeros(29,1);
    
    % Time domain
    feat(1) = mean(x);
    feat(2) = rms(x);
    feat(3) = skewness(x);
    feat(4) = kurtosis(x);
    feat(5) = range(x);
    feat(6) = max(x);
    feat(7) = min(x);
    feat(8) = std(x);
    feat(9) = peak2rms(x);
    feat(10:12) = covFeatures(x, fs);
    
    % Frequency domain
    feat(13:24) = spectralPeaksFeatures(x, fs);
    feat(25:29) = spectralPowerFeatures(x, fs);
end
 