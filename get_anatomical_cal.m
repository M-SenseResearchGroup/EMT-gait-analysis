function R = get_anatomical_cal(acc_s,g)
%Perform PCA to extract first principal component
coeff = pca(acc_s.','NumComponents',3); %1st principal comp will be body 
                                        %fixed and aligned approx with 
                                        %direction of travel during walking
acc_r = acc_s.' * coeff;
xp = coeff(:,1).';
if mean(acc_r(:,1)-acc_r(1,1)) < 0
    xp = -xp;
end
    
% Construct rotation using direction of gravity while standing and xp
Z = g;
Y = cross(Z,xp) ./ norm(cross(Z,xp)); %approx ML axis of foot
X = cross(Y,Z); %approx AP axis of foot
R = [X; Y; Z];

end

