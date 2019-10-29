%% Load data
load('Results.mat')

%% Assign column names
Results.Properties.VariableNames = {'sub','foot','bag','dir','sl','ld','sh','msv','faa','ct','st','cad'};

%% Create indices
A=Results.dir==1;
D=Results.dir==2;
AR=(Results.dir==1 & Results.foot==1);
AL=(Results.dir==1 & Results.foot==2);
DR=(Results.dir==2 & Results.foot==1);
DL=(Results.dir==2 & Results.foot==2);

%% Descriptive Statistics
mean_A=[];
mean_right=[];
mean_left=[];
std_A=[];
std_right=[];
std_left=[];
skew_A=[];
skew_right=[];
skew_left=[];
kurt_A=[];
kurt_right=[];
kurt_left=[];
mean_D=[];
mean_rightD=[];
mean_leftD=[];
std_D=[];
std_rightD=[];
std_leftD=[];
skew_D=[];
skew_rightD=[];
skew_leftD=[];
kurt_D=[];
kurt_rightD=[];
kurt_leftD=[];
for i=5:12
%step Length
mean_A=[mean_A; mean(Results.(cell2mat(Results.Properties.VariableNames(i)))(A))];
mean_right =[mean_right; mean(Results.(cell2mat(Results.Properties.VariableNames(i)))(AR))];
mean_left =[mean_left; mean(Results.(cell2mat(Results.Properties.VariableNames(i)))(AL))];
std_A =[std_A;std(Results.(cell2mat(Results.Properties.VariableNames(i)))(A))];
std_right =[std_right;std(Results.(cell2mat(Results.Properties.VariableNames(i)))(AR))];
std_left = [std_left;std(Results.(cell2mat(Results.Properties.VariableNames(i)))(AL))];
skew_A = [skew_A;skewness(Results.(cell2mat(Results.Properties.VariableNames(i)))(A))];
skew_right = [skew_right;skewness(Results.(cell2mat(Results.Properties.VariableNames(i)))(AR))];
skew_left = [skew_left;skewness(Results.(cell2mat(Results.Properties.VariableNames(i)))(AL))];
kurt_A = [kurt_A; kurtosis(Results.(cell2mat(Results.Properties.VariableNames(i)))(A))];
kurt_right = [kurt_right; kurtosis(Results.(cell2mat(Results.Properties.VariableNames(i)))(AR))];
kurt_left = [kurt_left; kurtosis(Results.(cell2mat(Results.Properties.VariableNames(i)))(AL))];

mean_D=[mean_D; mean(Results.(cell2mat(Results.Properties.VariableNames(i)))(D))];
mean_rightD =[mean_rightD; mean(Results.(cell2mat(Results.Properties.VariableNames(i)))(DR))];
mean_leftD =[mean_leftD; mean(Results.(cell2mat(Results.Properties.VariableNames(i)))(DL))];
std_D =[std_D;std(Results.(cell2mat(Results.Properties.VariableNames(i)))(D))];
std_rightD =[std_rightD;std(Results.(cell2mat(Results.Properties.VariableNames(i)))(DR))];
std_leftD = [std_leftD;std(Results.(cell2mat(Results.Properties.VariableNames(i)))(DL))];
skew_D = [skew_D;skewness(Results.(cell2mat(Results.Properties.VariableNames(i)))(D))];
skew_rightD = [skew_rightD;skewness(Results.(cell2mat(Results.Properties.VariableNames(i)))(DL))];
skew_leftD = [skew_leftD;skewness(Results.(cell2mat(Results.Properties.VariableNames(i)))(DL))];
kurt_D = [kurt_D; kurtosis(Results.(cell2mat(Results.Properties.VariableNames(i)))(D))];
kurt_rightD = [kurt_rightD; kurtosis(Results.(cell2mat(Results.Properties.VariableNames(i)))(DR))];
kurt_leftD = [kurt_leftD; kurtosis(Results.(cell2mat(Results.Properties.VariableNames(i)))(DL))];
end
Atab=table(mean_A, mean_right, mean_left, std_A, std_right, std_left, skew_A, skew_right, skew_left, kurt_A, kurt_right, kurt_left); 
Dtab=table(mean_D, mean_rightD, mean_leftD, std_D, std_rightD, std_leftD, skew_D, skew_rightD, skew_leftD, kurt_D, kurt_rightD, kurt_leftD); 

%% Lateral Dev
[R.a.ld,~] = boxcox(Results.ld(ascent_ind));
[R.d.ld,~] = boxcox(Results.ld(descent_ind));
%Step Height
[R.a.sh,~] = boxcox(Results.sh(ascent_ind)); 
[R.d.sh,~] = boxcox(Results.sh(descent_ind));
% max swing vel
[R.a.msv,~] = boxcox(Results.msv(ascent_ind));
[R.d.msv,~] = boxcox(Results.msv(descent_ind));
% contact time
[R.a.ct,~] = boxcox(Results.ct(ascent_ind));
[R.d.ct,~] = boxcox(Results.ct(descent_ind));
%step time
[R.a.st,~] = boxcox(Results.st(ascent_ind));
[R.d.st,~] = boxcox(Results.st(descent_ind));
%cadence
[R.a.cad,~] = boxcox(Results.cad(ascent_ind));
[R.d.cad,~] = boxcox(Results.cad(descent_ind));
% foot attack angle
%% Make Data Normal
%step Length
[R.a.sl,~] = boxcox(Results.sl(ascent_ind));
[R.d.sl,~] = boxcox(Results.sl(descent_ind));
% Lateral Dev
[R.a.ld,~] = boxcox(Results.ld(ascent_ind));
[R.d.ld,~] = boxcox(Results.ld(descent_ind));
%Step Height
[R.a.sh,~] = boxcox(Results.sh(ascent_ind)); 
[R.d.sh,~] = boxcox(Results.sh(descent_ind));
% max swing vel
[R.a.msv,~] = boxcox(Results.msv(ascent_ind));
[R.d.msv,~] = boxcox(Results.msv(descent_ind));
% contact time
[R.a.ct,~] = boxcox(Results.ct(ascent_ind));
[R.d.ct,~] = boxcox(Results.ct(descent_ind));
%step time
[R.a.st,~] = boxcox(Results.st(ascent_ind));
[R.d.st,~] = boxcox(Results.st(descent_ind));
%cadence
[R.a.cad,~] = boxcox(Results.cad(ascent_ind));
[R.d.cad,~] = boxcox(Results.cad(descent_ind));
%% foot attack angle
mffaa=min(Results.faa(ascent_ind))-1;
mffad=min(Results.faa(descent_ind))-1;
[R.a.faa,~] = boxcox(Results.faa(ascent_ind)- mffaa);
[R.d.faa,~] = boxcox(Results.faa(descent_ind)- mffad);
%% Step Length
figure;
normplot(Results.sl(ascent_ind))
[hsla,psla] = lillietest(Results.sl(ascent_ind));
[transdat,lambda] = boxcox(Results.sl(ascent_ind));
figure;
normplot(transdat)
[hslabc,pslabc] = lillietest(transdat);

%%



figure;
normplot(Results.sl(ascent_ind))
[hsla,psla] = lillietest(Results.sl(ascent_ind));

figure('name',' Ascent Step_length');
boxplot(Results.sl(ascent_ind),Results.bag(ascent_ind))
ylabel('Step Length (m)');
xlabel('Bag Type');
title('Ascent')

figure;
normplot(Results.sl(descent_ind))
[hsld,psld] = lillietest(Results.sl(descent_ind));

figure('name','Descent Step_length');
boxplot(Results.sl(descent_ind),Results.bag(descent_ind))
ylabel('Step Length (m)');
xlabel('Bag Type');
title('Descent')

%% Lateral Deviation
figure;
normplot(Results.lateral_dev(ascent_ind))
[hlda,plda] = lillietest(Results.lateral_dev(ascent_ind));

figure('name',' Ascent Lateral Deviation');
boxplot(Results.lateral_dev(ascent_ind),Results.Bag_Type(ascent_ind))
ylabel('Lateral Deviation (m)');
xlabel('Bag Type');
title('Ascent')

figure;
normplot(Results.lateral_dev(descent_ind))
[hldd,pldd] = lillietest(Results.lateral_dev(descent_ind));

figure('name','Descent Lateral Deviation');
boxplot(Results.lateral_dev(descent_ind),Results.Bag_Type(descent_ind))
ylabel('Lateral Deviation (m)');
xlabel('Bag Type');
title('Descent')

%% Step Height
figure;
normplot(Results.step_height(ascent_ind))
[hsha,psha] = lillietest(Results.step_height(ascent_ind));

figure('name',' Ascent Step Height');
boxplot(Results.step_height(ascent_ind),Results.Bag_Type(ascent_ind))
ylabel('Step Height (m)');
xlabel('Bag Type');
title('Ascent')

figure;
normplot(Results.step_height(descent_ind))
[hshd,pshd] = lillietest(Results.step_height(descent_ind));

figure('name','Descent Step Height');
boxplot(Results.step_height(descent_ind),Results.Bag_Type(descent_ind))
ylabel('Step Height (m)');
xlabel('Bag Type');
title('Descent')

%% Max Swing Velocity
figure;
normplot(Results.max_swing_vel(ascent_ind))
[hmsva,pmsva] = lillietest(Results.max_swing_vel(ascent_ind));

figure('name',' Ascent Max Swing Velocity');
boxplot(Results.max_swing_vel(ascent_ind),Results.Bag_Type(ascent_ind))
ylabel('Max Swing Velocity (m/s)');
xlabel('Bag Type');
title('Ascent')

figure;
normplot(Results.max_swing_vel(descent_ind))
[hmsvd,pmsvd] = lillietest(Results.max_swing_vel(descent_ind));

figure('name','Descent Max Swing Velocity');
boxplot(Results.max_swing_vel(descent_ind),Results.Bag_Type(descent_ind))
ylabel('Lateral Max Swing Velocity (m/s)');
xlabel('Bag Type');
title('Descent')

%% Foot Attack Angle
figure;
normplot(Results.foot_attack_angle(ascent_ind))
[hfaaa,pfaaa] = lillietest(Results.foot_attack_angle(ascent_ind));

figure('name',' Ascent Foot Attack Angle');
boxplot(Results.foot_attack_angle(ascent_ind),Results.Bag_Type(ascent_ind))
ylabel('Foot Attack Angle (deg/s)');
xlabel('Bag Type');
title('Ascent')

figure;
normplot(Results.foot_attack_angle(descent_ind))
[hfaad,pfaad] = lillietest(Results.foot_attack_angle(descent_ind));

figure('name','Descent Foot Attack Angle');
boxplot(Results.foot_attack_angle(descent_ind),Results.Bag_Type(descent_ind))
ylabel('Foot Attack Angle (deg/s)');
xlabel('Bag Type');
title('Descent')

%% Contact Time
figure;
normplot(Results.contact_time(ascent_ind))
[hcta,pcta] = lillietest(Results.contact_time(ascent_ind));

figure('name',' Ascent Contact Time');
boxplot(Results.contact_time(ascent_ind),Results.Bag_Type(ascent_ind))
ylabel('Contact Time (s)');
xlabel('Bag Type');
title('Ascent')

figure;
normplot(Results.contact_time(descent_ind))
[hctd,pctd] = lillietest(Results.contact_time(descent_ind));

figure('name','Descent Foot Attack Angle');
boxplot(Results.contact_time(descent_ind),Results.Bag_Type(descent_ind))
ylabel('Contact Time (s)');
xlabel('Bag Type');
title('Descent')

%% Step Time
figure;
normplot(Results.step_time(ascent_ind))
[hsta,psta] = lillietest(Results.step_time(ascent_ind));

figure('name',' Ascent Step_Time');
boxplot(Results.step_time(ascent_ind),Results.Bag_Type(ascent_ind))
ylabel('Step Time (s)');
xlabel('Bag Type');
title('Ascent')

figure;
normplot(Results.step_time(descent_ind))
[hstd,pstd] = lillietest(Results.step_time(descent_ind));

figure('name','Descent Step_Time');
boxplot(Results.step_time(descent_ind),Results.Bag_Type(descent_ind))
ylabel('Step Time (s)');
xlabel('Bag Type');
title('Descent')

%% Cadence
figure;
normplot(Results.cadence(ascent_ind))
[hcda,pcda] = lillietest(Results.cadence(ascent_ind));

figure('name',' Ascent Cadence');
boxplot(Results.cadence(ascent_ind),Results.Bag_Type(ascent_ind))
ylabel('Cadence (Steps/s)');
xlabel('Bag Type');
title('Ascent')

figure;
normplot(Results.cadence(descent_ind))
[hcdd,pcdd] = lillietest(Results.cadence(descent_ind));

figure('name','Descent Cadence');
boxplot(Results.cadence(descent_ind),Results.Bag_Type(descent_ind))
ylabel('Cadence (Steps/s)');
xlabel('Bag Type');
title('Descent')
