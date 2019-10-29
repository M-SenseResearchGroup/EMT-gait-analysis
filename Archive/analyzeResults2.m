%% Load Data
load('C:\Users\Laraw\Documents\UVM\Research\McGinnis\Rescue_Climb\DataProcessed\stair_table_20181128.mat')


%% Create indices
A=stepT.Type_stair=='A';
D=stepT.Type_stair=='D';
R=stepT.fside=='R';
L=stepT.fside=='L';

%% Assign column names
stepT.Properties.VariableNames = {'sub','foot','stair','bag','startt','endt','sl','ld','sh','msv','faa','ct','st','cad'};
rmAnovaTinput = stepT(A,[4,7:end]);
Meas = table([1:8]','VariableNames',{'Measurements'});
rm = fitrm(rmAnovaTinput,'sl-cad~bag','WithinDesign',Meas);
ranovatbl = ranova(rm);
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
mean_A=[mean_A; mean(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(A))];
mean_right =[mean_right; mean(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(A & R))];
mean_left =[mean_left; mean(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(A & L))];
std_A =[std_A;std(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(A))];
std_right =[std_right;std(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(A & R))];
std_left = [std_left;std(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(A & L))];
skew_A = [skew_A;skewness(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(A))];
skew_right = [skew_right;skewness(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(A & R))];
skew_left = [skew_left;skewness(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(A & L))];
kurt_A = [kurt_A; kurtosis(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(A))];
kurt_right = [kurt_right; kurtosis(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(A & R))];
kurt_left = [kurt_left; kurtosis(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(A & L))];

mean_D=[mean_D; mean(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(D))];
mean_rightD =[mean_rightD; mean(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(D & R))];
mean_leftD =[mean_leftD; mean(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(D & L))];
std_D =[std_D;std(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(D))];
std_rightD =[std_rightD;std(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(D & R))];
std_leftD = [std_leftD;std(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(D & L))];
skew_D = [skew_D;skewness(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(D))];
skew_rightD = [skew_rightD;skewness(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(D & L))];
skew_leftD = [skew_leftD;skewness(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(D & L))];
kurt_D = [kurt_D; kurtosis(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(D))];
kurt_rightD = [kurt_rightD; kurtosis(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(D & R))];
kurt_leftD = [kurt_leftD; kurtosis(stepT.(cell2mat(stepT.Properties.VariableNames(i)))(D & L))];
end
Atab=table(mean_A, mean_right, mean_left, std_A, std_right, std_left, skew_A, skew_right, skew_left, kurt_A, kurt_right, kurt_left); 
Dtab=table(mean_D, mean_rightD, mean_leftD, std_D, std_rightD, std_leftD, skew_D, skew_rightD, skew_leftD, kurt_D, kurt_rightD, kurt_leftD); 

%% Step length
figure;
normplot(stepT.sl(A))
[hsla,psla] = lillietest(stepT.sl(A));

figure('name',' Ascent Step_length');
boxplot(stepT.sl(A),stepT.bag(A))
ylabel('Step Length (m)');
xlabel('Bag Type');
title('Ascent')

figure;
normplot(stepT.sl(D))
[hsld,psld] = lillietest(stepT.sl(D));

figure('name','Descent Step_length');
boxplot(stepT.sl(D),stepT.bag(D))
ylabel('Step Length (m)');
xlabel('Bag Type');
title('Descent')

%% Lateral Deviation
figure;
normplot(stepT.ld(A))
[hlda,plda] = lillietest(stepT.ld(A));

figure('name',' Ascent Lateral Deviation');
boxplot(stepT.ld(A),stepT.bag(A))
ylabel('Lateral Deviation (m)');
xlabel('Bag Type');
title('Ascent')

figure;
normplot(stepT.ld(D))
[hldd,pldd] = lillietest(stepT.ld(D));

figure('name','Descent Lateral Deviation');
boxplot(stepT.ld(D),stepT.bag(D))
ylabel('Lateral Deviation (m)');
xlabel('Bag Type');
title('Descent')

%% Step Height
figure;
normplot(stepT.sh(A))
[hsha,psha] = lillietest(stepT.sh(A));

figure('name',' Ascent Step Height');
boxplot(stepT.sh(A),stepT.bag(A))
ylabel('Step Height (m)');
xlabel('Bag Type');
title('Ascent')

figure;
normplot(stepT.sh(D))
[hshd,pshd] = lillietest(stepT.sh);
figure('name','Descent Step Height');
boxplot(stepT.sh(D),stepT.bag(D))
ylabel('Step Height (m)');
xlabel('Bag Type');
title('Descent')

%% Max Swing Velocity
figure;
normplot(stepT.msv(A))
[hmsva,pmsva] = lillietest(stepT.msv(A));

figure('name',' Ascent Max Swing Velocity');
boxplot(stepT.msv(A),stepT.bag(A))
ylabel('Max Swing Velocity (m/s)');
xlabel('Bag Type');
title('Ascent')

figure;
normplot(stepT.msv(D))
[hmsvd,pmsvd] = lillietest(stepT.msv(D));

figure('name','Descent Max Swing Velocity');
boxplot(stepT.msv(D),stepT.bag(D))
ylabel('Lateral Max Swing Velocity (m/s)');
xlabel('Bag Type');
title('Descent')

%% Foot Attack Angle
figure;
normplot(stepT.faa(A))
[hfaaa,pfaaa] = lillietest(stepT.faa(A));

figure('name',' Ascent Foot Attack Angle');
boxplot(stepT.faa(A),stepT.bag(A))
ylabel('Foot Attack Angle (deg/s)');
xlabel('Bag Type');
title('Ascent')

figure;
normplot(stepT.faa(D))
[hfaad,pfaad] = lillietest(stepT.faa(D));

figure('name','Descent Foot Attack Angle');
boxplot(stepT.faa(D),stepT.bag(D))
ylabel('Foot Attack Angle (deg/s)');
xlabel('Bag Type');
title('Descent')

%% Contact Time
figure;
normplot(stepT.ct(A))
[hcta,pcta] = lillietest(stepT.ct(A));

figure('name',' Ascent Contact Time');
boxplot(stepT.ct(A),stepT.bag(A))
ylabel('Contact Time (s)');
xlabel('Bag Type');
title('Ascent')

figure;
normplot(stepT.ct(D))
[hctd,pctd] = lillietest(stepT.ct(D));

figure('name','Descent Foot Attack Angle');
boxplot(stepT.ct(D),stepT.bag(D))
ylabel('Contact Time (s)');
xlabel('Bag Type');
title('Descent')

%% Step Time
figure;
normplot(stepT.st(A))
[hsta,psta] = lillietest(stepT.st(A));

figure('name',' Ascent Step_Time');
boxplot(stepT.st(A),stepT.bag(A))
ylabel('Step Time (s)');
xlabel('Bag Type');
title('Ascent')

figure;
normplot(stepT.st(D))
[hstd,pstd] = lillietest(stepT.st(D));

figure('name','Descent Step_Time');
boxplot(stepT.st(D),stepT.bag(D))
ylabel('Step Time (s)');
xlabel('Bag Type');
title('Descent')

%% Cadence
figure;
normplot(stepT.cad(A))
[hcda,pcda] = lillietest(stepT.cad(A));

figure('name',' Ascent Cadence');
boxplot(stepT.cad(A),stepT.bag(A))
ylabel('Cadence (Steps/s)');
xlabel('Bag Type');
title('Ascent')

figure;
normplot(stepT.cad(D))
[hcdd,pcdd] = lillietest(stepT.cad(D));

figure('name','Descent Cadence');
boxplot(stepT.cad(D),stepT.bag(D))
ylabel('Cadence (Steps/s)');
xlabel('Bag Type');
title('Descent')
