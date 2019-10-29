%% Stats for EMT Gait Analysis
%Lara Weed
%02 FEB 2019
%% Load data
    load('C:\Users\Laraw\Documents\UVM\Research\McGinnis\Rescue_Climb\stepT20190202.mat')
    subject = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14','S15','S16','S17','S18','S19','S20','S21','S22','S23','S24','S25'}; %subjects

%% Create indices
    A=stepT.Type_Activity=='A'; %ascent
    D=stepT.Type_Activity=='D'; %descent
    Landing=stepT.Type_Activity=='L'; %Landing
    W=stepT.Type_Activity=='W'; %Landing
    R=stepT.fside=='R'; %Right leg
    L=stepT.fside=='L'; %Left Leg
    HB=sum(stepT.bag_type=='hb',2)==2; %Left Leg
    LB=sum(stepT.bag_type=='lb',2)==2; %Left Leg
    HS=sum(stepT.bag_type=='hs',2)==2; %Left Leg
    LS=sum(stepT.bag_type=='ls',2)==2; %Left Leg
    
%% Exclusion criteria:
        %steps longer than 2.25m and shoter than 0.25 m 
            sl_ind = stepT.step_length>=.25 &stepT.step_length<=2.25;
            
        %steps with heights greater than .75m
            sh_ind = stepT.step_height<=.75 ;
            
        %steps with contact time greater than 2 seconds
            ct_ind = stepT.contact_time<=2 ;
            
        %cadence greater than 20 and less than 70
            cad_ind = stepT.cadence>=20 &stepT.cadence<=70;
            
        %lateral devistion less than 0.75m
            ld_ind = stepT.lateral_dev<=.75;
        
        %combine    
            exclusion_ind = sl_ind & sh_ind & ct_ind & cad_ind & ld_ind;     
    
%% Trial Indices
    %WT
     wt_hb =  W & HB & exclusion_ind;
     wt_lb =  W & LB & exclusion_ind;
     wt_hs =  W & HS & exclusion_ind;
     wt_ls =  W & LS & exclusion_ind;
    
     %RC-Ascent 
     rca_hb =  A & HB & exclusion_ind;
     rca_lb =  A & LB & exclusion_ind;
     rca_hs =  A & HS & exclusion_ind;
     rca_ls =  A & LS & exclusion_ind;
     
     %RC-Descent
     rcd_hb =  D & HB & exclusion_ind;
     rcd_lb =  D & LB & exclusion_ind;
     rcd_hs =  D & HS & exclusion_ind;
     rcd_ls =  D & LS & exclusion_ind;
     
     ind_names={'wt_hb','wt_lb','wt_hs','wt_ls','rca_hb','rca_lb','rca_hs','rca_ls','rcd_hb','rcd_lb','rcd_hs','rcd_ls'};
%% Select 10 good steps from each subject for analysis for each activity type  
stepsR = zeros(23,12);
stepsL = zeros(23,12);
steps = zeros(23,12);
trials_ind = [];
for i=[1,2,4:23]
    sub_ind = strcmp(stepT.Subject,subject(:,i));
    for j = 1:12
        t_ind = eval(ind_names{j});
        % count number of steps
        stepsR(i,j) = sum(sub_ind & t_ind& R);
        stepsL(i,j) = sum(sub_ind & t_ind& L);
        steps(i,j) = sum(sub_ind & t_ind);
        
        % Try to get even number of steps per trial
        if j <5
            % 14 steps per trial for wt
            goodSteps = find(sub_ind & t_ind,16);
        elseif j>=5 && j<9
            % 6 per trial for rca
            goodSteps = find(sub_ind & t_ind,6);
        else
            % 3 steps per trial for rcd
            goodSteps = find(sub_ind & t_ind,3);
        end
        trials_ind = [trials_ind;goodSteps,j*ones(length(goodSteps),1)];% index into step table, which trial   
    end
end
%%
vars = stepT.Properties.VariableNames;
for z = 7:14
    clear p1 p2 p3 p4
    for p = 1:12
        concidered_steps_ind = trials_ind((trials_ind(:,2)==p),1);
        if p==1 || p==5 || p==9
            figure;
        end
        if p==1 || p==5 || p==9
            p1 = subplot(1,4,1);
        elseif p==2 || p==6 || p==10
            p2 = subplot(1,4,2);
        elseif p==3 || p==7 || p==11
            p3 = subplot(1,4,3);
        else
            p4 = subplot(1,4,4);
        end
        boxplot(table2array(stepT(concidered_steps_ind,z)))
        xlabel(ind_names{p})
        if p==1 || p==5 || p==9
            ylabel(vars{z});
        end
        hold on
        if exist('p1','var') && exist('p2','var') && exist('p3','var') && exist('p4','var')
            linkaxes([p1 p2 p3 p4],'y')
        end
    end
    
    
end
%%
vars = stepT.Properties.VariableNames;
for z = 7:14
    for p = [1,5,9]
        concidered_steps_ind_hb = trials_ind((trials_ind(:,2)==p),1);
        concidered_steps_ind_lb = trials_ind((trials_ind(:,2)==p+1),1);
        concidered_steps_ind_hs = trials_ind((trials_ind(:,2)==p+2),1);
        concidered_steps_ind_ls = trials_ind((trials_ind(:,2)==p+3),1);
        
       9;
        Meas = table([1 2 3 4 5 6 7 8]','VariableNames',{'Measurements'});
        rm = fitrm(T,'step_length-cadence~Trials','WithinDesign',Meas);
        [ranovatbl,A,C,D] = ranova(rm); 
    end
    linkaxes([p1 p2 p3 p4],'y')
end

%%
vars = stepT.Properties.VariableNames;
rmanova_results = [];
for z = 7:14
    for p = [1,5,9]
        concidered_steps_ind_hb = trials_ind((trials_ind(:,2)==p),1);
        concidered_steps_ind_lb = trials_ind((trials_ind(:,2)==p+1),1);
        concidered_steps_ind_hs = trials_ind((trials_ind(:,2)==p+2),1);
        concidered_steps_ind_ls = trials_ind((trials_ind(:,2)==p+3),1);
        
        val_hb = stepT(concidered_steps_ind_hb,[1,z]);
        val_hb.Properties.VariableNames = {'Subject','HB'};
        val_lb = stepT(concidered_steps_ind_lb,z);
        val_lb.Properties.VariableNames = {'LB'};
        val_hs = stepT(concidered_steps_ind_hs,z);
        val_hs.Properties.VariableNames = {'HS'};
        val_ls = stepT(concidered_steps_ind_ls,z);
        val_ls.Properties.VariableNames = {'LS'};
        
        T = [val_hb, val_lb, val_hs, val_ls];
        Meas = table([1 2 3 4]','VariableNames',{'Measurements'});
        rm = fitrm(T,'HB-LS~Subject','WithinDesign',Meas);
        [ranovatbl,A,C,D] = ranova(rm); 
        
        rmanova_results = [rmanova_results;z,p,rm.ranova.DF(2),rm.ranova.DF(3),rm.ranova.F(2),rm.ranova.pValue(2)];
        
    end
end

%% WT    
    %step length
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(W & HB,:).step_length)
    title('Step length WT')
    p2 = subplot(1,4,2);
    boxplot(stepT(W & LB,:).step_length)
    p3 = subplot(1,4,3);
    boxplot(stepT(W & HS,:).step_length)
    p4 = subplot(1,4,4);
    boxplot(stepT(W & LS,:).step_length)
    linkaxes([p1 p2 p3 p4],'y')
    
    %step length
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(W & HB,:).step_height)
    title('Step Height WT')
    p2 = subplot(1,4,2);
    boxplot(stepT(W & LB,:).step_height)
    p3 = subplot(1,4,3);
    boxplot(stepT(W & HS,:).step_height)
    p4 = subplot(1,4,4);
    boxplot(stepT(W & LS,:).step_height)
    linkaxes([p1 p2 p3 p4],'y')
    
    %Lateral Dev
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(W & HB,:).lateral_dev)
    title('Lateral Deviation WT')
    p2 = subplot(1,4,2);
    boxplot(stepT(W & LB,:).lateral_dev)
    p3 = subplot(1,4,3);
    boxplot(stepT(W & HS,:).lateral_dev)
    p4 = subplot(1,4,4);
    boxplot(stepT(W & LS,:).lateral_dev)
    linkaxes([p1 p2 p3 p4],'y')
    
    
    %Max Swing Velocity
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(W & HB,:).max_swing_vel)
    title('Max Swing Velocity WT')
    p2 = subplot(1,4,2);
    boxplot(stepT(W & LB,:).max_swing_vel)
    p3 = subplot(1,4,3);
    boxplot(stepT(W & HS,:).max_swing_vel)
    p4 = subplot(1,4,4);
    boxplot(stepT(W & LS,:).max_swing_vel)
    linkaxes([p1 p2 p3 p4],'y')
    
    %Foot Attack Angle
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(W & HB,:).foot_attack_angle)
    title('Foot Attack Angle WT')
    p2 = subplot(1,4,2);
    boxplot(stepT(W & LB,:).foot_attack_angle)
    p3 = subplot(1,4,3);
    boxplot(stepT(W & HS,:).foot_attack_angle)
    p4 = subplot(1,4,4);
    boxplot(stepT(W & LS,:).foot_attack_angle)
    linkaxes([p1 p2 p3 p4],'y')
    
    %contact time
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(W & HB,:).contact_time)
    title('Contact Time WT')
    p2 = subplot(1,4,2);
    boxplot(stepT(W & LB,:).contact_time)
    p3 = subplot(1,4,3);
    boxplot(stepT(W & HS,:).contact_time)
    p4 = subplot(1,4,4);
    boxplot(stepT(W & LS,:).contact_time)
    linkaxes([p1 p2 p3 p4],'y')
    
    %Step time
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(W & HB,:).step_time)
    title('Step Time WT')
    p2 = subplot(1,4,2);
    boxplot(stepT(W & LB,:).step_time)
    p3 = subplot(1,4,3);
    boxplot(stepT(W & HS,:).step_time)
    p4 = subplot(1,4,4);
    boxplot(stepT(W & LS,:).step_time)
    linkaxes([p1 p2 p3 p4],'y')
    
    
    %Step time
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(W & HB,:).cadence)
    title('Cadence WT')
    p2 = subplot(1,4,2);
    boxplot(stepT(W & LB,:).cadence)
    p3 = subplot(1,4,3);
    boxplot(stepT(W & HS,:).cadence)
    p4 = subplot(1,4,4);
    boxplot(stepT(W & LS,:).cadence)
    linkaxes([p1 p2 p3 p4],'y')
    
 %% RCA   
    %step length
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(A & HB,:).step_length)
    title('Step length RCA')
    p2 = subplot(1,4,2);
    boxplot(stepT(A & LB,:).step_length)
    p3 = subplot(1,4,3);
    boxplot(stepT(A & HS,:).step_length)
    p4 = subplot(1,4,4);
    boxplot(stepT(A & LS,:).step_length)
    linkaxes([p1 p2 p3 p4],'y')
    
    %step length
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(A & HB,:).step_height)
    title('Step Height RCA')
    p2 = subplot(1,4,2);
    boxplot(stepT(A & LB,:).step_height)
    p3 = subplot(1,4,3);
    boxplot(stepT(A & HS,:).step_height)
    p4 = subplot(1,4,4);
    boxplot(stepT(A & LS,:).step_height)
    linkaxes([p1 p2 p3 p4],'y')
    
    %Lateral Dev
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(A & HB,:).lateral_dev)
    title('Lateral Deviation RCA')
    p2 = subplot(1,4,2);
    boxplot(stepT(A & LB,:).lateral_dev)
    p3 = subplot(1,4,3);
    boxplot(stepT(A & HS,:).lateral_dev)
    p4 = subplot(1,4,4);
    boxplot(stepT(A & LS,:).lateral_dev)
    linkaxes([p1 p2 p3 p4],'y')
    
    
    %Max Swing Velocity
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(A & HB,:).max_swing_vel)
    title('Max Swing Velocity RCA')
    p2 = subplot(1,4,2);
    boxplot(stepT(A & LB,:).max_swing_vel)
    p3 = subplot(1,4,3);
    boxplot(stepT(A & HS,:).max_swing_vel)
    p4 = subplot(1,4,4);
    boxplot(stepT(A & LS,:).max_swing_vel)
    linkaxes([p1 p2 p3 p4],'y')
    
    %Foot Attack Angle
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(A & HB,:).foot_attack_angle)
    title('Foot Attack Angle RCA')
    p2 = subplot(1,4,2);
    boxplot(stepT(A & LB,:).foot_attack_angle)
    p3 = subplot(1,4,3);
    boxplot(stepT(A & HS,:).foot_attack_angle)
    p4 = subplot(1,4,4);
    boxplot(stepT(A & LS,:).foot_attack_angle)
    linkaxes([p1 p2 p3 p4],'y')
    
    %contact time
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(A & HB,:).contact_time)
    title('Contact Time RCA')
    p2 = subplot(1,4,2);
    boxplot(stepT(A & LB,:).contact_time)
    p3 = subplot(1,4,3);
    boxplot(stepT(A & HS,:).contact_time)
    p4 = subplot(1,4,4);
    boxplot(stepT(A & LS,:).contact_time)
    linkaxes([p1 p2 p3 p4],'y')
    
    %Step time
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(A & HB,:).step_time)
    title('Step Time RCA')
    p2 = subplot(1,4,2);
    boxplot(stepT(A & LB,:).step_time)
    p3 = subplot(1,4,3);
    boxplot(stepT(A & HS,:).step_time)
    p4 = subplot(1,4,4);
    boxplot(stepT(A & LS,:).step_time)
    linkaxes([p1 p2 p3 p4],'y')
    
    
    %Step time
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(A & HB,:).cadence)
    title('Cadence RCA')
    p2 = subplot(1,4,2);
    boxplot(stepT(A & LB,:).cadence)
    p3 = subplot(1,4,3);
    boxplot(stepT(A & HS,:).cadence)
    p4 = subplot(1,4,4);
    boxplot(stepT(A & LS,:).cadence)
    linkaxes([p1 p2 p3 p4],'y')
    
    
%% RCD   
    %step length
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(D & HB,:).step_length)
    title('Step length RCD')
    p2 = subplot(1,4,2);
    boxplot(stepT(D & LB,:).step_length)
    p3 = subplot(1,4,3);
    boxplot(stepT(D & HS,:).step_length)
    p4 = subplot(1,4,4);
    boxplot(stepT(D & LS,:).step_length)
    linkaxes([p1 p2 p3 p4],'y')
    
    %step length
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(D & HB,:).step_height)
    title('Step Height RCD')
    p2 = subplot(1,4,2);
    boxplot(stepT(D & LB,:).step_height)
    p3 = subplot(1,4,3);
    boxplot(stepT(D & HS,:).step_height)
    p4 = subplot(1,4,4);
    boxplot(stepT(D & LS,:).step_height)
    linkaxes([p1 p2 p3 p4],'y')
    
    %Lateral Dev
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(D & HB,:).lateral_dev)
    title('Lateral Deviation RCD')
    p2 = subplot(1,4,2);
    boxplot(stepT(D & LB,:).lateral_dev)
    p3 = subplot(1,4,3);
    boxplot(stepT(D & HS,:).lateral_dev)
    p4 = subplot(1,4,4);
    boxplot(stepT(D & LS,:).lateral_dev)
    linkaxes([p1 p2 p3 p4],'y')
    
    
    %Max Swing Velocity
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(D & HB,:).max_swing_vel)
    title('Max Swing Velocity RCD')
    p2 = subplot(1,4,2);
    boxplot(stepT(D & LB,:).max_swing_vel)
    p3 = subplot(1,4,3);
    boxplot(stepT(D & HS,:).max_swing_vel)
    p4 = subplot(1,4,4);
    boxplot(stepT(D & LS,:).max_swing_vel)
    linkaxes([p1 p2 p3 p4],'y')
    
    %Foot Attack Angle
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(D & HB,:).foot_attack_angle)
    title('Foot Attack Angle RCD')
    p2 = subplot(1,4,2);
    boxplot(stepT(D & LB,:).foot_attack_angle)
    p3 = subplot(1,4,3);
    boxplot(stepT(D & HS,:).foot_attack_angle)
    p4 = subplot(1,4,4);
    boxplot(stepT(D & LS,:).foot_attack_angle)
    linkaxes([p1 p2 p3 p4],'y')
    
    %contact time
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(D & HB,:).contact_time)
    title('Contact Time RCD')
    p2 = subplot(1,4,2);
    boxplot(stepT(D & LB,:).contact_time)
    p3 = subplot(1,4,3);
    boxplot(stepT(D & HS,:).contact_time)
    p4 = subplot(1,4,4);
    boxplot(stepT(D & LS,:).contact_time)
    linkaxes([p1 p2 p3 p4],'y')
    
    %Step time
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(D & HB,:).step_time)
    title('Step Time RCD')
    p2 = subplot(1,4,2);
    boxplot(stepT(D & LB,:).step_time)
    p3 = subplot(1,4,3);
    boxplot(stepT(D & HS,:).step_time)
    p4 = subplot(1,4,4);
    boxplot(stepT(D & LS,:).step_time)
    linkaxes([p1 p2 p3 p4],'y')
    
    
    %Step time
    figure;
    p1 = subplot(1,4,1);
    boxplot(stepT(D & HB,:).cadence)
    title('Cadence RCD')
    p2 = subplot(1,4,2);
    boxplot(stepT(D & LB,:).cadence)
    p3 = subplot(1,4,3);
    boxplot(stepT(D & HS,:).cadence)
    p4 = subplot(1,4,4);
    boxplot(stepT(D & LS,:).cadence)
    linkaxes([p1 p2 p3 p4],'y')
    
        
    
    
    
    
    
    
