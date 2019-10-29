%% Stats for EMT Gait Analysis
%Lara Weed
%07 APR 2019
%% Load data
    load('C:\Users\Laraw\OneDrive - University of Vermont\UVM\Research\McGinnis\Rescue_Climb\data\Processed\stepT20190407.mat')
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
for i=[1,2,4:23] % subjects
    sub_ind = strcmp(stepT.Subject,subject(:,i));
    for j = 1:12 % Trails
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
%% Plot Variable for all 12 trials
% vars = stepT.Properties.VariableNames;
% for z = 7:37
%     clear p1 p2 p3 p4
%     for p = 1:12
%         concidered_steps_ind = trials_ind((trials_ind(:,2)==p),1);
%         if p==1 || p==5 || p==9
%             figure;
%         end
%         if p==1 || p==5 || p==9
%             p1 = subplot(1,4,1);
%         elseif p==2 || p==6 || p==10
%             p2 = subplot(1,4,2);
%         elseif p==3 || p==7 || p==11
%             p3 = subplot(1,4,3);
%         else
%             p4 = subplot(1,4,4);
%         end
%         boxplot(table2array(stepT(concidered_steps_ind,z)))
%         xlabel(ind_names{p})
%         if p==1 || p==5 || p==9
%             ylabel(vars{z});
%         end
%         hold on
%         if p==4 || p==8 || p==12            
%             linkaxes([p1 p2 p3 p4],'y')
%         end
%     end    
% end

%% Friedman test and Post-Hoc Mann-Whitney U test
ft_results=[];
ft_results_comp=[];
vars = stepT.Properties.VariableNames;
trial = {'WT','','','','A','','','','D'};
rmanova_results_matrix = [];
phtt_matrix =[];
NT = [];
for z = 7:37
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
        
        %Test normality
%         [h_hb,p_hb] = lillietest(val_hb.HB);
%         [h_lb,p_lb] = lillietest(val_lb.LB);
%         [h_hs,p_hs] = lillietest(val_hs.HS);
%         [h_ls,p_ls] = lillietest(val_ls.LS);
%         NT = [NT;z,p,p_hb,p_lb,p_hs,p_ls];
        
        %Repeated measures anova
        T = [val_hb, val_lb, val_hs, val_ls];
        VarVals.(trial{p}).(vars{z})=T;
        Meas = table([1 2 3 4]','VariableNames',{'Measurements'});
        rm = fitrm(T,'HB-LS~Subject','WithinDesign',Meas);
        [ranovatbl,A,C,D] = ranova(rm); 
        rmanova_results_matrix = [rmanova_results_matrix;z,p,rm.ranova.DF(2),rm.ranova.DF(3),rm.ranova.F(2),rm.ranova.pValue(2)];
        
        % Post Hoc mann-Whitney U test
        [p1,h1,stats1] = ranksum(val_hb.HB,val_lb.LB);
        [p2,h2,stats2] = ranksum(val_hb.HB,val_hs.HS);
        [p3,h3,stats3] = ranksum(val_hb.HB,val_ls.LS);
        [p4,h4,stats4] = ranksum(val_lb.LB,val_hs.HS);
        [p5,h5,stats5] = ranksum(val_lb.LB,val_ls.LS);
        [p6,h6,stats6] = ranksum(val_ls.LS,val_hs.HS);
        phtt_matrix = [phtt_matrix;z,p,p1,p2,p3,p4,p5,p6,stats1.zval,stats2.zval,stats3.zval,stats4.zval,stats5.zval,stats6.zval];
        
        % Friedman's Test
        if p <5
            % 14 steps per trial for wt
            reps = 16;
        elseif p>=5 && p<9
            % 6 per trial for rca
            reps = 6;
        else
            % 3 steps per trial for rcd
            reps = 3;
        end
%         ft_prep = rmmissing([val_hb.HB,val_lb.LB,val_hs.HS,val_ls.LS]);
%         [pd,tbl,stats] = friedman(ft_prep,reps); 
%         ft_results =[ft_results; z,p,tbl{2,3},tbl{2,5},pd];
%         c = multcompare(stats);
%         ft_results_comp =[ft_results_comp; repmat(z,size(c,1),1),repmat(p,size(c,1),1),c];
        
    end
end
%%        
%nORMALITY table results
%  Variable = vars(NT(:,1))';
%  Trial = trial(NT(:,2))';
%  normp_hb = NT(:,3);
%  normp_lb = NT(:,4);
%  normp_hs = NT(:,5);
%  normp_ls = NT(:,6);
%  normT = table(Variable,Trial, normp_hb, normp_lb, normp_hs,normp_ls);

%rmAnova table results
 Variable = vars(rmanova_results_matrix(:,1))';
 Trial = trial(rmanova_results_matrix(:,2))';
 DoF_time = rmanova_results_matrix(:,3);
 DoF_error = rmanova_results_matrix(:,4);
 fstat = rmanova_results_matrix(:,5);
 pval = rmanova_results_matrix(:,6);
 rmanova_results = table(Variable,Trial,DoF_time,DoF_error,fstat,pval);
 
 %Posthoc Mann-Whitney U table
 Variable = repmat(vars(phtt_matrix(:,1))',6,1);
 Trial = repmat(trial(phtt_matrix(:,2))',6,1);
 bag = {'HB','LB','HS','LS'};
 Bags = [[repmat(bag{1},length(phtt_matrix)*3,1);repmat(bag{2},length(phtt_matrix)*2,1);repmat(bag{4},length(phtt_matrix),1)],[repmat(bag{2},length(phtt_matrix),1);repmat(bag{3},length(phtt_matrix),1);repmat(bag{4},length(phtt_matrix),1);repmat(bag{3},length(phtt_matrix),1);repmat(bag{4},length(phtt_matrix),1);repmat(bag{3},length(phtt_matrix),1)]];
 pval = [phtt_matrix(:,3);phtt_matrix(:,4);phtt_matrix(:,5);phtt_matrix(:,6);phtt_matrix(:,7);phtt_matrix(:,8)];
 zval = [phtt_matrix(:,9);phtt_matrix(:,10);phtt_matrix(:,11); phtt_matrix(:,12);phtt_matrix(:,13); phtt_matrix(:,14)];
 N = zeros(length(repmat(phtt_matrix(:,2),6,1)),1);
 N(repmat(phtt_matrix(:,2),6,1)==1)=16*22;
 N(repmat(phtt_matrix(:,2),6,1)==5)=6*22;
 N(repmat(phtt_matrix(:,2),6,1)==9)=3*22;
 r_abs = abs(zval./sqrt(N)); 
 phMWU_results = table(Variable, Trial,Bags, pval,zval,N,r_abs);
 


%% Figures for Paper
%relevant configurations only
for x = 1:size(phMWU_results,1)
    ind_b(x,1) = strcmp(phMWU_results.Bags(x,:),'HBLB');
    ind_s(x,1) = strcmp(phMWU_results.Bags(x,:),'LSHS');
    ind_h(x,1) = strcmp(phMWU_results.Bags(x,:),'HBHS');
    ind_l(x,1) = strcmp(phMWU_results.Bags(x,:),'LBLS');
end
bagind = logical(ind_b + ind_s + ind_h + ind_l);
ind_B = ind_b(bagind);
ind_S = ind_s(bagind);
ind_H = ind_h(bagind);
ind_L = ind_l(bagind);

G1 = phMWU_results(bagind,:);
%Vertical Parameters
for gind = 1:size(G1,1)
    xx1=cell2mat(strfind(G1.Variable(gind,:),'_v_'));
    xx2=cell2mat(strfind(G1.Variable(gind,:),'skew_'));
    xx3=cell2mat(strfind(G1.Variable(gind,:),'kurt_'));
    xx4=cell2mat(strfind(G1.Variable(gind,:),'path_length_torso'));
    xx5=cell2mat(strfind(G1.Variable(gind,:),'range_jerk_torso'));
    xx6=cell2mat(strfind(G1.Variable(gind,:),'cadence'));
    if isempty(xx1) && isempty(xx2) && isempty(xx3) && isempty(xx4)&& isempty(xx5)&& isempty(xx6)
        ind_V(gind,1) = 0;
    else
        ind_V(gind,1) = 1;
    end
end
G = G1(~ind_V,:);
ind_B = ind_B(~ind_V);
ind_S = ind_S(~ind_V);
ind_H = ind_H(~ind_V);
ind_L = ind_L(~ind_V);
%Torso Parameters
for gind = 1:size(G,1)
    xx=cell2mat(strfind(G.Variable(gind,:),'torso'));
    if isempty(xx)
        ind_T(gind,1) = 0;
    else
        ind_T(gind,1) = 1;
    end
end

%WT Parameters
for gind = 1:size(G,1)
    xx=cell2mat(strfind(G.Trial(gind,:),'WT'));
    if isempty(xx)
        ind_WT(gind,1) = 0;
    else
        ind_WT(gind,1) = 1;
    end
end

%A Parameters
for gind = 1:size(G,1)
    xx=cell2mat(strfind(G.Trial(gind,:),'A'));
    if isempty(xx)
        ind_A(gind,1) = 0;
    else
        ind_A(gind,1) = 1;
    end
end

%D Parameters
for gind = 1:size(G,1)
    xx=cell2mat(strfind(G.Trial(gind,:),'D'));
    if isempty(xx)
        ind_D(gind,1) = 0;
    else
        ind_D(gind,1) = 1;
    end
end



pind = G.pval<=0.05;


fprintf('Of the %d variables, %d were torso metrics and %d were foot parameters.\n',length(ind_T),sum(ind_T),length(ind_T)-sum(ind_T));
fprintf('Of the %d torso metrics, %d were significant.\n',sum(ind_T),sum(ind_T(pind)));
fprintf('Of the %d foot metrics, %d were significant.\n',length(ind_T)-sum(ind_T),sum(~ind_T(pind)));

%% Figure 1: What tasks and locations matter?
tor_tot = [sum(ind_T&ind_WT),sum(ind_T&ind_A),sum(ind_T&ind_D)];
tor_sig = [sum(ind_T&ind_WT&pind),sum(ind_T&ind_A&pind),sum(ind_T&ind_D&pind)];
foot_tot = [sum(~ind_T&ind_WT),sum(~ind_T&ind_A),sum(~ind_T&ind_D)];
foot_sig = [sum(~ind_T&ind_WT&pind),sum(~ind_T&ind_A&pind),sum(~ind_T&ind_D&pind)];
figure('Renderer', 'painters', 'Position', [100 100 1100 517])
bar([tor_tot;foot_tot]')
hold on
bar([tor_sig;foot_sig]')
legend({'Torso','Foot','Torso*','Foot*'}) 
xticklabels({'WT','A','D'})
ylim([0 43])
ylabel('Number of Significant Comparisons')
title('Significant Torso and Foot Comparisons')
set(gca,'fontweight','Bold','fontsize',18)
tor_sig = [mean(G.r_abs(ind_T&ind_WT&pind)),mean(G.r_abs(ind_T&ind_A&pind)),mean(G.r_abs(ind_T&ind_D&pind))];
foot_sig = [mean(G.r_abs(~ind_T&ind_WT&pind)),mean(G.r_abs(~ind_T&ind_A&pind)),mean(G.r_abs(~ind_T&ind_D&pind))];


%% Figure 2: What tasks and bags matter?
tor_tot = [sum(ind_T&ind_B),sum(ind_T&ind_S),sum(ind_T&ind_H),sum(ind_T&ind_L)];
tor_sig = [sum(ind_T&ind_B&pind),sum(ind_T&ind_S&pind),sum(ind_T&ind_H&pind),sum(ind_T&ind_L&pind)];
foot_tot = [sum(~ind_T&ind_B),sum(~ind_T&ind_S),sum(~ind_T&ind_H),sum(~ind_T&ind_L)];
foot_sig = [sum(~ind_T&ind_B&pind),sum(~ind_T&ind_S&pind),sum(~ind_T&ind_H&pind),sum(~ind_T&ind_L&pind)];
figure('Renderer', 'painters', 'Position', [100 100 1100 517])
bar([tor_tot;foot_tot]')
hold on
bar([tor_sig;foot_sig]')
legend({'Torso','Foot','Torso*','Foot*'}) 
xticklabels({'Backpacks','Shoulder Bags','Heavy Bags','Light Bags'})
ylabel('Number of Significant Comparisons')
title('Significant Torso and Foot Comparisons')
ylim([0 32])
set(gca,'fontweight','Bold','fontsize',18)
tor_sig = [mean(G.r_abs(ind_T&ind_B&pind)),mean(G.r_abs(ind_T&ind_S&pind)),mean(G.r_abs(ind_T&ind_H&pind)),mean(G.r_abs(ind_T&ind_L&pind))];
foot_sig = [mean(G.r_abs(~ind_T&ind_B&pind)),mean(G.r_abs(~ind_T&ind_S&pind)),mean(G.r_abs(~ind_T&ind_H&pind)),mean(G.r_abs(~ind_T&ind_L&pind))];

%% Create Table
sigT_WT = G(pind & ind_WT,[1,3,7,4]);
[~,ind_1] = sort(sigT_WT.r_abs);
sigT_WT = sigT_WT(flip(ind_1),:)

sigT_A = G(pind & ind_A,[1,3,7,4]);
[~,ind_2] = sort(sigT_A.r_abs);
sigT_A = sigT_A(flip(ind_2),:)

sigT_D = G(pind & ind_D,[1,3,7,4])
[~,ind_3] = sort(sigT_D.r_abs);
sigT_D = sigT_D(flip(ind_3),:)











