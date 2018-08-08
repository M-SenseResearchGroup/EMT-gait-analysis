load('..\Data\subj_data.mat');

%% Gait Code for Walk and Turn
dWT = fdesign.lowpass('Fp,Fst,Ap,Ast',2,14,10,60,62.5);
HdWT = design(dWT,'butter');

for i=1:length(subj)
    if i~=2
        if i~=3
    %filter design
           FTLHB = filter(HdWT,sqrt(sum(subj(i).wt.hb.dfl.g(1:end,2:end).^2,2)));
           FTLLB = filter(HdWT,sqrt(sum(subj(i).wt.lb.dfl.g(1:end,2:end).^2,2)));
           FTLHS = filter(HdWT,sqrt(sum(subj(i).wt.hs.dfl.g(1:end,2:end).^2,2)));
           FTLLS = filter(HdWT,sqrt(sum(subj(i).wt.ls.dfl.g(1:end,2:end).^2,2)));
           FTRHB = filter(HdWT,sqrt(sum(subj(i).wt.hb.dfr.g(1:end,2:end).^2,2)));
           FTRLB = filter(HdWT,sqrt(sum(subj(i).wt.lb.dfr.g(1:end,2:end).^2,2)));
           FTRHS = filter(HdWT,sqrt(sum(subj(i).wt.hs.dfr.g(1:end,2:end).^2,2)));
           FTRLS = filter(HdWT,sqrt(sum(subj(i).wt.ls.dfr.g(1:end,2:end).^2,2)));
        
            %WT Shank Left Side
            figure();
        
            subplot(4,2,1)
            [LHB_pk,LHB_loc] = findpeaks(FTLHB(1:(str2double(row(i).hb))),'MinPeakDistance',50,'MinPeakHeight',200); %peak detectiom
            LHB_mag(i) = mean(LHB_pk);%magnitude of acceleration
            LHB_cycles = diff(LHB_loc); %finding cycle times
            LHB_mc(i) = mean(LHB_cycles); %finding mean cycle time
            LHB_cad(i) = 1/LHB_mc(i); %finding cadence 
            plot(subj(i).wt.hb.dfl.g(1:end,1),FTLHB,subj(i).wt.hb.dfl.g(LHB_loc),LHB_pk,'or') %plotting and detailing
            title(strcat('Heavy Backpack Left Shank Filtered Walk Subject:',num2str(i))) %title
            ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex'); %y axis
            xlabel('t (s)'); %x axis
            
            
            subplot(4,2,3)
            [LLB_pk,LLB_loc] = findpeaks(FTLLB(1:str2double(row(i).lb)),'MinPeakDistance',50,'MinPeakHeight',200);
            LLB_mag(i) = mean(LLB_pk);
            LLB_cycles = diff(LLB_loc);
            LLB_mc(i) = mean(LLB_cycles);
            LLB_cad(i) = 1/LLB_mc(i);
            plot(subj(i).wt.lb.dfl.g(1:end,1),FTLLB,subj(i).wt.lb.dfl.g(LLB_loc),LLB_pk,'or')
            title(strcat('Light Backpack Left Shank Filtered Walk Subject:',num2str(i)))
            ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
            xlabel('t (s)');
            
            
            subplot(4,2,5)
            [LHS_pk,LHS_loc] = findpeaks(FTLHS(1:str2double(row(i).hs)),'MinPeakDistance',50,'MinPeakHeight',200);
            LHS_mag(i) = mean(LHS_pk);
            LHS_cycles = diff(LHS_loc);
            LHS_mc(i) = mean(LHS_cycles);
            LHS_cad(i) = 1/LHS_mc(i);
            plot(subj(i).wt.hs.dfl.g(1:end,1),FTLHS,subj(i).wt.hs.dfl.g(LHS_loc),LHS_pk,'or')
            title(strcat('Heavy Sidebag Left Shank Filtered Walk Subject:',num2str(i)))
            ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
            xlabel('t (s)');
            
            subplot(4,2,7)
            [LLS_pk,LLS_loc] = findpeaks(FTLLS(1:str2double(row(i).ls)),'MinPeakDistance',50,'MinPeakHeight',200);
            LLS_mag(i) = mean(LLS_pk);
            LLS_cycles = diff(LLS_loc);
            LLS_mc(i) = mean(LLS_cycles);
            LLS_cad(i) = 1/LLS_mc(i);
            plot(subj(i).wt.ls.dfl.g(1:end,1),FTLLS,subj(i).wt.ls.dfl.g(LLS_loc),LLS_pk,'or')
            title(strcat('Light Sidebag Left Shank Filtered Walk Subject:',num2str(i)))
            ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
            xlabel('t (s)');
            
            %WT Shank Right Side
            
            subplot(4,2,2)
            [RHB_pk,RHB_loc] = findpeaks(FTRHB(1:str2double(row(i).hb)),'MinPeakDistance',50,'MinPeakHeight',200);
            RHB_mag(i) = mean(RHB_pk);
            RHB_cycles = diff(RHB_loc);
            RHB_mc(i) = mean(RHB_cycles);
            RHB_cad(i) = 1/RHB_mc(i);
            plot(subj(i).wt.hb.dfr.g(1:end,1),FTRHB,subj(i).wt.hb.dfr.g(RHB_loc),RHB_pk,'or')
            title(strcat('Heavy Backpack Right Shank Filtered Walk Subject:',num2str(i)))
            ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
            xlabel('t (s)');
            
            subplot(4,2,4)
            [RLB_pk,RLB_loc] = findpeaks(FTRLB(1:str2double(row(i).lb)),'MinPeakDistance',50,'MinPeakHeight',200);
            RLB_mag(i) = mean(RLB_pk);
            RLB_cycles = diff(RLB_loc);
            RLB_mc(i) = mean(RLB_cycles);
            RLB_cad(i) = 1/RLB_mc(i);
            plot(subj(i).wt.lb.dfr.g(1:end,1),FTRLB,subj(i).wt.lb.dfr.g(RLB_loc),RLB_pk,'or')
            title(strcat('Light Backpack Right Shank Filtered Walk Subject:',num2str(i)))
            ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
            xlabel('t (s)');
            
            subplot(4,2,6)
            [RHS_pk,RHS_loc] = findpeaks(FTRHS(1:str2double(row(i).hs)),'MinPeakDistance',50,'MinPeakHeight',200);
            RHS_mag(i) = mean(RHS_pk);
            RHS_cycles = diff(RHS_loc);
            RHS_mc(i) = mean(RHS_cycles);
            RHS_cad(i) = 1/RHS_mc(i);
            plot(subj(i).wt.hs.dfr.g(1:end,1),FTRHS,subj(i).wt.hs.dfr.g(RHS_loc),RHS_pk,'or')
            title(strcat('Heavy Sidebag Right Shank Filtered Walk Subject:',num2str(i)))
            ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
            xlabel('t (s)');
            
            subplot(4,2,8)
            [RLS_pk,RLS_loc] = findpeaks(FTRLS(1:str2double(row(i).ls)),'MinPeakDistance',50,'MinPeakHeight',200);
            RLS_mag(i) = mean(RLS_pk);
            RLS_cycles = diff(RLS_loc);
            RLS_mc(i) = mean(RLS_cycles);
            RLS_cad(i) = 1/RLS_mc(i);
            plot(subj(i).wt.ls.dfr.g(1:end,1),FTRLS,subj(i).wt.ls.dfr.g(RLS_loc),RLS_pk,'or')
            title(strcat('Light Sidebag Right Shank Filtered Walk Subject:',num2str(i)))
            ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
            xlabel('t (s)');
            
            
            
        end
    end
end
%% table
figure();
%tables for values calculated above
%Subject_labels = strcat('subject'(1:length(subj)), num2str(1:length(subj)))
table_mc = table(LHB_mc', LLB_mc', LHS_mc', LLS_mc', RHB_mc', RLB_mc', RHS_mc', RLS_mc',...
    'VariableNames',{'LHB_mc' 'LLB_mc' 'LHS_mc' 'LLS_mc' 'RHB_mc' 'RLB_mc' 'RHS_mc' 'RLS_mc'})
Overallmeancycle_table = table(mean(LHB_mc'), mean(LLB_mc'), mean(LHS_mc'), mean(LLS_mc'), mean(RHB_mc'), mean(RLB_mc'), mean(RHS_mc'), mean(RLS_mc'),...
    'VariableNames',{'LHB_mc' 'LLB_mc' 'LHS_mc' 'LLS_mc' 'RHB_mc' 'RLB_mc' 'RHS_mc' 'RLS_mc'})
Peak_table = table(LHB_mag', LLB_mag', LHS_mag', LLS_mag', RHB_mag', RLB_mag', RHS_mag', RLS_mag',...
    'VariableNames',{'LHB_mag' 'LLB_mag' 'LHS_mag' 'LLS_mag' 'RHB_mag' 'RLB_mag' 'RHS_mag' 'RLS_mag'})
Overallmeanpeak_table = table(mean(LHB_mag'), mean(LLB_mag'), mean(LHS_mag'), mean(LLS_mag'), mean(RHB_mag'), mean(RLB_mag'), mean(RHS_mag'), mean(RLS_mag'),...
    'VariableNames',{'LHB_mag' 'LLB_mag' 'LHS_mag' 'LLS_mag' 'RHB_mag' 'RLB_mag' 'RHS_mag' 'RLS_mag'})
%% plot analysis

% gaitSignal_HBL
% gaitSignal_LBL
% gaitSignal_HSL
% gaitSignal_LSL
% gaitSignal_HBR
% gaitSignal_LBR
% gaitSignal_HSR
% gaitSignal_LSR

% Original Data Plotting
% for i=1:length(subj)
%     if i~=3
%         %WT Shank Left Side
%         figure();
%         
%         subplot(4,2,1)
%         plot(subj(i).wt.hb.dfl.g(1:end,1),sqrt(sum(subj(i).wt.hb.dfl.g(1:end,2:end).^2,2)))
%         title(strcat('Heavy Backpack Left Shank Original Walk Subject:', num2str(i)))
%         ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
%         xlabel('t (s)');
% 
%         subplot(4,2,3)
%         plot(subj(i).wt.lb.dfl.g(1:end,1),sqrt(sum(subj(i).wt.lb.dfl.g(1:end,2:end).^2,2)))
%         title(strcat('Light Backpack Left Shank Original Walk Subject:', num2str(i)))
%         ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
%         xlabel('t (s)');
%         
%         subplot(4,2,5)
%         plot(subj(i).wt.hs.dfl.g(1:end,1),sqrt(sum(subj(i).wt.hs.dfl.g(1:end,2:end).^2,2)))
%         title(strcat('Heavy Sidebag Left Shank Original Walk Subject:', num2str(i)))
%         ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
%         xlabel('t (s)');
%         
%         subplot(4,2,7)
%         plot(subj(i).wt.ls.dfl.g(1:end,1),sqrt(sum(subj(i).wt.ls.dfl.g(1:end,2:end).^2,2)))
%         title(strcat('Light Sidebag Left Shank Original Walk Subject:', num2str(i)))
%         ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
%         xlabel('t (s)');
%         
%         %WT Shank Right Side
%         
%         subplot(4,2,2)
%         plot(subj(i).wt.hb.dfr.g(1:end,1),sqrt(sum(subj(i).wt.hb.dfr.g(1:end,2:end).^2,2)))
%         title(strcat('Heavy Backpack Right Shank Original Walk Subject:', num2str(i)))
%         ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
%         xlabel('t (s)');
%         
%         subplot(4,2,4)
%         plot(subj(i).wt.lb.dfr.g(1:end,1),sqrt(sum(subj(i).wt.lb.dfr.g(1:end,2:end).^2,2)))
%         title(strcat('Light Backpack Right Shank Original Walk Subject:', num2str(i)))
%         ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
%         xlabel('t (s)');
%         
%         subplot(4,2,6)
%         plot(subj(i).wt.hs.dfr.g(1:end,1),sqrt(sum(subj(i).wt.hs.dfr.g(1:end,2:end).^2,2)))
%         title(strcat('Heavy Sidebag Right Shank Original Walk Subject:', num2str(i)))
%         ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
%         xlabel('t (s)');
%         
%         subplot(4,2,8)
%         plot(subj(i).wt.ls.dfr.g(1:end,1),sqrt(sum(subj(i).wt.ls.dfr.g(1:end,2:end).^2,2)))
%         title(strcat('Light Sidebag Right Shank Original Walk Subject:', num2str(i)))
%         ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
%         xlabel('t (s)');
%         
% 
%     end
% end
% 
% %% Importing Data
% %Set Up Metric Sets
% Stridelength_L= [];
% Stridelength_R= [];
% Stridecycle= [];
% Stridetime = [];
% Velocity = [];
% jerk_L = [];
% jerk_R = [];
%Perform Functions To Fill Sets
% for i=1:length(subj)
%     if i~=3
%         Stridelength_L = [Stridelength_L; ...
%                 SL(subj(i).wt.hb.dfl.g), SL(subj(i).wt.lb.dfl.g),...
%                 SL(subj(i).wt.hs.dfl.g), SL(subj(i).wt.ls.dfl.g)] ;
%            Stridelength_R = [Stridelength_R; ...
%                 SL(subj(i).wt.hb.dfr.g), SL(subj(i).wt.lb.dfr.g),...
%                 SL(subj(i).wt.hs.dfr.g), SL(subj(i).wt.ls.dfr.g)] ;
%              Stridecycle = [Stridecycle; ...
%                 SC(subj(i).wt.hb.dfr.g,subj(i).wt.hb.dfl.g), SC(subj(i).wt.lb.dfr.g,subj(i).wt.lb.dfl.g),...
%                 SC(subj(i).wt.hs.dfr.g,subj(i).wt.hs.dfl.g), SC(subj(i).wt.ls.dfr.g,subj(i).wt.ls.dfl.g)] ;
%             Stridetime = [Stridetime; ...
%                 ST(subj(i).wt.hb.dfr.g,subj(i).wt.hb.dfl.g), ST(subj(i).wt.lb.dfr.g,subj(i).wt.lb.dfl.g),...
%                 ST(subj(i).wt.hs.dfr.g,subj(i).wt.hs.dfl.g), ST(subj(i).wt.ls.dfr.g,subj(i).wt.ls.dfl.g)] ;
%              velocity = [velocity; ...
%                 V(subj(i).wt.hb.dfr.g,subj(i).wt.hb.dfl.g), V(subj(i).wt.lb.dfr.g,subj(i).wt.lb.dfl.g),...
%                 V(subj(i).wt.hs.dfr.g,subj(i).wt.hs.dfl.g), V(subj(i).wt.ls.dfr.g,subj(i).wt.ls.dfl.g)] ;
%             jerk_L = [jerk_L; ...
%                 JERK(subj(i).wt.hb.dfl.a), JERK(subj(i).wt.lb.dfl.a),...
%                 JERK(subj(i).wt.hs.dfl.a), JERK(subj(i).wt.ls.dfl.a)] ;
%            jerk_R = [jerk_R; ...
%                 JERK(subj(i).wt.hb.dfr.a), JERK(subj(i).wt.lb.dfr.a),...
%                 JERK(subj(i).wt.hs.dfr.a), JERK(subj(i).wt.ls.dfr.a)] ;
%     end
% end




%% functions

function J=JERK(x)
    for i= 2:length(x)
        delta_t=x(i,1)-x((i-1),1);
        diff_X(i-1)=(x(i,2)-x((i-1),2))/delta_t;
        diff_Z(i-1)=(x(i,4)-x((i-1),4))/delta_t;
        
    end
    J=0.5*sum(diff_X.*diff_X+diff_Z.*diff_Z)*delta_t;
end