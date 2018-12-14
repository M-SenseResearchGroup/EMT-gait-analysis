% Load data 
load('subject_data')

%% Gait Code for Walk and Turn
dWT = fdesign.lowpass('Fp,Fst,Ap,Ast',4,14,10,50,62.5)
HdWT = design(dWT,'butter');
for i=1:length(subj)
    if i~=3
        FTLHB = filter(HdWT,subj(i).wt.hb.lsl.g(1:end,3));
        FTLLB = filter(HdWT,subj(i).wt.lb.lsl.g(1:end,3));
        FTLHS = filter(HdWT,subj(i).wt.hs.lsl.g(1:end,3));
        FTLLS = filter(HdWT,subj(i).wt.ls.lsl.g(1:end,3));
        FTRHB = filter(HdWT,subj(i).wt.hb.lsr.g(1:end,3));
        FTRLB = filter(HdWT,subj(i).wt.lb.lsr.g(1:end,3));
        FTRHS = filter(HdWT,subj(i).wt.hs.lsr.g(1:end,3));
        FTRLS = filter(HdWT,subj(i).wt.ls.lsr.g(1:end,3));
        
        %WT Shank Left Side
        figure();
        
        subplot(4,2,1)
        plot(subj(i).wt.hb.lsl.g(1:end,1),FTLHB)
        title(strcat('Heavy Backpack Left Shank Filtered Walk Subject:',num2str(i)))
        ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
        xlabel('t (s)');

        subplot(4,2,3)
        plot(subj(i).wt.lb.lsl.g(1:end,1),FTLLB)
        title(strcat('Light Backpack Left Shank Filtered Walk Subject:',num2str(i)))
        ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
        xlabel('t (s)');
        
        subplot(4,2,5)
        plot(subj(i).wt.hs.lsl.g(1:end,1),FTLHS)
        title(strcat('Heavy Sidebag Left Shank Filtered Walk Subject:',num2str(i)))
        ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
        xlabel('t (s)');
        
        subplot(4,2,7)
        plot(subj(i).wt.ls.lsl.g(1:end,1),FTLLS)
        title(strcat('Light Sidebag Left Shank Filtered Walk Subject:',num2str(i)))
        ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
        xlabel('t (s)');
        
        %WT Shank Right Side
        
        subplot(4,2,2)
        plot(subj(i).wt.hb.lsr.g(1:end,1),FTRHB)
        title(strcat('Heavy Backpack Right Shank Filtered Walk Subject:',num2str(i)))
        ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
        xlabel('t (s)');
        
        subplot(4,2,4)
        plot(subj(i).wt.lb.lsr.g(1:end,1),FTRLB)
        title(strcat('Light Backpack Right Shank Filtered Walk Subject:',num2str(i)))
        ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
        xlabel('t (s)');
        
        subplot(4,2,6)
        plot(subj(i).wt.hs.lsr.g(1:end,1),FTRHS)
        title(strcat('Heavy Sidebag Right Shank Filtered Walk Subject:',num2str(i)))
        ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
        xlabel('t (s)');
        
        subplot(4,2,8)
        plot(subj(i).wt.ls.lsr.g(1:end,1),FTRLS)
        title(strcat('Light Sidebag Right Shank Filtered Walk Subject:',num2str(i)))
        ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
        xlabel('t (s)');
        

    end
end
%% Original Data Plotting
for i=1:length(subj)
    if i~=3
        %WT Shank Left Side
        figure();
        
        subplot(4,2,1)
        plot(subj(i).wt.hb.lsl.g(1:end,1),subj(i).wt.hb.lsl.g(1:end,3))
        title(strcat('Heavy Backpack Left Shank Original Walk Subject:', num2str(i)))
        ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
        xlabel('t (s)');

        subplot(4,2,3)
        plot(subj(i).wt.lb.lsl.g(1:end,1),subj(i).wt.lb.lsl.g(1:end,3))
        title(strcat('Light Backpack Left Shank Original Walk Subject:', num2str(i)))
        ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
        xlabel('t (s)');
        
        subplot(4,2,5)
        plot(subj(i).wt.hs.lsl.g(1:end,1),subj(i).wt.hs.lsl.g(1:end,3))
        title(strcat('Heavy Sidebag Left Shank Original Walk Subject:', num2str(i)))
        ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
        xlabel('t (s)');
        
        subplot(4,2,7)
        plot(subj(i).wt.ls.lsl.g(1:end,1),subj(i).wt.ls.lsl.g(1:end,3))
        title(strcat('Light Sidebag Left Shank Original Walk Subject:', num2str(i)))
        ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
        xlabel('t (s)');
        
        %WT Shank Right Side
        
        subplot(4,2,2)
        plot(subj(i).wt.hb.lsr.g(1:end,1),subj(i).wt.hb.lsr.g(1:end,3))
        title(strcat('Heavy Backpack Right Shank Original Walk Subject:', num2str(i)))
        ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
        xlabel('t (s)');
        
        subplot(4,2,4)
        plot(subj(i).wt.lb.lsr.g(1:end,1),subj(i).wt.lb.lsr.g(1:end,3))
        title(strcat('Light Backpack Right Shank Original Walk Subject:', num2str(i)))
        ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
        xlabel('t (s)');
        
        subplot(4,2,6)
        plot(subj(i).wt.hs.lsr.g(1:end,1),subj(i).wt.hs.lsr.g(1:end,3))
        title(strcat('Heavy Sidebag Right Shank Original Walk Subject:', num2str(i)))
        ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
        xlabel('t (s)');
        
        subplot(4,2,8)
        plot(subj(i).wt.ls.lsr.g(1:end,1),subj(i).wt.ls.lsr.g(1:end,3))
        title(strcat('Light Sidebag Right Shank Original Walk Subject:', num2str(i)))
        ylabel('$\omega$ ($\frac{rad}{s}$)','Interpreter','latex');
        xlabel('t (s)');
        

    end
end
%%
C_z_morse = cwt(FT,"morse"); %take the Continuous Wavelet Transform of the z-axis
% gyroscope data comparing to the Morse Wavelet

figure(2);%Create figure for contour plot
set(gcf,'pos',[100 100 1200 600]); 
[C2,h2] = contourf(abs(C_z_morse),'edgecolor','none'); % contour plot of CWT
ylabel('Scale');
xlabel('Data Point');
title('Morse Wavelet Contour Plot');

figure(1);
hold on;
plot(subj(10).wt.hb.lsl.g(1:end,1),abs(C_z_morse(5,:)),'DisplayName','Morse Scale 5'); %Add the 5th scale of the morse transform
legend('show');


C_z_morlet = cwt(FT,'amor'); %Take the CWT comparing to the Morlet Wavelet

figure(1); %add 5th scale to existing plot
plot(subj(10).wt.hb.lsl.g(1:end,1),abs(C_z_morlet(5,:)),'DisplayName','Morlet Scale 5'); % plot 5th morlet scale
legend('show');

figure(3); %create contour plot figure
set(gcf,'pos',[100 100 1200 600]); %update size and position
[C3,h3] = contourf(abs(C_z_morlet),'edgecolor','none'); %plot contour data
xlabel('Data Point');
ylabel('Scale');
title('Morlet Wavelet Contour Plot');

figure(1); %update plot again
plot(subj(10).wt.hb.lsl.g(1:end,1),abs(C_z_morlet(12,:)),'DisplayName','Morlet Scale 12'); %plot 12th morlet scale
legend('show');



%Importing Data
%Set Up Metric Sets
Stridelength_L= [];
Stridelength_R= [];
Stridecycle= [];
Stridetime = [];
Velocity = [];
jerk_L = [];
jerk_R = [];
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
%Gait Analysis Code


%Add Comment