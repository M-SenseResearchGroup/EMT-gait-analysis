%Data Division 
%% set up structures and filters

load('..\Data\subj_data.mat'); % import data

%below_threshold = (signal < threshold);
% first plot the sacrum data to find the point of turnd
WT2 = fdesign.lowpass('Fp,Fst,Ap,Ast',8,12,10,60,62.5);
WT3= designfilt('lowpassfir','PassbandFrequency',5, ...
         'StopbandFrequency',9,'PassbandRipple',1, ...
         'StopbandAttenuation',60,'DesignMethod','kaiserwin','SampleRate',62.5);
HdWT2 = design(WT2,'butter');

Turn = struct();
bag_type = {'hb','lb','hs','ls'};
row=struct();
%% Finding the Threshold Values Graphically
for i=1:length(subj)
    if i~=3
        %filter design inputs
        A1 = filter(HdWT2,sqrt(sum(subj(i).wt.hb.scm.g(1:end,3:end).^2,2)));
        B1 = filter(HdWT2,sqrt(sum(subj(i).wt.lb.scm.g(1:end,3:end).^2,2)));
        C1 = filter(HdWT2,sqrt(sum(subj(i).wt.hs.scm.g(1:end,3:end).^2,2)));
        D1 = filter(HdWT2,sqrt(sum(subj(i).wt.ls.scm.g(1:end,3:end).^2,2)));
        X1 =horzcat(subj(i).wt.hb.scm.g(:,1),A1);
        X2 =horzcat(subj(i).wt.lb.scm.g(:,1),B1);
        X3 =horzcat(subj(i).wt.hs.scm.g(:,1),C1);
        X4 =horzcat(subj(i).wt.ls.scm.g(:,1),D1);
        A = A1(1:(end/2)); % using the first half to find the turn
        B = B1(1:(end/2));
        C = C1(1:(end/2));
        D = D1(1:(end/2+100)); %creating a margin for turn detection
       
        figure();
        %plotting the neck (mdc) data to find turn point detection
        P1 = filter(HdWT2,sqrt(sum(subj(i).wt.hb.scm.g(1:end,3:end).^2,2)));
        subplot(4,2,1)
        plot(subj(i).wt.hb.scm.g(1:end,1),P1) 
        hold on
        mu0 = max(A);
        row1 = find(X1(:,2)==mu0);
        row(i).hb=string(row1-50);
        mu1 = X1(row1-50,1);
        hline = refline([0 mu0]);
        hline.Color = 'r';
        hold on
        plot(mu1,mu0,'r*')
        title(strcat('Heavy Backpack MDC Filtered Walk Subject:',num2str(i)))
        
        
        P2 = filter(HdWT2,sqrt(sum(subj(i).wt.lb.scm.g(1:end,3:end).^2,2)));
        subplot(4,2,3)
        plot(subj(i).wt.lb.scm.g(1:end,1),P2)
        hold on
        mu2 = max(B);
        row2 = find(X2(:,2)==mu2);
        row(i).lb=string(row2-50);
        mu3 = X2(row2-50,1);
        hline = refline([0 mu2]);
        hline.Color = 'r';
        hold on
        plot(mu3,mu2,'r*')
        title(strcat('Light Backpack MDC Filtered Walk Subject:',num2str(i)))
        
        P3 = filter(HdWT2,sqrt(sum(subj(i).wt.hs.scm.g(1:end,3:end).^2,2)));
        subplot(4,2,2)
        plot(subj(i).wt.hs.scm.g(1:end,1),P3)
        hold on
        mu4 = max(C);
        row3 = find(X3(:,2)==mu4);
        row(i).hs=string(row3-50);
        mu5 = X3(row3-50,1);
        hline = refline([0 mu4]);
        hline.Color = 'r';
        hold on
        plot(mu5,mu4,'r*')
        title(strcat('Heavy Side MDC Filtered Walk Subject:',num2str(i)))
        
        P4 = filter(HdWT2,sqrt(sum(subj(i).wt.ls.scm.g(1:end,3:end).^2,2)));
        subplot(4,2,4)
        plot(subj(i).wt.ls.scm.g(1:end,1),P4)
        hold on
        mu6 = max(D);
        row4 = find(X4(:,2)==mu6);
        row(i).ls=string(row4-50);
        mu7 = X4(row4-50,1);
        hline = refline([0 mu6]);
        hline.Color = 'r';
        hold on
        plot(mu7,mu6,'r*')
        title(strcat('Light Side MDC Filtered Walk Subject:',num2str(i)))
        

    end
end

  