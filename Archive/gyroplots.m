load('Lsubj_data_1min.mat')

for i=13%[1,2,4:15]%subject
for bagtype=1:4

        activitytype=2;
        bag={'hb','lb','hs','ls'};
        activity={'bsb','rc','wt'};

        %Sacrum Data
        tsa = subjL(i).(cell2mat(activity(activitytype))).(cell2mat(bag(bagtype))).scm.a(:,1)/1000;
        gsa = subjL(i).(cell2mat(activity(activitytype))).(cell2mat(bag(bagtype))).scm.g(:,2:4);

        %lowpass filter
        lpsa = fdesign.lowpass('Fp,Fst,Ap,Ast',3,10,40,90,62.5);
        lpfiltsa = design(lpsa,'butter');
        gmaglpsa = filter(lpfiltsa,gmagsa);

        %Moving Median Filter
        gmaglpmm1sa = movmedian(gmaglpsa,200);
        gmaglpmmsa = movmedian(gmaglpmm1sa,300);

        % Indexes for activity
        standind=gmaglpmmsa<6;
        
        %Find Times Standing
        standtime=tsa(standind);
        standstart=standtime(1);
        standend=[];
        for kp=2:length(standtime)
            if standtime(kp)-standtime(kp-1)>.5
                standstart=[standstart;standtime(kp)];
            end
        end
        for kp=1:length(standtime)-1
            if standtime(kp+1)-standtime(kp)>.5
                standend=[standend;standtime(kp)];
            end
        end
        standend=[standend;standtime(end)];
        standtimes=[standstart,standend,standend-standstart];

        %Determine longest still time for kalman filter start
        gcalind=find(max(standtimes(:,3))==standtimes(:,3));
        
        for foot=1:2
            if foot==1 %Right
                aF = subjL(i).(cell2mat(activity(activitytype))).(cell2mat(bag(bagtype))).dfr.a(:,2:4)*9.8; % Accelerations in sensor frame (m/s^2).
                gF = subjL(i).(cell2mat(activity(activitytype))).(cell2mat(bag(bagtype))).dfr.g(:,2:4); % Rates of turn in sensor frame.
                tF = subjL(i).(cell2mat(activity(activitytype))).(cell2mat(bag(bagtype))).dfr.a(:,1)/1000; % Timestamps of measurements.
            else %Left  
                tF = subjL(i).(cell2mat(activity(activitytype))).(cell2mat(bag(bagtype))).dfl.a(:,1)/1000; % Timestamps of measurements.
                aF = subjL(i).(cell2mat(activity(activitytype))).(cell2mat(bag(bagtype))).dfl.a(:,2:4)*9.8; % Accelerations in sensor frame (m/s^2).
                gF = subjL(i).(cell2mat(activity(activitytype))).(cell2mat(bag(bagtype))).dfl.g(:,2:4); % Rates of turn in sensor frame.
            end
        end
        %determine gravity of right foot during calibration
        g = mean(aF(tF>=standtimes(gcalind,1) & tF<=standtimes(gcalind,2),:));

        plot(sqrt(sum((aF(:,:)-g).^2,2)))
        hold on
end
end