%% Stair Climb Task
%% Load Data
load('C:\Users\Laraw\OneDrive - University of Vermont\UVM\Research\McGinnis\Rescue_Climb\data\Preprocessed\EMTdata.mat')
load('C:\Users\Laraw\OneDrive - University of Vermont\UVM\Research\McGinnis\Rescue_Climb\data\Activity Segmentation\SegmentTimes.mat')
locations ={'anterior_thigh_right','proximal_lateral_shank_left','dorsal_foot_left','proximal_lateral_shank_right','dorsal_foot_right','sacrum','anterior_thigh_left','medial_chest'}; 
subject = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14','S15','S16','S17','S18','S19','S20','S21','S22','S23','S24','S25'};
foot={'dorsal_foot_right','dorsal_foot_left'};
shin ={'proximal_lateral_shank_right', 'proximal_lateral_shank_left'};  
foot_side = {'R','L'};
bag={'hb','lb','hs','ls'};
%% Set thresholds
stepT = [];
DCheck_HZ=[];
DCheck=[];
NumSteps_HB_RC = [];
NumSteps_LB_RC = [];
NumSteps_HS_RC = [];
NumSteps_LS_RC = []; 
NumSteps_HB_WAT = [];
NumSteps_LB_WAT = [];
NumSteps_HS_WAT = [];
NumSteps_LS_WAT = []; 
NumSteps_Subject = [];
NumSteps_foot = [];
error_1 = [];
error_2 = [];
%% Processing Loop
for i=1%[1,2,4:23]%subject  
    startTime = str2double(data.(subject{i}).annotations.StartTimestamp_ms_(1))/1000;
    calendTime = str2double(data.(subject{i}).annotations.StopTimestamp_ms_(1))/1000;
    wt_cal_start = str2double(data.(subject{i}).annotations.StartTimestamp_ms_(find(strcmp(data.S01.annotations.EventType,'Walk and Turn'),1)))/1000;
    wt_cal_end = str2double(data.(subject{i}).annotations.StopTimestamp_ms_(find(strcmp(data.S01.annotations.EventType,'Walk and Turn'),1)))/1000;    
    endTime = str2double(data.(subject{i}).annotations.StopTimestamp_ms_(end))/1000;
    subTrials = T(strcmp(T.Subject,subject{i}),:);
    
     %% Torso
        % Load data
         tT_all = data.(subject{i}).medial_chest.time/1000; % Timestamps of measurements (seconds)
         aT_all = data.(subject{i}).medial_chest.accel*9.8; % Accelerations in sensor frame (m/s^2).
         gT_all = data.(subject{i}).medial_chest.gyro; % Rates of turn in sensor frame.

        % Truncate Torso
         aT = aT_all(tT_all>=startTime & tT_all<=endTime,:);
         gT = gT_all(tT_all>=startTime & tT_all<=endTime,:);
         tT = tT_all(tT_all>=startTime & tT_all<=endTime,:); 
         dtT = [diff(tT);0];
         mean_cal = mean(aT(tT<=calendTime,:));
         
         aT_RMS = sqrt(sum((aT-mean_cal).^2,2));
    
        % Rotate Torso
            % Type 1- rotate grave vec during standing, pca during walking
             rotTor_aa = vrrotvec(mean_cal,[0 0 1]);
             rotTor_mat = vrrotvec2mat(rotTor_aa);
             aTp = (rotTor_mat*aT')'; %aTp(:,3) = vertical, aTp(:,1:2) = hz plane. Note: do not use hz axes independently          
             
%             figure
%             p1 = subplot(2,1,1);
%             plot(tT,aT_RMS)
%             title(sprintf('Subject %d Torso',i)) 
% %             p2 = subplot(3,1,2);
% %             %plot(tT,aTp);%_1)
% %             %hold on
% %             plot(tT,aT_1)
%             p3 = subplot(2,1,2);
%             plot(tT,aTp)
%             linkaxes([p1 p3],'x')
            
             
             
    
    for j=1:2 %Foot
        fprintf('%s::%s\n',subject{i},foot_side{j})
        zed_thresh = 1.25; %sigs(sigs(:,1)==i & sigs(:,2)==j,3); %1.25;%sigs(sigs(:,1)==i & sigs(:,2)==j,3); 
        
        %% set variables
            %foot
            if (i==6 && j==1) %sensors placed on wrong side for feet
                tF_all = data.(subject{i}).dorsal_foot_left.time/1000; % Timestamps of measurements (seconds)
                aF_all = data.(subject{i}).dorsal_foot_left.accel*9.8; % Accelerations in sensor frame (m/s^2).
                gF_all = data.(subject{i}).dorsal_foot_left.gyro; % Rates of turn in sensor frame.
            elseif (i==6 && j==2)
                tF_all = data.(subject{i}).dorsal_foot_right.time/1000; % Timestamps of measurements (seconds)
                aF_all = data.(subject{i}).dorsal_foot_right.accel*9.8; % Accelerations in sensor frame (m/s^2).
                gF_all = data.(subject{i}).dorsal_foot_right.gyro; % Rates of turn in sensor frame.
            else % sensors placed where theyre supposed to be
                tF_all = data.(subject{i}).(foot{j}).time/1000; % Timestamps of measurements (seconds)
                aF_all = data.(subject{i}).(foot{j}).accel*9.8; % Accelerations in sensor frame (m/s^2).
                gF_all = data.(subject{i}).(foot{j}).gyro; % Rates of turn in sensor frame.
            end
            
            %shin
            tS_all = data.(subject{i}).(shin{j}).time/1000; % Timestamps of measurements (seconds)
            aS_all = data.(subject{i}).(shin{j}).accel*9.8; % Accelerations in sensor frame (m/s^2).
            gS_all = data.(subject{i}).(shin{j}).gyro; % Rates of turn in sensor frame.
        
        %% Truncate
            %foot
            aF = aF_all(tF_all>=startTime & tF_all<=endTime,:);
            gF = gF_all(tF_all>=startTime & tF_all<=endTime,:);
            tF = tF_all(tF_all>=startTime & tF_all<=endTime,:);
            
            %shin
            aS = aS_all(tS_all>=startTime & tS_all<=endTime,:);
            gS = gS_all(tS_all>=startTime & tS_all<=endTime,:);
            tS = tS_all(tS_all>=startTime & tS_all<=endTime,:);
        
        
        %%  sampling rate for foot
            fs = 1/mean(diff(tF));
   
        %% Data Segmentation
        % uses mediolateral shank gyro to find HS and TO
            %% lowpass filter
                fc = 7; % Cut off frequency
                [b,a] = butter(2,fc/(fs/2)); % Butterworth filter of order 2
                if j==1
                    gSfilt = filtfilt(b,a,gS(:,3));
                else
                    gSfilt = filtfilt(b,a,-gS(:,3));
                end

                [pkss,locss] = findpeaks(gSfilt,tS,'MinPeakHeight',85); %swing Peaks
                [pksh,locsh] = findpeaks(-gSfilt,tS,'MinPeakProminence',3); %HS TO Peaks

        %% find HS and TO
                %set variables
                HS=[];
                TO=[];
                HS_pks=[];
                TO_pks=[];
                
                %TO before first swing
                TO = [TO; max(locsh(locsh<locss(1)))];
                TO_pks = [TO_pks;-pksh(find(locsh<locss(1),1,'last'))];
                
                %Between first and last swing
                for kk=2:length(locss)
                    candidates = locsh(locsh>locss(kk-1) & locsh<locss(kk));
                    candidates_pks = -pksh(locsh>locss(kk-1) & locsh<locss(kk));
                    if length(candidates)>1
                        HS = [HS;candidates(1)];
                        HS_pks = [HS_pks;candidates_pks(1)];
                        TO = [TO;candidates(end)];
                        TO_pks = [TO_pks;candidates_pks(end)];
                    end
                end
                
                %HS after last swing
                HS = [HS; min(locsh(locsh>locss(end)))];
                HS_pks = [HS_pks;-pksh(find(locsh>locss(end),1,'first'))];
                temp_ge = [HS,ones(length(HS),1);TO,zeros(length(TO),1)];
                gait_events = sortrows(temp_ge);%1=HS, 0=TO
                

        %% Swing Phase Logical
                swingPhase = zeros(length(tF),1);
                for hh = 1:length(TO)
                    if hh<length(TO)
                        swingPhase = swingPhase + (tF>= TO(hh) & tF<=HS(hh));
                    else 
                        if HS(end)>TO(end)
                            swingPhase = swingPhase + (tF>= TO(hh) & tF<=HS(hh));
                        else
                            swingPhase = swingPhase + (tF>= TO(hh));
                        end
                    end
                end
                swingPhase = logical(swingPhase);
                
                %Check with plot
                %figure
                %plot(tS,gSfilt)
                %hold on
                %plot(tF(swingPhase),ones(1,length(tF(swingPhase))),'*')
            
               
        %% ZUPT Detection
        %Zero-velocity update detection using foot accleration and gyro
            %% "Euclidian Distance" from still  
                % set gravity for foot
                grav = mean(aF(tF<=calendTime,:));
                
                % Gyroscope bias for foot      
                gyro_bias= mean(gF(tF<=calendTime,:))'*pi/180;
                
                %z-score accel and gyro metrics so they have equal scale on
                %all axes
                zscorea=zscore([aF(:,1)-grav(1);aF(:,2)-grav(2);aF(:,3)-grav(3)]);
                zscoreg=zscore([gF(:,1)-gyro_bias(1);gF(:,2)-gyro_bias(2);gF(:,3)-gyro_bias(3)]);
                
                %resize 
                len=length(aF(:,1));
                za=[zscorea(1:len),zscorea(1+len:len+len),zscorea(1+len+len:len+len+len)];
                zg=[zscoreg(1:len),zscoreg(1+len:len+len),zscoreg(1+len+len:len+len+len)];
                
                %RMS
                zed=sqrt((sum((zg(:,:)).^2,2))+(sum((za(:,:)).^2,2)));

                % baseline zupt ind
                zind=zed<zed_thresh & zed>-zed_thresh & ~swingPhase;%-min(zpks);
                
        %% Pedestrian Tracking Kalman Filter
            acc_s=aF';
            gyro_s=gF'*pi/180; %(rad/s)
            timestamp=tF';
            cal_still = tF<=calendTime;
            sigma_v = 0.001; %sigs(sigs(:,1)==i & sigs(:,2)==j,4);%.001;%sigs(sigs(:,1)==i & sigs(:,2)==j,4);%
    
            [ acc_n, vel_n, pos_n, heading, distance ] = KF_6dof( acc_s, gyro_s, timestamp, zind, sigma_v, cal_still );

%             figure;
%             plot3(pos_n(:,1),pos_n(:,2),pos_n(:,3))
        %% check that KF is working
            WATHB_sind = find(tF>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_WAT(1),1);
            WATHB_eind = max(find(tF<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_WAT(4)));
            WATHB_HZ = sqrt(sum((pos_n(1:2,WATHB_eind)-pos_n(1:2,WATHB_sind)).^2));
            WATHB = sqrt(sum((pos_n(:,WATHB_eind)-pos_n(:,WATHB_sind)).^2));

            WATLB_sind = find(tF>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_WAT(1),1);
            WATLB_eind = max(find(tF<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_WAT(4)));
            WATLB_HZ = sqrt(sum((pos_n(1:2,WATLB_eind)-pos_n(1:2,WATLB_sind)).^2));
            WATLB = sqrt(sum((pos_n(:,WATLB_eind)-pos_n(:,WATLB_sind)).^2));

            WATHS_sind = find(tF>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_WAT(1),1);
            WATHS_eind = max(find(tF<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_WAT(4)));
            WATHS_HZ = sqrt(sum((pos_n(1:2,WATHS_eind)-pos_n(1:2,WATHS_sind)).^2));
            WATHS = sqrt(sum((pos_n(:,WATHS_eind)-pos_n(:,WATHS_sind)).^2));

            WATLS_sind = find(tF>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_WAT(1),1);
            WATLS_eind = max(find(tF<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_WAT(4)));
            WATLS_HZ = sqrt(sum((pos_n(1:2,WATLS_eind)-pos_n(1:2,WATLS_sind)).^2));
            WATLS = sqrt(sum((pos_n(:,WATLS_eind)-pos_n(:,WATLS_sind)).^2));

            RCHB_sind = find(tF>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(1),1);
            RCHB_eind = max(find(tF<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(4)));
            RCHB_HZ = sqrt(sum((pos_n(1:2,RCHB_eind)-pos_n(1:2,RCHB_sind)).^2));
            RCHB = sqrt(sum((pos_n(:,RCHB_eind)-pos_n(:,RCHB_sind)).^2));

            RCLB_sind = find(tF>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(1),1);
            RCLB_eind = max(find(tF<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(4)));
            RCLB_HZ = sqrt(sum((pos_n(1:2,RCLB_eind)-pos_n(1:2,RCLB_sind)).^2));
            RCLB = sqrt(sum((pos_n(:,RCLB_eind)-pos_n(:,RCLB_sind)).^2));

            RCHS_sind = find(tF>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(1),1);
            RCHS_eind = max(find(tF<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(4)));
            RCHS_HZ = sqrt(sum((pos_n(1:2,RCHS_eind)-pos_n(1:2,RCHS_sind)).^2));
            RCHS = sqrt(sum((pos_n(:,RCHS_eind)-pos_n(:,RCHS_sind)).^2));


            RCLS_sind = find(tF>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(1),1);
            RCLS_eind = max(find(tF<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(4)));
            RCLS_HZ = sqrt(sum((pos_n(1:2,RCLS_eind)-pos_n(1:2,RCLS_sind)).^2));
            RCLS = sqrt(sum((pos_n(:,RCLS_eind)-pos_n(:,RCLS_sind)).^2));


            DCheck_HZ = [DCheck_HZ; i,j,WATHB_HZ, WATLB_HZ, WATHS_HZ, WATLS_HZ, RCHB_HZ, RCLB_HZ, RCHS_HZ ,RCLS_HZ];
            DCheck = [DCheck; i,j,WATHB, WATLB, WATHS, WATLS, RCHB, RCLB, RCHS ,RCLS];


        %% Perform offline step detection and computation of gait parameters
            %Initialize variables
                %foot
                step_length = [];
                lateral_dev = [];
                step_height = [];
                max_swing_vel = [];
                foot_angle_strike = [];
                foot_attack_angle = [];
                contact_time = [];
                step_time = [];
                cadence = [];
                foot_strike_time=[];
                fside=[];
                start_time = [];
                end_time = [];
                
                % Torso- RMS
                mean_accel_torso = [];
                max_accel_torso = [];
                min_accel_torso = [];
                std_accel_torso = [];
                range_accel_torso = [];
                skew_accel_torso = [];
                kurt_accel_torso = [];
                
                % Torso- vertical
                mean_accel_v_torso = [];
                max_accel_v_torso = [];
                min_accel_v_torso = [];
                std_accel_v_torso = [];
                range_accel_v_torso = [];
                skew_accel_v_torso = [];
                kurt_accel_v_torso = [];
                
                %Torso- Horizontal
                mean_jerk_torso = [];
                max_jerk_torso = [];
                std_jerk_torso = [];
                range_jerk_torso = [];
                skew_jerk_torso = [];
                kurt_jerk_torso = [];
                path_length_torso =[];
                
                %Torso- raw
                path_length_3d_torso =[];
                 
                % Other
                bag_type=[];
                Type_Activity = [];
                Subject = [];

                  figure;
            hold on; grid on; axis equal; view([0,-1,0]); xlabel('X'); ylabel('Y'); zlabel('Z');
            
            %divide into steps
            for s=2:length(HS)
                ind_step = tF>=HS(s-1) & tF<=HS(s);
                ind_step_T = tT>=HS(s-1) & tT<=HS(s);

                % Rotate so that step is primarily in the x-direction
                coeff = pca(pos_n(1:2,ind_step).');
                pos_r = [pos_n(1:2,ind_step).'*coeff, pos_n(3,ind_step).'];
                vel_r = [vel_n(1:2,ind_step).'*coeff, vel_n(3,ind_step).'];

                % Correct step so that x,y position begins at origin and min z
                % position is 0
                pos_r = pos_r - ones(length(pos_r(:,1)),1) * [pos_r(1,1:2), min(pos_r(:,3))];

                % Ensure that step proceeds in positive x-direction
                if pos_r(end,1) < 0
                    % Rotate by 180 deg
                    R = [cos(pi), sin(pi), 0; -sin(pi), cos(pi), 0; 0, 0, 1];
                    pos_r = pos_r*R;
                    vel_r = vel_r*R;
                end

                % Extract spatial foot variables
                start_time = [start_time; HS(s)];
                end_time = [end_time; HS(s)];
                step_length = [step_length; pos_r(end,1)];
                lateral_dev = [lateral_dev; range(pos_r(:,2))];
                step_height = [step_height; max(pos_r(:,3))];
                max_swing_vel = [max_swing_vel; max(sqrt(sum(vel_r.^2,2)))];
                foot_attack_angle = [foot_attack_angle; 180*atan2(vel_r(1,3),vel_r(1,1))/pi];

                % Extract temporal variables
                contact_time = [contact_time; TO(s)-HS(s-1)];
                step_time = [step_time; HS(s)-HS(s-1)];
                cadence = [cadence; 60/(HS(s)-HS(s-1))];
                
                
                % Extract Torso Acceleration Variables
                 feat = torsoSigFeatsRMS(aT_RMS(ind_step_T));
                 mean_accel_torso = [mean_accel_torso; feat(1)];
                 max_accel_torso = [max_accel_torso; feat(2)];
                 min_accel_torso = [min_accel_torso; feat(3)];
                 std_accel_torso = [std_accel_torso; feat(4)];
                 range_accel_torso = [range_accel_torso; feat(5)];
                 skew_accel_torso = [skew_accel_torso; feat(6)];
                 kurt_accel_torso = [kurt_accel_torso; feat(7)];
                 
                 feat = torsoSigFeatsRMS(aTp(ind_step_T,3));
                 mean_accel_v_torso = [mean_accel_v_torso; feat(1)];
                 max_accel_v_torso = [max_accel_v_torso; feat(2)];
                 min_accel_v_torso = [min_accel_v_torso; feat(3)];
                 std_accel_v_torso = [std_accel_v_torso; feat(4)];
                 range_accel_v_torso = [range_accel_v_torso; feat(5)];
                 skew_accel_v_torso = [skew_accel_v_torso; feat(6)];
                 kurt_accel_v_torso = [kurt_accel_v_torso; feat(7)];
                  
                 feat = torsoSigFeatsRMS(.5*cumtrapz(((aTp(ind_step_T,1)./dtT(ind_step_T)).^2)+((aTp(ind_step_T,2)./dtT(ind_step_T)).^2)));
                 mean_jerk_torso = [mean_jerk_torso; feat(1)];
                 max_jerk_torso = [max_jerk_torso; feat(2)];
                 std_jerk_torso = [std_jerk_torso; feat(4)];
                 range_jerk_torso = [range_jerk_torso; feat(5)];
                 skew_jerk_torso = [skew_jerk_torso; feat(6)];
                 kurt_jerk_torso = [kurt_jerk_torso; feat(7)];
                 path_length_torso =[path_length_torso;sqrt(sum(diff(aTp(ind_step_T,1)).^2+diff(aTp(ind_step_T,2)).^2))];
                 path_length_3d_torso =[path_length_3d_torso;sqrt(sum(diff(aTp(ind_step_T,1)).^2+diff(aTp(ind_step_T,2)).^2+diff(aTp(ind_step_T,3)).^2))];
              

                %Label section
                fside=[fside;foot_side{j}];
                Subject = [Subject;subject{i}];
                %WAT HB
                if start_time(end)>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_WAT(1) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_WAT(4)% Walk and turn
                    bag_type=[bag_type;bag{1}];

                    if start_time(end)>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_WAT(1) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_WAT(2)%Walk
                        Type_Activity = [Type_Activity ; 'W'];%W=walk
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_WAT(2) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_WAT(3)%turn
                        Type_Activity = [Type_Activity ; 'T'];%W=walk and turn
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_WAT(3) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_WAT(4)%Walk
                        Type_Activity = [Type_Activity ; 'W'];%W=walk
                    end

                    plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'g','linewidth',2)
                %WAT LB
                elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_WAT(1) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_WAT(4)% Walk and turn
                    bag_type=[bag_type;bag{2}];

                    if start_time(end)>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_WAT(1) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_WAT(2)%Walk
                        Type_Activity = [Type_Activity ; 'W'];%W=walk
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_WAT(2) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_WAT(3)%turn
                        Type_Activity = [Type_Activity ; 'T'];%W=walk and turn
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_WAT(3) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_WAT(4)%Walk
                        Type_Activity = [Type_Activity ; 'W'];%W=walk
                    end

                    plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'g','linewidth',2)    

                %WAT HS
                elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_WAT(1) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_WAT(4)% Walk and turn
                    bag_type=[bag_type;bag{3}];

                    if start_time(end)>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_WAT(1) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_WAT(2)%Walk
                        Type_Activity = [Type_Activity ; 'W'];%W=walk
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_WAT(2) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_WAT(3)%turn
                        Type_Activity = [Type_Activity ; 'T'];%W=walk and turn
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_WAT(3) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_WAT(4)%Walk
                        Type_Activity = [Type_Activity ; 'W'];%W=walk
                    end

                    plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'g','linewidth',2)    

               %WAT LS
                elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_WAT(1) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_WAT(4)% Walk and turn
                    bag_type=[bag_type;bag{4}];

                    if start_time(end)>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_WAT(1) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_WAT(2)%Walk
                        Type_Activity = [Type_Activity ; 'W'];%W=walk
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_WAT(2) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_WAT(3)%turn
                        Type_Activity = [Type_Activity ; 'T'];%W=walk and turn
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_WAT(3) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_WAT(4)%Walk
                        Type_Activity = [Type_Activity ; 'W'];%W=walk
                    end

                    plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'g','linewidth',2)         
                    
                %RC HB
                elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(1) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(8)% Rescue Climb
                    bag_type=[bag_type;bag{1}];

                    if start_time(end)>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(1) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(2)%Ascent
                        Type_Activity = [Type_Activity ; 'A'];%
                        plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'r','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(2) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(3)%Landing
                        Type_Activity = [Type_Activity ; 'L'];%
                        %plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'b','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(3) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(4)%Ascent
                        Type_Activity = [Type_Activity ; 'A'];%
                        plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'r','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(4) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(5)%Landing
                        Type_Activity = [Type_Activity ; 'L'];% 
                        %plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'b','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(5) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(6)%Descent
                        Type_Activity = [Type_Activity ; 'D'];%  
                        plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'m','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(6) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(7)%Descent
                        Type_Activity = [Type_Activity ; 'L'];% 
                        %plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'b','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(7) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(8)%Descent
                        Type_Activity = [Type_Activity ; 'D'];%  
                        plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'m','linewidth',2)
                    end

                %RC LB
                elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(1) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(8)% Rescue Climb
                    bag_type=[bag_type;bag{2}];

                    if start_time(end)>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(1) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(2)%Ascent
                        Type_Activity = [Type_Activity ; 'A'];%
                        plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'r','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(2) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(3)%Landing
                        Type_Activity = [Type_Activity ; 'L'];%
                        %plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'b','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(3) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(4)%Ascent
                        Type_Activity = [Type_Activity ; 'A'];%
                        plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'r','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(4) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(5)%Landing
                        Type_Activity = [Type_Activity ; 'L'];%  
                        %plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'b','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(5) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(6)%Descent
                        Type_Activity = [Type_Activity ; 'D'];% 
                        plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'m','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(6) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(7)%landing
                        Type_Activity = [Type_Activity ; 'L'];%   
                        %plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'b','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(7) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(8)%Descent
                        Type_Activity = [Type_Activity ; 'D'];%    
                        plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'m','linewidth',2)
                    end
                   
                %RC HS
                elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(1) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(8)% Rescue Climb
                    bag_type=[bag_type;bag{3}];

                    if start_time(end)>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(1) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(2)%Ascent
                        Type_Activity = [Type_Activity ; 'A'];%
                        plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'r','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(2) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(3)%Landing
                        Type_Activity = [Type_Activity ; 'L'];%
                        %plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'b','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(3) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(4)%Ascent
                        Type_Activity = [Type_Activity ; 'A'];%
                        plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'r','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(4) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(5)%Landing
                        Type_Activity = [Type_Activity ; 'L'];%   
                        %plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'b','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(5) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(6)%Descent
                        Type_Activity = [Type_Activity ; 'D'];%  
                        plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'m','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(6) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(7)%landing
                        Type_Activity = [Type_Activity ; 'L'];%   
                       % plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'b','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(7) && end_time(end)<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(8)%Descent
                        Type_Activity = [Type_Activity ; 'D'];%   
                        plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'m','linewidth',2)
                    end    

                %RC LS
                elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(1) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(8)% Rescue Climb
                    bag_type=[bag_type;bag{4}];

                    if start_time(end)>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(1) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(2)%Ascent
                        Type_Activity = [Type_Activity ; 'A'];%
                        plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'r','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(2) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(3)%Landing
                        Type_Activity = [Type_Activity ; 'L'];%
                        %plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'b','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(3) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(4)%Ascent
                        Type_Activity = [Type_Activity ; 'A'];%
                        plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'r','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(4) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(5)%Landing
                        Type_Activity = [Type_Activity ; 'L'];%    
                        %plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'b','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(5) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(6)%Descent
                        Type_Activity = [Type_Activity ; 'D'];%  
                        plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'m','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(6) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(7)%landing
                        Type_Activity = [Type_Activity ; 'L'];%   
                        %plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'b','linewidth',2)
                    elseif start_time(end)>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(7) && end_time(end)<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(8)%Descent
                        Type_Activity = [Type_Activity ; 'D'];%  
                        plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'m','linewidth',2)
                    end

                else
                    Type_Activity = [Type_Activity ; 'N'];
                    bag_type=[bag_type;'NA'];
                    %plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'g','linewidth',2)
                end
            end
                            
            %% Create table for step variables
            step_table = table(Subject,fside,Type_Activity,bag_type,start_time,end_time,step_length,lateral_dev,step_height,max_swing_vel,...
                foot_attack_angle,contact_time,step_time,cadence, path_length_torso, path_length_3d_torso,mean_accel_torso, max_accel_torso, min_accel_torso, std_accel_torso,...
                range_accel_torso, skew_accel_torso, kurt_accel_torso, mean_accel_v_torso, max_accel_v_torso, min_accel_v_torso, std_accel_v_torso,...
                range_accel_v_torso, skew_accel_v_torso, kurt_accel_v_torso, mean_jerk_torso, max_jerk_torso, std_jerk_torso,...
                range_jerk_torso, skew_jerk_torso, kurt_jerk_torso);
            fprintf('     Step Table Done\n')
            title(sprintf('%s::%s\n',subject{i},foot_side{j}))
            
            %% Plot Kinematics
            if j==1
                col = 'r';
                figure(i+500);
            else
                col = 'b';
                figure(i+500); 
                hold on;
            end
            box on;
            hold on; 
            plot3(pos_n(1,:),pos_n(2,:),pos_n(3,:),'LineWidth',2,'Color',col);
            plot3(pos_n(1,zind),pos_n(2,zind),pos_n(3,zind),'+k');
            start = plot3(pos_n(1,1),pos_n(2,1),pos_n(3,1),'Marker','^','LineWidth',2,'LineStyle','none');
            stop = plot3(pos_n(1,end),pos_n(2,end),pos_n(3,end),'Marker','o','LineWidth',2,'LineStyle','none');
            xlabel('x (m)');
            ylabel('y (m)');
            ylabel('z (m)');
            title(sprintf('Subject %d',i))
%             title('Estimated 3D path');
            legend([start;stop],'Start','End');
            axis equal;
            grid on;
            hold off;
    
%             figure; 
%             subplot(311)
%             plot(timestamp',acc_n.');
%             xlabel('time (s)'); ylabel('acceleration (m/s^2)');
% 
%             subplot(312)
%             plot(timestamp.',vel_n.');
%             xlabel('time (s)'); ylabel('velocity (m/s)');
% 
%             subplot(313)
%             plot(timestamp.',pos_n.');
%             xlabel('time (s)'); ylabel('position (m)');

%             figure; 
%             subplot(311)
%             plot(timestamp.',sqrt(sum(acc_n.'.^2,2)));
%             xlabel('time (s)'); ylabel('|acceleration| (m/s^2)');
% 
%             subplot(312)
%             plot(timestamp.',sqrt(sum(vel_n.'.^2,2)));
%             xlabel('time (s)'); ylabel('speed (m/s)');
% 
%             subplot(313)
%             plot(timestamp.',sqrt(sum(pos_n(1:2,:).'.^2,2)));
%             xlabel('time (s)'); ylabel('distance (m)');

            % Plot altitude estimates.
%             figure;
%             box on;
%             hold on;
%             plot(distance,pos_n(3,:),'Linewidth',2, 'Color','b');
%             xlabel('Distance Travelled (m)');
%             ylabel('z (m)');
%             title('Estimated altitude');
%             grid;
% 
%             % Display lines representing true altitudes of each floor.
%             floor_colour = [0 0.5 0]; % Colour for lines representing floors.
%             floor_heights = [0 3.6 7.2 10.8]; % Altitude of each floor measured from the ground floor.
%             floor_names = {'A' 'B' 'C' 'D'};
%             lim = xlim;
%             for floor_idx = 1:length(floor_heights)
%                 line(lim, [floor_heights(floor_idx) floor_heights(floor_idx)], 'LineWidth', 2, 'LineStyle', '--', 'Color', floor_colour);
%             end
%             ax1=gca; % Save handle to main axes.
%             axes('YAxisLocation','right','Color','none','YTickLabel', floor_names, 'YTick', floor_heights,'XTickLabel', {});
%             ylim(ylim(ax1));
%             ylabel('Floor');
%             hold off;

 stepT = [stepT;step_table];   

    end
end
NumSteps = table(NumSteps_Subject,NumSteps_foot,NumSteps_HB_RC,NumSteps_LB_RC,NumSteps_HS_RC,NumSteps_LS_RC,NumSteps_HB_WAT,NumSteps_LB_WAT,NumSteps_HS_WAT,NumSteps_LS_WAT);
%disp(DCheck)
figure
p1 = subplot(1,2,1);
boxplot(DCheck(:,3:end))
title('DCheck: All Dimensions')
p2 = subplot(1,2,2);
boxplot(DCheck_HZ(:,3:end))
title('DCheck: HZ Dimensions')
linkaxes([p1 p2],'y')



