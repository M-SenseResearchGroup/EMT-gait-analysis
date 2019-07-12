%% Stair Climb Task
%% Load Data
load('C:\Users\Laraw\Documents\UVM\Research\McGinnis\Rescue_Climb\EMTdata.mat')
load('C:\Users\Laraw\Documents\UVM\Research\McGinnis\Rescue_Climb\SegmentTimes.mat')
locations ={'anterior_thigh_right','proximal_lateral_shank_left','dorsal_foot_left','proximal_lateral_shank_right','dorsal_foot_right','sacrum','anterior_thigh_left','medial_chest'}; 
subject = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14','S15','S16','S17','S18','S19','S20','S21','S22','S23','S24','S25'};
foot={'dorsal_foot_right','dorsal_foot_left'};
shin ={'proximal_lateral_shank_right', 'proximal_lateral_shank_left'};  
foot_side = {'R','L'};
bag={'hb','lb','hs','ls'};
%% Set thresholds
%zed_thresh = 1.25;
    HOT_sig = [];
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
    for i=[1,2,4:23]%subject  
        tic;
        startTime = str2double(data.(subject{i}).annotations.StartTimestamp_ms_(1))/1000;
        calendTime = str2double(data.(subject{i}).annotations.StopTimestamp_ms_(1))/1000;
        endTime = str2double(data.(subject{i}).annotations.StopTimestamp_ms_(end))/1000;
        subTrials = T(strcmp(T.Subject,subject{i}),:);
        for j=1:2 %Foot
            fprintf('%s::%s\n',subject{i},foot_side{j})
            for zed_thresh = 1.25
                for sigma_v = 0.001:0.004:0.1%.001;%1e-3;% ZUPT measurement noise covariance matrix.
                    fprintf('sigma_v: %s\n',sigma_v)

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

                    %Torso
                    tT_all = data.(subject{i}).medial_chest.time/1000; % Timestamps of measurements (seconds)
                    aT_all = data.(subject{i}).medial_chest.accel*9.8; % Accelerations in sensor frame (m/s^2).
                    gT_all = data.(subject{i}).medial_chest.gyro; % Rates of turn in sensor frame.
                %% Truncate
                    %foot
                    aF = aF_all(tF_all>=startTime & tF_all<=endTime,:);
                    gF = gF_all(tF_all>=startTime & tF_all<=endTime,:);
                    tF = tF_all(tF_all>=startTime & tF_all<=endTime,:);

                    %shin
                    aS = aS_all(tS_all>=startTime & tS_all<=endTime,:);
                    gS = gS_all(tS_all>=startTime & tS_all<=endTime,:);
                    tS = tS_all(tS_all>=startTime & tS_all<=endTime,:);

                    %Torso
                    aT = aT_all(tT_all>=startTime & tT_all<=endTime,:);
                    gT = gT_all(tT_all>=startTime & tT_all<=endTime,:);
                    tT = tT_all(tT_all>=startTime & tT_all<=endTime,:);

                %% set gravity for foot
                    g = mean(sqrt(sum(aF(tF<=calendTime,:).^2,2))); 
                    grav = mean(aF(tF<=calendTime,:));

                %% Gyroscope bias for foot      
                    gyro_bias= mean(gF(tF<=calendTime,:))'*pi/180;

                %%  sampling rate for foot
                    fs = 1/mean(diff(tF));

                %% Data Segmentation
                %uses mediolateral shank gyro to find HS and TO
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


                        % check with plot
                        %figure;
                        %plot(tS,gSfilt)
                        %hold on
                        %plot(HS,HS_pks,'*m')        
                        %plot(TO,TO_pks,'*g')
                        %% check number of steps
        % %                 NumSteps_HB_RC = [NumSteps_HB_RC;sum(HS>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(1)&HS<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(4)),sum(HS>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(5)&HS<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(8))];
        % %                 NumSteps_LB_RC = [NumSteps_LB_RC;sum(HS>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(1)&HS<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(4)),sum(HS>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(5)&HS<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(8))];
        % %                 NumSteps_HS_RC = [NumSteps_HS_RC;sum(HS>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(1)&HS<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(4)),sum(HS>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(5)&HS<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(8))];
        % %                 NumSteps_LS_RC = [NumSteps_LS_RC;sum(HS>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(1)&HS<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(4)),sum(HS>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(5)&HS<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(8))];
        % % 
        % %                 NumSteps_HB_WAT = [NumSteps_HB_WAT;sum(HS>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_WAT(1)&HS<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_WAT(2)),sum(HS>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_WAT(3)&HS<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_WAT(4))];
        % %                 NumSteps_LB_WAT = [NumSteps_LB_WAT;sum(HS>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_WAT(1)&HS<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_WAT(2)),sum(HS>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_WAT(3)&HS<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_WAT(4))];
        % %                 NumSteps_HS_WAT = [NumSteps_HS_WAT;sum(HS>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_WAT(1)&HS<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_WAT(2)),sum(HS>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_WAT(3)&HS<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_WAT(4))];
        % %                 NumSteps_LS_WAT = [NumSteps_LS_WAT;sum(HS>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_WAT(1)&HS<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_WAT(2)),sum(HS>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_WAT(3)&HS<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_WAT(4))];
        % % 
        % %                 NumSteps_Subject = [NumSteps_Subject;i];
        % %                 NumSteps_foot = [NumSteps_foot;j];
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

                %% Transpose Foot Data for Kalman Filter
                %Right Foot
                acc_s=aF';
                gyro_s=gF'*pi/180; %(rad/s)
                timestamp=tF';

                % Set data size for kalman filter
                data_size = length(acc_s); % Data size from foot

                % Orientation from accelerometers. Sensor is assumed to be stationary.
                r = vrrotvec(acc_s(:,1),[0 0 sqrt(sum(acc_s(:,1).^2))]);
                C_prev = vrrotvec2mat(r);

                % Preallocate storage for heading estimate. Different from direction of
                % travel, the heading indicates the direction that the sensor, and therefore
                % the pedestrian, is facing.
                heading = nan(1, data_size);
                heading(1) = 0;

                % Preallocate storage for accelerations in navigation frame.
                acc_n = nan(3, data_size);
                acc_n(:,1) = C_prev*acc_s(:,1);

                % Preallocate storage for velocity (in navigation frame).
                % Initial velocity assumed to be zero.
                vel_n = nan(3, data_size);
                vel_n(:,1) = [0 0 0]';

                % Preallocate storage for position (in navigation frame).
                % Initial position arbitrarily set to the origin.
                pos_n = nan(3, data_size);
                pos_n(:,1) = [0 0 0]';

                % Preallocate storage for distance travelled used for altitude plots.
                distance = nan(1,data_size-1);
                distance(1) = 0;

                % Error covariance matrix.
                P = zeros(9);

                % Process noise parameter, gyroscope and accelerometer noise.
                sigma_omega = std(sqrt(sum(gF(tF<=calendTime,:).^2,2))*pi/180); 
                sigma_a = std(sqrt(sum(aF(tF<=calendTime,:).^2,2)));
                %sigma_v = std(sqrt(sum(aF(tF<=calendTime,:).^2,2)));

                % ZUPT measurement matrix.
                H = [zeros(3) zeros(3) eye(3)];

                R = diag([sigma_v sigma_v sigma_v]).^2;

        %% Main Loop
                    for t = 2:data_size
                        %%% Start INS (transformation, double integration) %%%
                        dt = timestamp(t) - timestamp(t-1);%ms

                        % Remove bias from gyro measurements.
                        gyro_s1 = gyro_s(:,t) - gyro_bias;

                        % Skew-symmetric matrix for angular rates
                        ang_rate_matrix = [0   -gyro_s1(3)   gyro_s1(2);
                                            gyro_s1(3)  0   -gyro_s1(1);
                                            -gyro_s1(2)  gyro_s1(1)  0];

                        % orientation estimation
                        C = C_prev*(2*eye(3)+(ang_rate_matrix*dt))/(2*eye(3)-(ang_rate_matrix*dt));

                        % Transforming the acceleration from sensor frame to navigation frame.
                        acc_n(:,t) = 0.5*(C + C_prev)*acc_s(:,t);

                        % Velocity and position estimation using trapeze integration.
                        vel_n(:,t) = vel_n(:,t-1) + ((acc_n(:,t) - [0; 0; g] )+(acc_n(:,t-1) - [0; 0; g]))*dt/2;
                        pos_n(:,t) = pos_n(:,t-1) + (vel_n(:,t) + vel_n(:,t-1))*dt/2;

                        % Skew-symmetric cross-product operator matrix formed from the n-frame accelerations.
                        S = [0  -acc_n(3,t)  acc_n(2,t);
                            acc_n(3,t)  0  -acc_n(1,t);
                            -acc_n(2,t) acc_n(1,t) 0];

                        % State transition matrix.
                        F = [eye(3)  zeros(3,3)    zeros(3,3);
                            zeros(3,3)   eye(3)  dt*eye(3);
                            -dt*S  zeros(3,3)    eye(3) ];

                        % Compute the process noise covariance Q.
                        Q = diag([sigma_omega sigma_omega sigma_omega 0 0 0 sigma_a sigma_a sigma_a]*dt).^2;

                        % Propagate the error covariance matrix.
                        P = F*P*F' + Q;
                        %%% End INS %%%

                        % Stance phase detection and zero-velocity updates.
                        if zind(t)==1 %zupt(t)==1 %zind(t)==1  

                            K = (P*(H)')/((H)*P*(H)' + R); % Kalman gain.

                            % Update the filter state.
                            delta_x = K*(vel_n(:,t));%+[0; 0; Vzexpected(k)]);

                            % Update the error covariance matrix.
                            P = (eye(9) - K*H)*P; % Simplified covariance update found in most books.

                            % Extract errors from the KF state.
                            attitude_error = delta_x(1:3);
                            pos_error = delta_x(4:6);
                            vel_error = delta_x(7:9);
                            %%% End Kalman filter zero-velocity update %%%

                            %%% Apply corrections to INS estimates. %%%
                            % Skew-symmetric matrix for small angles to correct orientation.
                            ang_matrix = -[0   -attitude_error(3,1)   attitude_error(2,1);
                                attitude_error(3,1)  0   -attitude_error(1,1);
                                -attitude_error(2,1)  attitude_error(1,1)  0];

                            % Correct orientation.
                            C = (2*eye(3)+(ang_matrix))/(2*eye(3)-(ang_matrix))*C;

                            % Correct position and velocity based on Kalman error estimates.
                            vel_n(:,t)=vel_n(:,t)-vel_error;
                            pos_n(:,t)=pos_n(:,t)-pos_error;
                        end

                        heading(t) = atan2(C(2,1), C(1,1)); % Estimate and save the yaw of the sensor (different from the direction of travel). Unused here but potentially useful for orienting a GUI correctly.
                        C_prev = C; % Save orientation estimate, required at start of main loop.

                        % Compute horizontal distance.
                        distance(1,t) = distance(1,t-1) + sqrt((pos_n(1,t)-pos_n(1,t-1))^2 + (pos_n(2,t)-pos_n(2,t-1))^2);
                    end

        fprintf('     KF Done\n')
        % figure;
        % plot(pos_n')
        % title(sprintf('%d::%d',i,j))

        % figure;
        % plot3(pos_n(1,:),pos_n(2,:),pos_n(3,:));
        % title(sprintf('Subject %d :: Foot %d',i,j))
        % axis equal

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
        cost = sum([WATHB, WATLB, WATHS, WATLS, RCHB, RCLB, RCHS ,RCLS]);

        HOT_sig = [HOT_sig; i, j, zed_thresh, sigma_v, cost];
        toc;
    %     %% Plot Kinematics
    %             if j==1
    %                 col = 'r';
    %                 figure();
    %             else
    %                 col = 'b';
    %                 figure(); 
    %                 hold on;
    %             end
    %             box on;
    %             hold on; 
    %             plot3(pos_n(1,:),pos_n(2,:),pos_n(3,:),'LineWidth',2,'Color',col);
    %             plot3(pos_n(1,zind),pos_n(2,zind),pos_n(3,zind),'+k');
    %             start = plot3(pos_n(1,1),pos_n(2,1),pos_n(3,1),'Marker','^','LineWidth',2,'LineStyle','none');
    %             stop = plot3(pos_n(1,end),pos_n(2,end),pos_n(3,end),'Marker','o','LineWidth',2,'LineStyle','none');
    %             xlabel('x (m)');
    %             ylabel('y (m)');
    %             ylabel('z (m)');
    %             title(sprintf('Subject %d',i))
    % %             title('Estimated 3D path');
    %             legend([start;stop],'Start','End');
    %             axis equal;
    %             grid on;
    %             hold off;
    %     
                end
            end
        end
    end


    
    