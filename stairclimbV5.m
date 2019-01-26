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
zed_thresh = 1.25;%1.25;
sigma_v = .001;%1e-3;% ZUPT measurement noise covariance matrix.

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
%% Processing Loop
for i=[1,2,4:23]%subject  
    startTime = str2double(data.(subject{i}).annotations.StartTimestamp_ms_(1))/1000;
    calendTime = str2double(data.(subject{i}).annotations.StopTimestamp_ms_(1))/1000;
    endTime = str2double(data.(subject{i}).annotations.StopTimestamp_ms_(end))/1000;
    subTrials = T(strcmp(T.Subject,subject{i}),:);
    for j=1:2 %Foot
        fprintf('%s::%s\n',subject{i},foot_side{j})
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
                
                % check with plot
                %figure;
                %plot(tS,gSfilt)
                %hold on
                %plot(HS,HS_pks,'*m')        
                %plot(TO,TO_pks,'*g')
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
            
                % check number of steps
                NumSteps_HB_RC = [NumSteps_HB_RC;sum(HS>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(1)&HS<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(4)),sum(HS>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(5)&HS<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_RC(8))];
                NumSteps_LB_RC = [NumSteps_LB_RC;sum(HS>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(1)&HS<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(4)),sum(HS>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(5)&HS<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_RC(8))];
                NumSteps_HS_RC = [NumSteps_HS_RC;sum(HS>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(1)&HS<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(4)),sum(HS>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(5)&HS<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_RC(8))];
                NumSteps_LS_RC = [NumSteps_LS_RC;sum(HS>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(1)&HS<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(4)),sum(HS>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(5)&HS<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_RC(8))];

                NumSteps_HB_WAT = [NumSteps_HB_WAT;sum(HS>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_WAT(1)&HS<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_WAT(2)),sum(HS>=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_WAT(3)&HS<=subTrials(find(sum(subTrials.Load=='HB',2)==2),:).SegmentTimes_WAT(4))];
                NumSteps_LB_WAT = [NumSteps_LB_WAT;sum(HS>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_WAT(1)&HS<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_WAT(2)),sum(HS>=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_WAT(3)&HS<=subTrials(find(sum(subTrials.Load=='LB',2)==2),:).SegmentTimes_WAT(4))];
                NumSteps_HS_WAT = [NumSteps_HS_WAT;sum(HS>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_WAT(1)&HS<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_WAT(2)),sum(HS>=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_WAT(3)&HS<=subTrials(find(sum(subTrials.Load=='HS',2)==2),:).SegmentTimes_WAT(4))];
                NumSteps_LS_WAT = [NumSteps_LS_WAT;sum(HS>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_WAT(1)&HS<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_WAT(2)),sum(HS>=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_WAT(3)&HS<=subTrials(find(sum(subTrials.Load=='LS',2)==2),:).SegmentTimes_WAT(4))];

                NumSteps_Subject = [NumSteps_Subject;i];
                NumSteps_foot = [NumSteps_foot;j];
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

%                 figure;
%                 p1 = subplot(1,2,1);
%                 boxplot(zed(swingPhase))
%                 title(sprintf('%s::%s',subject{i},foot_side{j}))
%                 xlabel('Swing')
%                 p2 = subplot(1,2,2);
%                 boxplot(zed(~swingPhase))
%                 xlabel('Not Swing')
%                 linkaxes([p1 p2],'y')
                
                
%                 if j==2
%                     p1 = subplot(2,2,1);
%                     hist(zed(swingPhase))
%                     title(sprintf('%s::%s',subject{i},foot_side{j}))
%                     ylabel('Swing')
%                     p2 = subplot(2,2,3);
%                     hist(zed(~swingPhase))
%                     ylabel('Not Swing')
%                     linkaxes([p1 p2 p3 p4],'xy')
%                 else
%                     figure;
%                     p3 = subplot(2,2,2);
%                     hist(zed(swingPhase))
%                     title(sprintf('%s::%s',subject{i},foot_side{j}))
%                     ylabel('Swing')
%                     p4 = subplot(2,2,4);
%                     hist(zed(~swingPhase))
%                     ylabel('Not Swing')
%                 end
%                 

        %%
        % Plot for visual
% % %          figure;
% % %          plot(tS,gSfilt)
% % %          plot(tF,zed)
% % %         hold on
% % %         plot(tF(zind),zed(zind),'*')
% % %         fprintf('ZUPTs\n')

%            figure;
%           plot(tS,gSfilt)
%           hold on
%           plot(HS,HS_pks,'*m')        
%           plot(TO,TO_pks,'*g')
%           plot(tF,gF(:,3))
%           plot(tF(zind),zed(zind),'*k')
%           plot(tF,zed)
%           hold on
%         %plot(tF(zind),zed(zind),'*')
%        %fprintf('ZUPTs\n')
%            title(sprintf('Subject %d :: Foot %d',i,j))
%            legend('Swing','HS','TO')
%        
         
     % % figure;
% % p1 = subplot(1,4,1);
% % boxplot(NumSteps.NumSteps_HB_RC)
% % ylabel('Steps RC')
% % title('HB')
% % xticklabels({'A','D'})
% % p2 = subplot(1,4,2);
% % boxplot(NumSteps.NumSteps_LB_RC)
% % title('LB')
% % xticklabels({'A','D'})
% % p3 = subplot(1,4,3);
% % boxplot(NumSteps.NumSteps_HS_RC)
% % title('HS')
% % xticklabels({'A','D'})
% % p4 = subplot(1,4,4);
% % boxplot(NumSteps.NumSteps_LS_RC)
% % title('LS')
% % xticklabels({'A','D'})
% % linkaxes([p1 p2 p3 p4],'y')
% % 
% % figure;
% % p1 = subplot(1,4,1);
% % boxplot(NumSteps.NumSteps_HB_WAT)
% % ylabel('Steps WAT')
% % title('HB')
% % xticklabels({'L1','L2'})
% % p2 = subplot(1,4,2);
% % boxplot(NumSteps.NumSteps_LB_WAT)
% % title('LB')
% % xticklabels({'L1','L2'})
% % p3 = subplot(1,4,3);
% % boxplot(NumSteps.NumSteps_HS_WAT)
% % title('HS')
% % xticklabels({'L1','L2'})
% % p4 = subplot(1,4,4);
% % boxplot(NumSteps.NumSteps_LS_WAT)
% % title('LS')
% % xticklabels({'L1','L2'})
% % linkaxes([p1 p2 p3 p4],'y')

        %% Transpose Foot Data for Kalman Filter
        %Right Foot
        acc_s=aF';
        gyro_s=gF'*pi/180; %(rad/s)
        timestamp=tF';
    
        % Set data size for kalman filter
        data_size = length(acc_s); % Data size from foot

        %% Orientation from accelerometers. Sensor is assumed to be stationary.
        pitch = mean(-asin(acc_s(1,1)/g));
        roll = mean(atan(acc_s(2,1)./acc_s(3,1)));
        yaw = 0;

        C=[cos(pitch) sin(roll)*sin(pitch) cos(roll)*sin(pitch);...
            0 cos(roll) -sin(roll);...
            -sin(pitch) sin(roll)*cos(pitch) cos(roll)*cos(pitch)];

        C_prev = C;

        % Preallocate storage for heading estimate. Different from direction of
        % travel, the heading indicates the direction that the sensor, and therefore
        % the pedestrian, is facing.
        heading = nan(1, data_size);
        heading(1) = yaw;

        % Preallocate storage for accelerations in navigation frame.
        acc_n = nan(3, data_size);
        acc_n(:,1) = C*acc_s(:,1);

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
                if  zind(t)==1

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
             
fprintf('KF Done\n')
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



   
%%


            %% Perform offline step detection and computation of gait parameters
            %Initialize variables
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
            bag_type=[];
            fside=[];
            start_time = [];
            end_time = [];
            Type_Activity = [];
            Subject = [];

                  figure;
            hold on; grid on; axis equal; view([0,-1,0]); xlabel('X'); ylabel('Y'); zlabel('Z');
            
            %divide into steps
            for s=2:length(HS)
                ind_step = tF>=HS(s-1) & tF<=HS(s);

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

                % Extract spatial variables
                start_time = [start_time; HS(s)];
                end_time = [end_time; HS(s)];
                step_length = [step_length; pos_r(end,1)];
                lateral_dev = [lateral_dev; range(pos_r(:,2))];
                step_height = [step_height; max(pos_r(:,3))];
                max_swing_vel = [max_swing_vel; max(sqrt(sum(vel_r.^2,2)))];
                foot_attack_angle = [foot_attack_angle; 180*atan2(vel_r(1,3),vel_r(1,1))/pi];

                % Extract temporal variables
                contact_time = [contact_time; TO(s-1)-HS(s-1)];
                step_time = [step_time; HS(s)-HS(s-1)];
                cadence = [cadence; 1/(HS(s)-HS(s-1))];
                
                % Extract Torso Acceleration Variables
%                 torso_accel_HS = 
%                 torso_accel_TO = 
%                 torso_max_accel = 
%                 torso_min_accel = 

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
                foot_attack_angle,contact_time,step_time,cadence);
            fprintf('Step Table Done\n')
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
