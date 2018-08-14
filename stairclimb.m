%% Stair Climb Task
%Load Data
%load('C:\Users\Lukas Adamowicz\Documents\Study Data\EMT\Lsubj_data\Lsubj_data.mat')
load('C:\Users\Lukas Adamowicz\Dropbox\Masters\Project - EMT Study\Gait Analysis\Data\subj_data.mat', 'subj*')
subjL = subj;

i=1;%subject
zed_thresh=1.25;
min_stance_time = 0.03; %min seconds in stance
swing_thresh = 2;
%looks good: 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15
%look at again: 1(left bad, right good)
%looks bad: 2, 7(No Standing),
bagt=[];
direction=[];
steplength=[];
lateraldev=[];
stepheight=[];
maxswingvel=[];
footattackangle=[];
contacttime=[];
steptime=[];
cad=[];
footside=[];

for bagtype=1%:4
    
    activitytype=3;
    bag={'hb','lb','hs','ls'};
    activity={'bsb','rc','wt'};
    
    % Load Data
    % time
    tsa = subjL(i).(cell2mat(activity(activitytype))).(cell2mat(bag(bagtype))).scm.a(:,1)/1000;
    % gyro sacrum data
    gsa = subjL(i).(cell2mat(activity(activitytype))).(cell2mat(bag(bagtype))).scm.g(:,2:4);
    
    % Determine when calibration occurs using sacrum gyro data
    % Sacrum RMS Gyro
    gmagsa=sqrt(sum((gsa(:,:)).^2,2));
    %
%     figure;
%     plot(tsa(tsa>=(tsa(1)+100)),gsa(tsa>=(tsa(1)+100),2))
%     hold on
%     plot(tsa(tsa>=(tsa(1)+100)),gmagsa(tsa>=(tsa(1)+100)))
    %RCtimes=ginput(8);
    
    %lowpass filter
    lpsa = fdesign.lowpass('Fp,Fst,Ap,Ast',3,10,40,90,62.5);
    lpfiltsa = design(lpsa,'butter');
    gmaglpsa = filter(lpfiltsa,gmagsa);
    
    %Moving Median Filter
    gmaglpmm1sa = movmedian(gmaglpsa,200);
    gmaglpmmsa = movmedian(gmaglpmm1sa,300);
    
    %Alternate filter
%     wo = 2/31.25;
%     bw = wo/1;
%     [b, a] = iirnotch(wo, .7);

%     [b, a] = ellip(5, 2, 20, 1.5/31.25);
%     freqz(b, a);
%     alt = filtfilt(b, a, gmagsa);
    alt1 = movmean(gmaglpsa, 200);
    alt = movmean(alt1, 300);
    
    % Indexes for activity
    standind=gmaglpmmsa<6;
    
    
    figure;
    plot(tsa,gmagsa)
    hold on
%     plot(tsa,gmaglpsa, 'linewidth', 3)
    plot(tsa, gmaglpmmsa, 'linewidth', 3)
    plot(tsa, alt, 'linewidth', 3)
    legend('mag', 'Lara filt', 'alt filt')
    hold on
    plot(tsa(standind),gmagsa(standind),'b*')
    hold on
%     plot(tsa(walkind),gmagsa(walkind),'r*')
end

%%
for i=1
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
    
    % Load Foot data
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
        %determine gravity of foot during calibration
        g = median(sqrt(sum((aF(tF>=standtimes(gcalind,1) & tF<=standtimes(gcalind,2),:)).^2,2)));
        grav=median(aF(tF>=standtimes(gcalind,1) & tF<=standtimes(gcalind,2),:));
        
        % Gyroscope bias, to be determined for each sensor.
        gyro_bias= median(gF(tF>=standtimes(gcalind,1) & tF<=standtimes(gcalind,2),:))'*pi/180;
        
        % sampling rate
        fs = 1/mean(diff(tF));
        
        %% Cut Data to begining of calibraiton
        
        %Foot
        aF=aF(tF>=standtimes(gcalind,1),:);
        gF=gF(tF>=standtimes(gcalind,1),:);
        tF=tF(tF>=standtimes(gcalind,1),:);
        
        
        %% "Euclidian Distance" from still
        zscorea=zscore([aF(:,1)-grav(1);aF(:,2)-grav(2);aF(:,3)-grav(3)]);
        zscoreg=zscore([gF(:,1)-gyro_bias(1);gF(:,2)-gyro_bias(2);gF(:,3)-gyro_bias(3)]);
        
        len=length(aF(:,1));
        
        za=[zscorea(1:len),zscorea(1+len:len+len),zscorea(1+len+len:len+len+len)];
        zg=[zscoreg(1:len),zscoreg(1+len:len+len),zscoreg(1+len+len:len+len+len)];
        
        zed=sqrt((sum((zg(:,:)).^2,2))+(sum((za(:,:)).^2,2)));
        
        %lowpass filter
        lpzrf = fdesign.lowpass('Fp,Fst,Ap,Ast',2,10,10,50,62.5);%4,15,20,30,62.5);
        lpzfiltrf = design(lpzrf,'butter');
        zedlp = filter(lpzfiltrf,zed);
        
        
        figure;
        plot(tF,zed)
        hold on
        plot(tF,-zedlp)
        
        zind=zed<zed_thresh & zed>-zed_thresh;%-min(zpks);
        plot(tF(zind),zed(zind),'*')
        
        %% Transpose Foot Data for Kalman Filter
        %Right Foot
        acc_s=aF';
        gyro_s=gF'*pi/180; %(rad/s)
        timestamp=tF';
        
        % Set data size for kalman filter
        data_size = length(acc_s); % Data size from foot
        
        state = ones(1,data_size);
        state(1) = 0;
        
        % Preallocate storage for smoothed state
        statef = state;
        
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
        sigma_omega = 1e-2; sigma_a = 1e-2;
        
        % ZUPT measurement matrix.
        H = [zeros(3) zeros(3) eye(3)];
        
        % ZUPT measurement noise covariance matrix.
        sigma_v = 1e-3;
        R = diag([sigma_v sigma_v sigma_v]).^2;
        
        wmax = 0;
        start_swing = 1;
        start_stance = 1;
        
        
        
        %%
        %Vzexpected=-.1:.001:.1;
        %Cost=zeros(length(Vzexpected),1);
        %for k=1:length(Vzexpected)
        % Main Loop
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
            if  zind(t)==1 || (sum(tF(t)>=standtimes(:,1) & tF(t)<=standtimes(:,2))==1)
                %still_ind(t)==1
                %astill_ind(t)==1 && gstill_ind(t)==1
                %sum(timers(t)>=HSTO(:,1) & timers(t)<=HSTO(:,2))==1
                %norm(gyro_srf(:,t)) < gyro_threshold
                %%% Start Kalman filter zero-velocity update %%%
                
                state(t) = 0; %Stance
                
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
            
            % Smooth states
            statef(t) = state(t);
            if state(t)-state(t-1)==1 %0(stance) -> 1(swing)
                % Deal with false positive stance
                if (t-start_stance)/fs <= min_stance_time
                    statef(start_stance:t)=1;
                end
                % Start swing tracking
                start_swing = t;
                wmax = 0;
            end
            
            if norm(gyro_s(:,t))>wmax
                wmax = norm(gyro_s(:,t));
            end
            
            if state(t)-state(t-1)==-1 %1(swing) -> 0(stance)
                % Deal with false positive swing
                if wmax < swing_thresh
                    statef(start_swing:t)=0;
                end
                
                %Start stance tracking
                start_stance = t;
            end
            
            heading(t) = atan2(C(2,1), C(1,1)); % Estimate and save the yaw of the sensor (different from the direction of travel). Unused here but potentially useful for orienting a GUI correctly.
            C_prev = C; % Save orientation estimate, required at start of main loop.
            
            % Compute horizontal distance.
            distance(1,t) = distance(1,t-1) + sqrt((pos_n(1,t)-pos_n(1,t-1))^2 + (pos_n(2,t)-pos_n(2,t-1))^2);
        end
        %    Cost(k)=(pos_n(3,end)-pos_n(3,1)).^2;
        %end
        
        %% Perform offline step detection and computation of gait parameters
        
        % Initialize variables
        state_buf = [];
        tran_buf = [];
        gait_inds = [];
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
        sectionstart=[];
        sectionend=[];
        bag_type=[];
        fside=[];
        
        figure;
        hold on; grid on; axis equal; view([0,-1,0]); xlabel('X'); ylabel('Y'); zlabel('Z');
        
        
        
        % Detect steps and compute gait parameters
        for t = 2:data_size
            % Buffer state transitions
            if statef(t)~=statef(t-1)
                state_buf = [state_buf, statef(t)];
                tran_buf = [tran_buf,t];
            end
            
            % Identify step if there have been enough transitions
            if length(state_buf)==4
                %Identify potential step
                if isequal(state_buf,[1,0,1,0]) % Step = Footstrike-Footstrike
                    % Find foot strike and toe-off times within states
                    t_temp = timestamp(tran_buf(1):tran_buf(4)).';
                    w_temp = sqrt(sum(gyro_s(:,tran_buf(1):tran_buf(4)).'.^2,2));
                    a_temp = sqrt(sum(acc_s(:,tran_buf(1):tran_buf(4)).'.^2,2));
                    j_temp = diff(a_temp)./diff(timestamp(tran_buf(1):tran_buf(4)).');
                    trans_temp = tran_buf-tran_buf(1)+1;
                    fs = 1/mean(diff(t_temp));
                    
                    if max(w_temp(trans_temp(1):trans_temp(2))) > swing_thresh && ...
                            max(w_temp(trans_temp(3):trans_temp(4))) > swing_thresh %Detect step
                        
                        % Identify temporal characteristics of step
                        ind_swing1 = round(mean(trans_temp(1:2)));
                        ind_swing2 = floor(mean(trans_temp(3:4)));%round(mean(trans_temp(3:4)));
                        ind_fs = find(j_temp==max(j_temp(ind_swing1:trans_temp(2))),1,'first'); %first foot strike
                        ind_fo = find(j_temp==max(j_temp(trans_temp(3):ind_swing2)),1,'first'); %foot off
                        ind_fs2 = find(j_temp==max(j_temp(ind_swing2:end)),1,'first'); %second foot strike
                        gait_ind_temp = [ind_fs, ind_fo, ind_fs2];
                        gait_inds = [gait_inds; gait_ind_temp+tran_buf(1)];
                        
                        % Extract gait parameters
                        start_ind = gait_inds(end,1); %1st Foot Strike
                        stop_ind = gait_inds(end,3); %2nd Foot Strike
                        fs_ind = gait_inds(end,3);
                        fs_ind_rel = fs_ind - start_ind + 1; %location of fs
                        ind_step = start_ind:stop_ind;
                        
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
                        
                        if tF(start_ind)>=(subjL(i).(cell2mat(activity(activitytype))).(cell2mat(bag(bagtype))).dfr.a(1,1)/1000)+100
                            
                            % Plot
                            plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'linewidth',2)
                            
                            % Compute orientation of long axis of foot for foot strike angle
                            %foot_ap_n = R_s(:,:,fs_ind) * [1; 0; 0];
                            
                            % Extract spatial variables
                            step_length = [step_length; pos_r(end,1)];
                            lateral_dev = [lateral_dev; range(pos_r(:,2))];
                            step_height = [step_height; max(pos_r(:,3))];
                            max_swing_vel = [max_swing_vel; max(sqrt(sum(vel_r.^2,2)))];
                            %foot_angle_strike = [foot_angle_strike; atan2(foot_ap_n(3)/norm(foot_ap_n(1:2)))];
                            foot_attack_angle = [foot_attack_angle; 180*atan2(vel_r(fs_ind_rel,3),vel_r(fs_ind_rel,1))/pi];
                            
                            % Extract temporal variables
                            contact_time = [contact_time; timestamp(gait_inds(end,2))-timestamp(gait_inds(end,1))];
                            step_time = [step_time; timestamp(gait_inds(end,3))-timestamp(gait_inds(end,1))];
                            cadence = [cadence; 1/(timestamp(gait_inds(end,3))-timestamp(gait_inds(end,1)))];
                            
                            %Label section
                            bag_type=[bag_type; bagtype];
                            fside=[fside;foot];
                            if tF(start_ind)>=RCtimes(1) && tF(start_ind)<RCtimes(2)%starts in ascent 1 is ascent
                                sectionstart=[sectionstart; 1];%Ascent 1
                                
                            elseif tF(start_ind)>=RCtimes(2) && tF(start_ind)<RCtimes(3)
                                sectionstart=[sectionstart; 0];%Landing
                                
                            elseif tF(start_ind)>=RCtimes(3) && tF(start_ind)<RCtimes(4)
                                sectionstart=[sectionstart; 1];%Ascent 2
                                
                            elseif tF(start_ind)>=RCtimes(4) && tF(start_ind)<RCtimes(5)
                                sectionstart=[sectionstart; 0];%Top Landing
                                
                            elseif tF(start_ind)>=RCtimes(5) && tF(start_ind)<RCtimes(6)
                                sectionstart=[sectionstart; 2];%Descent 1
                                
                            elseif tF(start_ind)>=RCtimes(6) && tF(start_ind)<RCtimes(7)
                                sectionstart=[sectionstart; 0];%Landing
                                
                            elseif tF(start_ind)>=RCtimes(7)
                                sectionstart=[sectionstart; 2];%Descent 2
                            else
                                sectionstart=[sectionstart; 0];%Transition
                            end
                            
                            if tF(stop_ind)>=RCtimes(1) && tF(stop_ind)<RCtimes(2)%starts in ascent 1 is ascent
                                sectionend=[sectionend; 1];%Ascent 1
                                
                            elseif tF(stop_ind)>=RCtimes(2) && tF(stop_ind)<RCtimes(3)
                                sectionend=[sectionend; 0];%Landing
                                
                            elseif tF(stop_ind)>=RCtimes(3) && tF(stop_ind)<RCtimes(4)
                                sectionend=[sectionend; 1];%Ascent 2
                                
                            elseif tF(stop_ind)>=RCtimes(4) && tF(stop_ind)<RCtimes(5)
                                sectionend=[sectionend; 0];%Top Landing
                                
                            elseif tF(stop_ind)>=RCtimes(5) && tF(stop_ind)<RCtimes(6)
                                sectionend=[sectionend; 2];%Descent 1
                                
                            elseif tF(stop_ind)>=RCtimes(6) && tF(stop_ind)<RCtimes(7)
                                sectionend=[sectionend; 0];%Landing
                                
                            elseif tF(stop_ind)>=RCtimes(7)
                                sectionend=[sectionend; 2];%Descent 2
                            else
                                sectionend=[sectionend; 0];%Transition
                            end
                        end
                    end
                end
                
                %Pop off oldest state from buffers
                state_buf = state_buf(2:4);
                tran_buf = tran_buf(2:4);
            end
        end
        sectionstart(1)=0;
        sectionend(end)=0;
        %% Create table for step variables
        step_table = table(fside,sectionstart,bag_type,sectionend,step_length,lateral_dev,step_height,max_swing_vel,...
            foot_attack_angle,contact_time,step_time,cadence);
        
        stairs=(step_table.sectionstart>0 & step_table.sectionend>0);
        stair_table=step_table(stairs,:);
        ascent_table=stair_table(stair_table.sectionstart==1,:);
        descent_table=stair_table(stair_table.sectionstart==2,:);
        
        
        % Display table in command window
        %disp(step_table);
        
        %% Plot Kinematics
        ind_state = state==0;
        
        figure;
        box on;
        hold on;
        plot3(pos_n(1,:),pos_n(2,:),pos_n(3,:),'LineWidth',2,'Color','r');
        %plot3(pos_n(1,:),pos_n(2,:),pos_n(3,:),'LineWidth',2,'Color','r');
        plot3(pos_n(1,ind_state),pos_n(2,ind_state),pos_n(3,ind_state),'+k');
        start = plot3(pos_n(1,1),pos_n(2,1),pos_n(3,1),'Marker','^','LineWidth',2,'LineStyle','none');
        stop = plot3(pos_n(1,end),pos_n(2,end),pos_n(3,end),'Marker','o','LineWidth',2,'LineStyle','none');
        xlabel('x (m)');
        ylabel('y (m)');
        ylabel('z (m)');
        title('Estimated 3D path');
        legend([start;stop],'Start','End');
        axis equal;
        grid on;
        hold off;
        %%
        figure;
        subplot(311)
        plot(timestamp',acc_n.');
        xlabel('time (s)'); ylabel('acceleration (m/s^2)');
        
        subplot(312)
        plot(timestamp.',vel_n.');
        xlabel('time (s)'); ylabel('velocity (m/s)');
        
        subplot(313)
        plot(timestamp.',pos_n.');
        xlabel('time (s)'); ylabel('position (m)');
        
        figure;
        subplot(311)
        plot(timestamp.',sqrt(sum(acc_n.'.^2,2)));
        xlabel('time (s)'); ylabel('|acceleration| (m/s^2)');
        
        subplot(312)
        plot(timestamp.',sqrt(sum(vel_n.'.^2,2)));
        xlabel('time (s)'); ylabel('speed (m/s)');
        
        subplot(313)
        plot(timestamp.',sqrt(sum(pos_n(1:2,:).'.^2,2)));
        xlabel('time (s)'); ylabel('distance (m)');
        
        %% Plot altitude estimates.
        figure;
        box on;
        hold on;
        plot(distance,pos_n(3,:),'Linewidth',2, 'Color','b');
        xlabel('Distance Travelled (m)');
        ylabel('z (m)');
        title('Estimated altitude');
        grid;
        
        % Display lines representing true altitudes of each floor.
        floor_colour = [0 0.5 0]; % Colour for lines representing floors.
        floor_heights = [0 3.6 7.2 10.8]; % Altitude of each floor measured from the ground floor.
        floor_names = {'A' 'B' 'C' 'D'};
        lim = xlim;
        for floor_idx = 1:length(floor_heights)
            line(lim, [floor_heights(floor_idx) floor_heights(floor_idx)], 'LineWidth', 2, 'LineStyle', '--', 'Color', floor_colour);
        end
        ax1=gca; % Save handle to main axes.
        axes('YAxisLocation','right','Color','none','YTickLabel', floor_names, 'YTick', floor_heights,'XTickLabel', {});
        ylim(ylim(ax1));
        ylabel('Floor');
        hold off;
        
        bagt=[bagt;bag_type(stairs)];
        footside=[footside;fside(stairs)];
        direction=[direction;sectionend(stairs)];
        steplength=[steplength;step_length(stairs)];
        lateraldev=[lateraldev;lateral_dev(stairs)];
        stepheight=[stepheight;step_height(stairs)];
        maxswingvel=[maxswingvel;max_swing_vel(stairs)];
        footattackangle=[footattackangle;foot_attack_angle(stairs)];
        contacttime=[contacttime;contact_time(stairs)];
        steptime=[steptime;step_time(stairs)];
        cad=[cad;cadence(stairs)];
        sub=i*ones(length(cad),1);
    end
    
end
T = table(sub,footside,bagt,direction,steplength,lateraldev,stepheight,maxswingvel,footattackangle,contacttime,steptime,cad);
disp(T);