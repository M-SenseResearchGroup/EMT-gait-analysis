%% Strikefoot for Opal
%Strikefoot_main_script_v2.m by Dr. Ryan McGinnis adapted for opal sensors
%by Lara Weed
%% Load data from Opal 
load('subject_data');
bags={'hb','lb','hs','ls'};
for i=1 %:length(subj)% subjects
    for W=1:length(bags)%weights
        %% from 1 foot
        %Right foot
        ac = subj(i).rc.(cell2mat(bags(W))).dfr.a(:,2:4);
        gy = subj(i).rc.(cell2mat(bags(W))).dfr.g(:,2:4);
        q_s = MSdata.(cell2mat(subjects(i))).ObstacleWalk.(cell2mat(bags(W))).RF.orientation;
        timestamp=subj(i).rc.(cell2mat(bags(W))).dfr.a(:,1);

        %%
        % Define calibration period for example (includes standing and walking)
        figure;
        plot(timestamp,ac)
        gcal=ginput(2);
        ind_cal = timestamp>gcal(1) & timestamp<gcal(2);

        % determine sampling rate
        fs = 1/mean(diff(timestamp));

        % Extract filtered acceleration and angular velocity magnitude
        acc_s = ac; %m/s^2
        data_size = length(acc_s);

        % Transpose data
        acc_s = acc_s.'; %m/s^2
        gy = gy.'; %rad/s

        %% Apply basic 'anatomical' calibration to IMU data
        g = find_grav_dir(ac(ind_cal,:));
        R = get_anatomical_cal(acc_s,g);

        % Apply rotation to data and orientation estimate
        acc_s = (acc_s.'*R.').';
        gy = (gy.'*R.').';
        R_s = rotRot(quat2rotm(q_s),R.');

        %% Compute foot trajectory using zero velocity updates and Kalman filter

        % Estimate magnitude of acceleration due to gravity
        g = median(sqrt(sum((ac(ind_cal,:)).^2,2)));

        % Define acceleration in world frame
        acc_n = dcmRot(R_s,acc_s.').';

        % Preallocate storage for velocity and position
        vel_n = nan(size(acc_n));
        pos_n = nan(size(acc_n));

        % Set initial velocity and position - assumed to be zero and at origin
        vel_n(:,1) = [0; 0; 0];
        pos_n(:,1) = [0; 0; 0];

        % Preallocate storage for state (0=stance, 1=not stance)
        % And set initial state - assumed to be stance
        state = ones(1,data_size);
        state(1) = 0;

        % Gyroscope stance phase detection threshold.
        gyro_threshold = .35;

        % Define KF parameters
        % Initial error covariance matrix.
        P = zeros(6);

        % Process noise parameter - accelerometer noise
        sigma_a = 0.05; %accelerometer measurement noise std deviation 
                        %(std dev of accel signal during stationary calibration)

        % Define zero velocity update measurement matrix
        H = [zeros(3) eye(3)];

        % Define zero velocity update measurement noise covariance matrix
        sigma_v = 1e-3; %zero velocity measurement noise std deviation
        R = diag([sigma_v sigma_v sigma_v]).^2;

        % Define and initialize state smoothing parameters
        wmax = 0; 
        start_swing = 1;
        start_stance = 1;
        swing_thresh = 2; 
        min_stance_time = 0.05; %min seconds in stance 

        % Preallocate storage for smoothed state
        statef = state;

        %% Apply KF to data
        for t = 2:data_size
            % Extract time step
            dt = timestamp(t) - timestamp(t-1);

            % Estimate velocity and position (error prone) 
            vel_n(:,t) = vel_n(:,t-1) + dt*(acc_n(:,t) - [0; 0; g] );
            pos_n(:,t) = pos_n(:,t-1) + dt*vel_n(:,t);

            % Define state transition matrix
            F = [eye(3)         dt*eye(3); %pos
                 zeros(3,3)      eye(3) ]; %vel

            % Define process noise covariance matrix
            Q = diag([0 0 0 sigma_a sigma_a sigma_a]*dt).^2;

            % Propagate the error covariance matrix.
            P = F*P*F' + Q;

            % Stance phase detection and zero-velocity updates.
            if norm(gy(:,t)) < gyro_threshold %Stance detection
                %Define state as stance
                state(t) = 0;

                % Update error estimate
                K = (P*H')/(H*P*H' + R);       % Compute kalman gain
                P = (eye(6) - K*H)*P;          % Update error covariance matrix
                delta_x = K*vel_n(:,t);        % Update filter state

                % Correct position and velocity based on error estimates.
                vel_n(:,t)=vel_n(:,t)-delta_x(4:6);
                pos_n(:,t)=pos_n(:,t)-delta_x(1:3);
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
            if norm(gy(:,t))>wmax
                wmax = norm(gy(:,t));
            end
            if state(t)-state(t-1)==-1 %1(swing) -> 0(stance)
                % Deal with false positive swing
                if wmax < swing_thresh
                    statef(start_swing:t)=0;
                end

                %Start stance tracking
                start_stance = t;
            end
        end

        %% Perform offline step detection and computation of gait parameters

        % Initialize variables
        state_buf = [];
        tran_buf = [];
        swing_thresh = 2;
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

        figure;
        hold on; grid on; axis equal; view([0,-1,0]); xlabel('X'); ylabel('Y'); zlabel('Z');

        %% Detect steps and compute gait parameters
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
                    w_temp = sqrt(sum(gy(:,tran_buf(1):tran_buf(4)).'.^2,2));
                    a_temp = sqrt(sum(acc_s(:,tran_buf(1):tran_buf(4)).'.^2,2));
                    j_temp = diff(a_temp)./diff(timestamp(tran_buf(1):tran_buf(4)).');
                    trans_temp = tran_buf-tran_buf(1)+1; 
                    fs = 1/mean(diff(t_temp));

                    if max(w_temp(trans_temp(1):trans_temp(2))) > swing_thresh && ...
                        max(w_temp(trans_temp(3):trans_temp(4))) > swing_thresh %Detect step

                        % Identify temporal characteristics of step
                        ind_swing1 = round(mean(trans_temp(1:2)));
                        ind_swing2 = round(mean(trans_temp(3:4)));
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

                        % Plot
                        plot3(pos_r(:,1), pos_r(:,2), pos_r(:,3),'linewidth',2)

                        % Compute orientation of long axis of foot for foot strike angle
                        foot_ap_n = R_s(:,:,fs_ind) * [1; 0; 0];

                        % Extract spatial variables
                        step_length = [step_length; pos_r(end,1)];
                        lateral_dev = [lateral_dev; range(pos_r(:,2))];
                        step_height = [step_height; max(pos_r(:,3))];
                        max_swing_vel = [max_swing_vel; max(sqrt(sum(vel_r.^2,2)))];
                        foot_angle_strike = [foot_angle_strike; atand(foot_ap_n(3)/norm(foot_ap_n(1:2)))];
                        foot_attack_angle = [foot_attack_angle; atand(vel_r(fs_ind_rel,3)/vel_r(fs_ind_rel,1))];

                        % Extract temporal variables
                        contact_time = [contact_time; timestamp(gait_inds(end,2))-timestamp(gait_inds(end,1))];
                        step_time = [step_time; timestamp(gait_inds(end,3))-timestamp(gait_inds(end,1))];
                        cadence = [cadence; 1/(timestamp(gait_inds(end,3))-timestamp(gait_inds(end,1)))];

                    end
                end

                %Pop off oldest state from buffers
                state_buf = state_buf(2:4);
                tran_buf = tran_buf(2:4);
            end
        end

        %% Create table for step variables
        step_table = table(step_length,lateral_dev,step_height,max_swing_vel,...
            foot_angle_strike,foot_attack_angle,contact_time,step_time,cadence);

        % Display table in command window
        disp(step_table);

        %% Plot Kinematics
        ind_state = state==0;

        figure;
        box on;
        hold on; 
        pos_r = pos_n(1:2,:);
        plot(pos_r(1,:),pos_r(2,:),'LineWidth',2,'Color','r');
        plot(pos_r(1,ind_state),pos_r(2,ind_state),'+k');
        start = plot(pos_r(1,1),pos_r(2,1),'Marker','^','LineWidth',2,'LineStyle','none');
        stop = plot(pos_r(1,end),pos_r(2,end),'Marker','o','LineWidth',2,'LineStyle','none');
        xlabel('x (m)');
        ylabel('y (m)');
        title('Estimated 2D path');
        legend([start;stop],'Start','End');
        axis equal;
        grid on;
        hold off;

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
    end
end