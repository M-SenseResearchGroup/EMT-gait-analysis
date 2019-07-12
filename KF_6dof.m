function [ acc_n, vel_n, pos_n, heading, distance ] = KF_6dof( acc_s, gyro_s, timestamp, zupt, sigma_v, cal_still )
%KF_6dof Pedestrian tracking kalman filter for 6 dof imu 
%   Inputs:
%          acc_s - 3xn acceleration in the sensor frame that starts still [m/s^2] 
%           
%          gyro_s - 3xn gyro in the sensor frame [rad/s]
%
%          timestamp - timestamp [s]
%
%          zupt - zero-vleocity update flag
%
%          sigma_v -  
%
%          cal_still - logical of time when sensor is still to to set 
%                      acceleration and gyro covariances
%
%   Outputs:
%          acc_n - 3xn acceleration in the navigational frame estimate [m/s^2]
%
%          vel_n - 3xn velocity in the navigational frame estimate [m/s]
%
%          pos_n - 3xn position in the navigational frame estimate [m]
%
%          heading - 
%
%          distance - 
%
% Created by Lara Weed 24 Mar 2019
% Adapted from Fischer et. al.
%----------------------------------------------------------------------------
    
    % determine gravity
    g = mean(sqrt(sum(acc_s(:,cal_still).^2))); 
    
    % detemine gyro bias
    gyro_bias= mean(gyro_s(:,cal_still),2);
    
    % Set data size for kalman filter
    data_size = size(acc_s,2); % Data size from foot
    
    % Determine first rotation such that gravity is pointing down
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
    sigma_omega = std(sqrt(sum(gyro_s(:,cal_still).^2))); 
    sigma_a = std(sqrt(sum(acc_s(:,cal_still).^2)));

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
        if zupt(t)==1

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

end

