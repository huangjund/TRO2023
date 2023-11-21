rotx = @(t) [1 0 0; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)] ;
roty = @(t) [cos(t) 0 sin(t) ; 0 1 0 ; -sin(t) 0  cos(t)] ;
rotz = @(t) [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1] ;

RB2_euler = quat2eul(RB2(:,4:end),'ZYX')*180/pi;
RB3_euler = quat2eul(RB3(:,4:end),'ZYX')*180/pi;
R_imu1 = [];
R_rb2 = [];
R_imu2 = [];
R_rb3 = [];
R_rb2_imu1 = [];
R_rb2_imu2 = [];
for i = 1:length(IMU1)
    R_imu1(:,:,i) = rotz(IMU1(i,1)/180*pi)*roty(IMU1(i,2)/180*pi)*rotx(IMU1(i,3)/180*pi);
    R_rb2(:,:,i) = rotz(RB2_euler(i,1))*roty(RB2_euler(i,2))*rotx(RB2_euler(i,3));
    R_imu2(:,:,i) = rotz(IMU2(i,1)/180*pi)*roty(IMU2(i,2)/180*pi)*rotx(IMU2(i,3)/180*pi);
    R_rb3(:,:,i) = rotz(RB3_euler(i,1))*roty(RB3_euler(i,2))*rotx(RB3_euler(i,3));
    
    %R_rb2_temp = R_rb2(:,:,i); R_rb3_temp = R_rb3(:,:,i);
    R_rb2_imu1(:,:,i) = R_imu1(:,:,i)*R_rb2(:,:,i)';
    R_rb3_imu2(:,:,i) = R_imu2(:,:,i)*R_rb3(:,:,i)';
    
    Euler_rb2_imu1(i,:) = rotm2eul(R_rb2_imu1(:,:,i),'ZYX');
    Euler_rb3_imu2(i,:) = rotm2eul(R_rb3_imu2(:,:,i),'ZYX');
end

figure(5);
subplot(2,3,1); plot(LowPassFilter(Euler_rb2_imu1(:,1)/pi*180));
subplot(2,3,2); plot(LowPassFilter(Euler_rb2_imu1(:,2)/pi*180));
subplot(2,3,3); plot(LowPassFilter(Euler_rb2_imu1(:,3)/pi*180));
subplot(2,3,4); plot(LowPassFilter(Euler_rb3_imu2(:,1)/pi*180));
subplot(2,3,5); plot(LowPassFilter(Euler_rb3_imu2(:,2)/pi*180));
subplot(2,3,6); plot(LowPassFilter(Euler_rb3_imu2(:,3)/pi*180));


% figure(5);
% subplot(2,3,1); plot((Euler_rb2_imu1(:,1)/pi*180));
% subplot(2,3,2); plot((Euler_rb2_imu1(:,2)/pi*180));
% subplot(2,3,3); plot((Euler_rb2_imu1(:,3)/pi*180));
% subplot(2,3,4); plot((Euler_rb3_imu2(:,1)/pi*180));
% subplot(2,3,5); plot((Euler_rb3_imu2(:,2)/pi*180));
% subplot(2,3,6); plot((Euler_rb3_imu2(:,3)/pi*180));

function filtered = LowPassFilter(data)
    % Filter parameters
    Fs = 5;  % Sampling frequency (e.g., 1000 Hz)
    Fc = 0.1; % Cutoff frequency for the low-pass filter (e.g., 100 Hz)

    % Normalize the frequency
    Wn = Fc/(Fs/2);

    % Design the filter
    [n, Wn] = buttord(Wn, Wn+0.1, 3, 40); % Filter order and parameters (these can be adjusted)
    [b, a] = butter(n, Wn, 'low');

    % Apply the filter
    filtered = filter(b, a, data);

%     % Plot original and filtered data (optional)
%     figure();
%     plot(data);
%     hold on;
%     plot(filtered_data);
%     legend('Original Data', 'Filtered Data');

end