clear all;
clc;
close all;

rotx = @(t) [1 0 0; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)] ;
rotz = @(t) [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1] ;
roty = @(t) [cos(t) 0 sin(t) ; 0 1 0 ; -sin(t) 0  cos(t)] ;
%% plot data1
% [IMU1, IMU2, IMU3, IMU4, RB1, RB2, RB3] = readMatrix('./data/20231214SUMdata1.csv');
% % Extract quaternion part (last four columns)
% quaternions1 = RB1(:, 4:7);
% quaternions2 = RB2(:, 4:7);
% quaternions3 = RB3(:, 4:7);
% euler_angles1 = quat2eul(quaternions1,'ZYX')*180/pi;
% euler_angles2 = quat2eul(quaternions2,'ZYX')*180/pi;
% euler_angles3 = quat2eul(quaternions3,'ZYX')*180/pi;
%  
% figure();
% 
% % Plot Euler angles in the second row of subplots
% for i = 1:3
%     subplot(7, 3, i);
%     plot(IMU1(:, i));
%     title(['Euler Angle ' num2str(i)]);
% end
% for i = 1:3
%     subplot(7, 3, i+3);
%     plot(IMU2(:, i));
%     title(['Euler Angle ' num2str(i)]);
% end
% for i = 1:3
%     subplot(7, 3, i+6);
%     plot(IMU3(:, i));
%     title(['Euler Angle ' num2str(i)]);
% end
% for i = 1:3
%     subplot(7, 3, i+9);
%     plot(IMU4(:, i));
%     title(['Euler Angle ' num2str(i)]);
% end
% for i = 1:3
%     subplot(7, 3, i+12);
%     plot(euler_angles1(:, i));
%     title(['Euler Angle ' num2str(i)]);
% end
% for i = 1:3
%     subplot(7, 3, i+15);
%     plot(euler_angles2(:, i));
%     title(['Euler Angle ' num2str(i)]);
% end
% for i = 1:3
%     subplot(7, 3, i+18);
%     plot(euler_angles3(:, i));
%     title(['Euler Angle ' num2str(i)]);
% end

%% plot data2
% [IMU1, IMU2, IMU3, IMU4, RB1, RB2, RB3] = readMatrix('./data/20231214SUMdata2.csv');
% % Extract quaternion part (last four columns)
% quaternions1 = RB1(:, 4:7);
% quaternions2 = RB2(:, 4:7);
% quaternions3 = RB3(:, 4:7);
% euler_angles1 = quat2eul(quaternions1,'ZYX')*180/pi;
% euler_angles2 = quat2eul(quaternions2,'ZYX')*180/pi;
% euler_angles3 = quat2eul(quaternions3,'ZYX')*180/pi;
%  
% figure();
% 
% % Plot Euler angles in the second row of subplots
% for i = 1:3
%     subplot(7, 3, i);
%     plot(IMU1(:, i));
%     title(['Euler Angle ' num2str(i)]);
% end
% for i = 1:3
%     subplot(7, 3, i+3);
%     plot(IMU2(:, i));
%     title(['Euler Angle ' num2str(i)]);
% end
% for i = 1:3
%     subplot(7, 3, i+6);
%     plot(IMU3(:, i));
%     title(['Euler Angle ' num2str(i)]);
% end
% for i = 1:3
%     subplot(7, 3, i+9);
%     plot(IMU4(:, i));
%     title(['Euler Angle ' num2str(i)]);
% end
% for i = 1:3
%     subplot(7, 3, i+12);
%     plot(euler_angles1(:, i));
%     title(['Euler Angle ' num2str(i)]);
% end
% for i = 1:3
%     subplot(7, 3, i+15);
%     plot(euler_angles2(:, i));
%     title(['Euler Angle ' num2str(i)]);
% end
% for i = 1:3
%     subplot(7, 3, i+18);
%     plot(euler_angles3(:, i));
%     title(['Euler Angle ' num2str(i)]);
% end

%% plot data3
% [IMU1, IMU2, IMU3, IMU4, RB1, RB2] = readMatrix2('./data/20231214SUMdata3.csv');
% % Extract quaternion part (last four columns)
% quaternions1 = RB1(:, 4:7);
% quaternions2 = RB2(:, 4:7);
% euler_angles1 = quat2eul(quaternions1,'ZYX')*180/pi;
% euler_angles2 = quat2eul(quaternions2,'ZYX')*180/pi;
%  
% figure();
% 
% % Plot Euler angles in the second row of subplots
% for i = 1:3
%     subplot(6, 3, i);
%     plot(IMU1(10:end, i));
%     title(['Euler Angle ' num2str(i)]);
% end
% for i = 1:3
%     subplot(6, 3, i+3);
%     plot(IMU2(10:end, i));
%     title(['Euler Angle ' num2str(i)]);
% end
% for i = 1:3
%     subplot(6, 3, i+6);
%     plot(IMU3(10:end, i));
%     title(['Euler Angle ' num2str(i)]);
% end
% for i = 1:3
%     subplot(6, 3, i+9);
%     plot(IMU4(10:end, i));
%     title(['Euler Angle ' num2str(i)]);
% end
% for i = 1:3
%     subplot(6, 3, i+12);
%     plot(euler_angles1(10:end, i));
%     title(['Euler Angle ' num2str(i)]);
% end
% for i = 1:3
%     subplot(6, 3, i+15);
%     plot(euler_angles2(10:end, i));
%     title(['Euler Angle ' num2str(i)]);
% end

%% plot data5
% [IMU1, IMU2, IMU3, IMU4, RB1, RB2] = readMatrix2('./data/20231214SUMdata5.csv');
% % Extract quaternion part (last four columns)
% quaternions1 = RB1(:, 4:7);
% quaternions2 = RB2(:, 4:7);
% euler_angles1 = quat2eul(quaternions1,'ZYX')*180/pi;
% euler_angles2 = quat2eul(quaternions2,'ZYX')*180/pi;
% 
% for i = 1:length(IMU1)
%     [RIMU1_omghat(i,:), RIMU1_theta(i)] = AxisAng3(so3ToVec(MatrixLog3(rotz(IMU1(i,1)/180*pi)*roty(IMU1(i,2)/180*pi)*rotx(IMU1(i,3)/180*pi))));
%     [RIMU2_omghat(i,:), RIMU2_theta(i)] = AxisAng3(so3ToVec(MatrixLog3(rotz(IMU2(i,1)/180*pi)*roty(IMU2(i,2)/180*pi)*rotx(IMU2(i,3)/180*pi))));
%     [Rrb1_omghat(i,:), Rrb1_theta(i)] = AxisAng3(so3ToVec(MatrixLog3(rotz(euler_angles1(i,1)/180*pi)*roty(euler_angles1(i,2)/180*pi)*rotx(euler_angles1(i,3)/180*pi))));
% end
% 
% figure();
% 
% % Plot Euler angles in the second row of subplots
% for i = 1:3
%     subplot(6, 3, i);
%     plot(IMU1(10:end, i));
%     title(['Euler Angle ' num2str(i)]);
%     subplot(6, 3, i+3);
%     plot(IMU2(10:end, i));
%     title(['Euler Angle ' num2str(i)]);
%     subplot(6, 3, i+6);
%     plot(IMU3(10:end, i));
%     title(['Euler Angle ' num2str(i)]);
%     subplot(6, 3, i+9);
%     plot(IMU4(10:end, i));
%     title(['Euler Angle ' num2str(i)]);
%     subplot(6, 3, i+12);
%     plot(euler_angles1(10:end, i));
%     title(['Euler Angle ' num2str(i)]);
%     subplot(6, 3, i+15);
%     plot(euler_angles2(10:end, i));
%     title(['Euler Angle ' num2str(i)]);
% end
% 
% figure();
% plot(Rrb1_theta/pi*180);

function [m1,m2,m3,m4,m5,m6,m7] = readMatrix(filename)
% Read the matrix from the CSV file
fullMatrix = readmatrix(filename);

% Check if the matrix has enough columns
totalColumnsNeeded = 3*4 + 7 + 7 + 7;
if size(fullMatrix, 2) < totalColumnsNeeded
    error('The matrix does not have enough columns.');
end

% Divide the matrix into five matrices
m1 = fullMatrix(:, 1:3);         % Columns 1 to 4
m2 = fullMatrix(:, 4:6);         % Columns 5 to 8
m3 = fullMatrix(:, 7:9);        % Columns 9 to 15
m4 = fullMatrix(:, 10:12);       % Columns 16 to 22
m5 = fullMatrix(:, 14:20);       % Columns 23 to 29
m6 = fullMatrix(:, 21:27);       % Columns 23 to 29
m7 = fullMatrix(:, 28:34);       % Columns 23 to 29

end


function [m1,m2,m3,m4,m5,m6,m7] = readMatrix2(filename)
% Read the matrix from the CSV file
fullMatrix = readmatrix(filename);

% Check if the matrix has enough columns
totalColumnsNeeded = 3*4 + 7 + 7;
if size(fullMatrix, 2) < totalColumnsNeeded
    error('The matrix does not have enough columns.');
end

% Divide the matrix into five matrices
m1 = fullMatrix(:, 1:3);         % Columns 1 to 4
m2 = fullMatrix(:, 4:6);         % Columns 5 to 8
m3 = fullMatrix(:, 7:9);        % Columns 9 to 15
m4 = fullMatrix(:, 10:12);       % Columns 16 to 22
m5 = fullMatrix(:, 14:20);       % Columns 23 to 29
m6 = fullMatrix(:, 21:27);       % Columns 23 to 29

end

function filtered = LowPassFilter(data)
    % Filter parameters
    Fs = 50;  % Sampling frequency (e.g., 1000 Hz)
    Fc = 1; % Cutoff frequency for the low-pass filter (e.g., 100 Hz)

    % Normalize the frequency
    Wn = Fc/(Fs/2);

    % Design the filter
    [n, Wn] = buttord(Wn, Wn+0.1, 3, 40); % Filter order and parameters (these can be adjusted)
    [b, a] = butter(n, Wn, 'low');

    % Apply the filter
    filtered = filter(b, a, data);
end

