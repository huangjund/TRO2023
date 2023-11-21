clear all;
clc;
close all;

[IMU1, IMU2, RB1, RB2, RB3] = readMatrix('./data/SUMdata2.csv');

%% normalization
RB2(:,1:3) = RB2(:,1:3)-RB1(:,1:3);
RB3(:,1:3) = RB3(:,1:3)-RB1(:,1:3);

%% data filter
RB2(:,4) = LowPassFilter(RB2(:,4));
RB2(:,5) = LowPassFilter(RB2(:,5));
RB2(:,6) = LowPassFilter(RB2(:,6));
RB2(:,7) = LowPassFilter(RB2(:,7));
RB3(:,4) = LowPassFilter(RB3(:,4));
RB3(:,5) = LowPassFilter(RB3(:,5));
RB3(:,6) = LowPassFilter(RB3(:,6));
RB3(:,7) = LowPassFilter(RB3(:,7));

%% plot
subplotIMU(IMU1,'IMU1');
subplotIMU(IMU2,'IMU2');
subplotRB(RB2,'RB2');
subplotRB(RB3,'RB3');



function [m1,m2,m3,m4,m5] = readMatrix(filename)
% Read the matrix from the CSV file
fullMatrix = readmatrix(filename);

% Check if the matrix has enough columns
totalColumnsNeeded = 4 + 4 + 7 + 7 + 7;
if size(fullMatrix, 2) < totalColumnsNeeded
    error('The matrix does not have enough columns.');
end

% Divide the matrix into five matrices
m1 = fullMatrix(:, 1:4);         % Columns 1 to 4
m2 = fullMatrix(:, 5:8);         % Columns 5 to 8
m3 = fullMatrix(:, 9:15);        % Columns 9 to 15
m4 = fullMatrix(:, 16:22);       % Columns 16 to 22
m5 = fullMatrix(:, 23:29);       % Columns 23 to 29

end

function subplotRB(RB,titlename)
% Assuming matrix3 is your 7-column matrix
% Extract quaternion part (last four columns)
quaternions = RB(:, 4:7);

% Convert quaternion to Euler angles (roll, pitch, yaw)
euler_angles = quat2eul(quaternions,'ZYX')*180/pi;

% Plotting
figure('Name', titlename);

% Plot the first three columns in the first row of subplots
for i = 1:3
    subplot(2, 3, i);
    plot(RB(:, i));
    title(['Position ' num2str(i)]);
end

% Plot Euler angles in the second row of subplots
for i = 1:3
    subplot(2, 3, i + 3);
    plot(euler_angles(:, i));
    title(['Euler Angle ' num2str(i)]);
end

end

function subplotIMU(IMU,titlename)
% Convert quaternion to Euler angles (roll, pitch, yaw)
euler_angles = IMU(:, 1:3);

% Plotting
figure('Name', titlename);

% Plot the first three columns in the first row of subplots
subplot(2, 3, 5);
plot(IMU(:, 4));
title('Position Angles');

% Plot Euler angles in the second row of subplots
for i = 1:3
    subplot(2, 3, i);
    plot(euler_angles(:, i));
    title(['Euler Angle ' num2str(i)]);
end

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

