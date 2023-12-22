clear all;
clc;
close all;

[IMU1, RB1] = readMatrix('./data/20231123SUMdata5.csv');

% %% data filter
% RB1(:,4) = LowPassFilter(RB1(:,4));
% RB1(:,5) = LowPassFilter(RB1(:,5));
% RB1(:,6) = LowPassFilter(RB1(:,6));
% RB1(:,7) = LowPassFilter(RB1(:,7));

%% plot
% Plotting imu
figure('Name', 'IMU1');

% Plot Euler angles in the second row of subplots
for i = 1:3
    subplot(4, 3, i);
    plot(IMU1(:, i)); grid on;
    title(['Euler Angle ' num2str(i)]);
end
% plot Rigidbody
quaternions = RB1(:, 4:7);

% Convert quaternion to Euler angles (roll, pitch, yaw)
euler_angles = quat2eul(quaternions,'ZYX')*180/pi;


% Plot Euler angles in the second row of subplots
for i = 1:3
    subplot(4, 3, i + 3);
    plot(euler_angles(:, i));grid on;hold on;
    title(['Euler Angle ' num2str(i)]);
end



%% plot rotation matrix
rotx = @(t) [1 0 0; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)] ;
roty = @(t) [cos(t) 0 sin(t) ; 0 1 0 ; -sin(t) 0  cos(t)] ;
rotz = @(t) [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1] ;

% R_imu1 = [];
% R_rb1 = [];
% R_rb1_imu1 = [];
% for i = 1:length(IMU1)
%     R_imu1(:,:,i) = rotz(IMU1(i,1)/180*pi)...
%         *roty(IMU1(i,2)/180*pi)...
%         *rotx(IMU1(i,3)/180*pi);
%     R_rb1(:,:,i) = rotz(euler_angles(i,1)/180*pi)...
%         *roty(euler_angles(i,2)/180*pi)...
%         *rotx(euler_angles(i,3)/180*pi);
% 
%     %R_rb2_temp = R_rb2(:,:,i); R_rb3_temp = R_rb3(:,:,i);
%     R_rb1_imu1(:,:,i) = R_imu1(:,:,i)*R_rb1(:,:,i)';
% 
%     Euler_rb1_imu1(i,:) = rotm2eul(R_rb1_imu1(:,:,i),'ZYX');
% end
% 
% subplot(3,3,7); plot((Euler_rb1_imu1(:,1)/pi*180));grid on;
% hold on;
% subplot(3,3,8); plot((Euler_rb1_imu1(:,2)/pi*180));grid on;
% hold on;
% subplot(3,3,9); plot((Euler_rb1_imu1(:,3)/pi*180));
% hold on;grid on;

% subplot(3,3,7); plot(LowPassFilter(Euler_rb1_imu1(:,1)/pi*180));
% subplot(3,3,8); plot(LowPassFilter(Euler_rb1_imu1(:,2)/pi*180));
% subplot(3,3,9); plot(LowPassFilter(Euler_rb1_imu1(:,3)/pi*180));

%% test
% solution below: A = [0.7394 0.3867 0.5509; -0.4898 -0.2519 0.8346;0.4618 -0.8871 0.0031];
Rg1 = rotz(IMU1(187,1)/180*pi)*roty(IMU1(187,2)/180*pi)*rotx(IMU1(187,3)/180*pi);
Rrb1 =rotz(euler_angles(187,1)/180*pi)*roty(euler_angles(187,2)/180*pi)*rotx(euler_angles(187,3)/180*pi); 
Rg2 = rotz(IMU1(234,1)/180*pi)*roty(IMU1(234,2)/180*pi)*rotx(IMU1(234,3)/180*pi);
Rrb2 =rotz(euler_angles(234,1)/180*pi)*roty(euler_angles(234,2)/180*pi)*rotx(euler_angles(234,3)/180*pi);

% solution below: A = [-0.1386 -0.8 -0.585; 0.3655 0.5074 -0.7804; 0.9204 -0.3218 0.2219];
% Rg1 = rotz(IMU1(649,1)/180*pi)*roty(IMU1(649,2)/180*pi)*rotx(IMU1(649,3)/180*pi);
% Rrb1 =rotz(euler_angles(649,1)/180*pi)*roty(euler_angles(649,2)/180*pi)*rotx(euler_angles(649,3)/180*pi); 
% Rg2 = rotz(IMU1(736,1)/180*pi)*roty(IMU1(736,2)/180*pi)*rotx(IMU1(736,3)/180*pi);
% Rrb2 =rotz(euler_angles(736,1)/180*pi)*roty(euler_angles(736,2)/180*pi)*rotx(euler_angles(736,3)/180*pi);

lhs = Rg1'*Rg2; rhs = Rrb1'*Rrb2;

% Define known matrices B and C (example matrices) A*B = C*A
B = rhs; % Example B, must be a valid SO(3) matrix
C = lhs; % Example C, must be a valid SO(3) matrix

% Objective function to minimize
objectiveFunction = @(A) norm(reshape(A, [3, 3]) * B - C * reshape(A, [3, 3]), 'fro');

% Constraint function to ensure A is in SO(3)
constraintFunction = @(A) deal([], ... % Inequality constraints (none)
    [norm(reshape(A, [3, 3])' * reshape(A, [3, 3]) - eye(3), 'fro'); ... % Orthogonality constraint
     abs(det(reshape(A, [3, 3])) - 1)]); % Determinant constraint

% Initial guess for A
initialGuess = randomSO3();
initialGuess = initialGuess(:);

% Options for the optimizer
options = optimoptions('fmincon', 'Algorithm', 'sqp', ...
    'Display', 'iter', 'SpecifyObjectiveGradient', false, ...
    'SpecifyConstraintGradient', false, 'CheckGradients', false);

% Run the optimization
[A_optimized, fval] = fmincon(objectiveFunction, initialGuess, ... % Objective function and initial guess
    [], [], [], [], ... % Linear inequality and equality constraints (none)
    [], [], ... % Bounds (none)
    constraintFunction, ... % Nonlinear constraints
    options); % Options

% Extract optimized matrix A
A_optimized = reshape(A_optimized, [3, 3]);

% Display the optimized matrix A
disp('Optimized matrix A:');
disp(A_optimized);
disp(reshape(initialGuess, [3, 3]))
disp(A_optimized*B*A_optimized.'-C);

%% plot corrected
A = [0.7394 0.3867 0.5509; -0.4898 -0.2519 0.8346;0.4618 -0.8871 0.0031];
% A = [-0.1386 -0.8 -0.585; 0.3655 0.5074 -0.7804; 0.9204 -0.3218 0.2219];
%A = [0.234699733668458 -0.207986430826956 -0.949556712668646; 0.598790006344904 0.800439598704894 -0.027322841777775; 0.765745641617386 -0.562172576353252 0.312403229644010];
E = Rg1*A*Rrb1';

for i = 1:length(IMU1)
    estimation_matrix(:,:,i) = E'*eul2rotm(IMU1(i,:)/180*pi,'ZYX')*A;
    estimation_euler(i,:) = rotm2eul(estimation_matrix(:,:,i),'ZYX');
end

for i = 1:3
    subplot(4, 3, i+3);
    plot(estimation_euler(:, i)*180/pi); grid on;
    %title(['Estimated Euler Angle ' num2str(i)]);

    subplot(4, 3, i+9);
    plot(estimation_euler(:, i)*180/pi-euler_angles(:,i)); grid on;
    title(['difference' num2str(i)]);
end

%% functions

function m = randomSO3()
% Generate a random rotation matrix in SO(3)

% Step 1: Create a random rotation axis
axis = randn(3, 1);
axis = axis / norm(axis); % Normalize the axis

% Step 2: Generate a random rotation angle
angle = 2 * pi * rand; % Random angle between 0 and 2π

% Step 3: Create the corresponding rotation matrix using Rodrigues' rotation formula
K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0]; % Skew-symmetric matrix
m = eye(3) + sin(angle) * K + (1 - cos(angle)) * K^2; % Rotation matrix

end

function m = v2matrix(v)
    m = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

function [m1,m3] = readMatrix(filename)
% Read the matrix from the CSV file
fullMatrix = readmatrix(filename);

% Check if the matrix has enough columns
totalColumnsNeeded = 3 + 7;
if size(fullMatrix, 2) < totalColumnsNeeded
    error('The matrix does not have enough columns.');
end

% Divide the matrix into five matrices
m1 = fullMatrix(500:end, 1:3);         % Columns 1 to 4
m3 = fullMatrix(500:end, 4:10);        % Columns 9 to 15

end

function subplotRB(RB,titlename, mode)
% Assuming matrix3 is your 7-column matrix
% Extract quaternion part (last four columns)
quaternions = RB(:, 4:7);

% Convert quaternion to Euler angles (roll, pitch, yaw)
euler_angles = quat2eul(quaternions,mode)*180/pi;

% Plotting
figure('Name', titlename);
% 
% % Plot the first three columns in the first row of subplots
% for i = 1:3
%     subplot(2, 3, i);
%     plot(RB(:, i));
%     title(['Position ' num2str(i)]);
% end

% Plot Euler angles in the second row of subplots
for i = 1:3
    subplot(1, 3, i);
    plot(euler_angles(:, i));
    title(['Euler Angle ' num2str(i)]);
end

end

function subplotIMU(IMU,titlename)
% Convert quaternion to Euler angles (roll, pitch, yaw)
euler_angles = IMU(:, 1:3);

% Plotting
figure('Name', titlename);

% Plot Euler angles in the second row of subplots
for i = 1:3
    subplot(1, 3, i);
    plot(euler_angles(:, i));
    title(['Euler Angle ' num2str(i)]);
end

end

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
