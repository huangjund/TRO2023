close all;

%% plot rotation matrix data3 process
% rotx = @(t) [1 0 0; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)] ;
% roty = @(t) [cos(t) 0 sin(t) ; 0 1 0 ; -sin(t) 0  cos(t)] ;
% rotz = @(t) [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1] ;
% 
% Rg1 = rotz(IMU2(81,1)/180*pi)*roty(IMU2(81,2)/180*pi)*rotx(IMU2(81,3)/180*pi);
% Rrb1 =rotz(euler_angles1(81,1)/180*pi)*roty(euler_angles1(81,2)/180*pi)*rotx(euler_angles1(81,3)/180*pi); 
% Rg2 = rotz(IMU2(61,1)/180*pi)*roty(IMU2(61,2)/180*pi)*rotx(IMU2(61,3)/180*pi);
% Rrb2 =rotz(euler_angles1(61,1)/180*pi)*roty(euler_angles1(61,2)/180*pi)*rotx(euler_angles1(61,3)/180*pi);
% 
% lhs = Rg1'*Rg2; rhs = Rrb1'*Rrb2;
% 
% % Define known matrices B and C (example matrices) A*B = C*A
% B = rhs; % Example B, must be a valid SO(3) matrix
% C = lhs; % Example C, must be a valid SO(3) matrix
% 
% % Objective function to minimize
% objectiveFunction = @(A) norm(reshape(A, [3, 3]) * B - C * reshape(A, [3, 3]), 'fro')*10;
% 
% % Constraint function to ensure A is in SO(3)
% constraintFunction = @(A) deal([], ... % Inequality constraints (none)
%     [norm(reshape(A, [3, 3])' * reshape(A, [3, 3]) - eye(3), 'fro'); ... % Orthogonality constraint
%      abs(det(reshape(A, [3, 3])) - 1)]); % Determinant constraint
% 
% % Initial guess for A
% initialGuess = randomSO3();
% initialGuess = initialGuess(:);
% 
% % Options for the optimizer
% options = optimoptions('fmincon', 'Algorithm', 'sqp', ...
%     'Display', 'iter', 'SpecifyObjectiveGradient', false, ...
%     'SpecifyConstraintGradient', false, 'CheckGradients', false,'MaxFunctionEvaluations',2000);
% 
% % Run the optimization
% [A_optimized, fval] = fmincon(objectiveFunction, initialGuess, ... % Objective function and initial guess
%     [], [], [], [], ... % Linear inequality and equality constraints (none)
%     [], [], ... % Bounds (none)
%     constraintFunction, ... % Nonlinear constraints
%     options); % Options
% 
% % Extract optimized matrix A
% A_optimized = reshape(A_optimized, [3, 3]);
% 
% % Display the optimized matrix A
% disp('Optimized matrix A:');
% disp(A_optimized);
% disp(reshape(initialGuess, [3, 3]))
% disp(A_optimized*B*A_optimized.'-C);
% 
% E = Rg1*A_optimized*Rrb1';
% 
% for i = 1:length(IMU2)
%     estimation_matrix(:,:,i) = E'*eul2rotm(IMU2(i,:)/180*pi,'ZYX')*A_optimized;
%     estimation_euler(i,:) = rotm2eul(estimation_matrix(:,:,i),'ZYX');
% end
% 
% % Plot Euler angles in the second row of subplots
% for i = 1:3
%     subplot(4, 3, i);
%     plot(IMU1(10:end, i));
%     title(['Euler Angle ' num2str(i)]);
% end
% for i = 1:3
%     subplot(4, 3, i+3);
%     plot(IMU2(10:end, i));
%     title(['Euler Angle ' num2str(i)]);
% end
% for i = 1:3
%     subplot(4, 3, i+6);
%     plot(euler_angles1(10:end, i)); hold on;
%     title(['Euler Angle ' num2str(i)]);
% end
% 
% for i = 1:3
%     subplot(4, 3, i+6);
%     plot(estimation_euler(10:end, i)*180/pi); grid on;
%     %title(['Estimated Euler Angle ' num2str(i)]);
% 
%     subplot(4, 3, i+9);
%     plot(estimation_euler(10:end, i)*180/pi-euler_angles1(10:end,i)); grid on;
%     title(['difference' num2str(i)]);
% end

%% plot rotation matrix data5 process
rotx = @(t) [1 0 0; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)] ;
roty = @(t) [cos(t) 0 sin(t) ; 0 1 0 ; -sin(t) 0  cos(t)] ;
rotz = @(t) [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1] ;

Rg1 = rotz(IMU1(25,1)/180*pi)*roty(IMU1(25,2)/180*pi)*rotx(IMU1(25,3)/180*pi);
Rrb1 =rotz(euler_angles1(25,1)/180*pi)*roty(euler_angles1(25,2)/180*pi)*rotx(euler_angles1(25,3)/180*pi); 
Rg2 = rotz(IMU1(41,1)/180*pi)*roty(IMU1(41,2)/180*pi)*rotx(IMU1(41,3)/180*pi);
Rrb2 =rotz(euler_angles1(41,1)/180*pi)*roty(euler_angles1(41,2)/180*pi)*rotx(euler_angles1(41,3)/180*pi);

lhs = Rg1'*Rg2; rhs = Rrb1'*Rrb2;

% Define known matrices B and C (example matrices) A*B = C*A
B = rhs; % Example B, must be a valid SO(3) matrix
C = lhs; % Example C, must be a valid SO(3) matrix

% Objective function to minimize
objectiveFunction = @(A) norm(reshape(A, [3, 3]) * B - C * reshape(A, [3, 3]), 'fro')*10;

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
    'SpecifyConstraintGradient', false, 'CheckGradients', false,'MaxFunctionEvaluations',2000);

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

E = Rg1*A_optimized*Rrb1';

for i = 1:length(IMU1)
    estimation_matrix(:,:,i) = E'*eul2rotm(IMU1(i,:)/180*pi,'ZYX')*A_optimized;
    estimation_euler(i,:) = rotm2eul(estimation_matrix(:,:,i),'ZYX');
end

% Plot Euler angles in the second row of subplots
for i = 1:3
    subplot(4, 3, i);
    plot(IMU1(10:end, i));
    title(['Euler Angle ' num2str(i)]);
end
for i = 1:3
    subplot(4, 3, i+3);
    plot(IMU2(10:end, i));
    title(['Euler Angle ' num2str(i)]);
end
for i = 1:3
    subplot(4, 3, i+6);
    plot(euler_angles1(10:end, i)); hold on;
    title(['Euler Angle ' num2str(i)]);
end

for i = 1:3
    subplot(4, 3, i+6);
    plot(estimation_euler(10:end, i)*180/pi); grid on;
    %title(['Estimated Euler Angle ' num2str(i)]);

    subplot(4, 3, i+9);
    plot(estimation_euler(10:end, i)*180/pi-euler_angles1(10:end,i)); grid on;
    title(['difference' num2str(i)]);
end

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