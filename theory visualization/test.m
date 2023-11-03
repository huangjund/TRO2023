% Define circle_center and tangent_at_P2 as before...

% P1 and T1
P1 = [1, 2, 3];
T1 = [0, 0, 1];

% Generate a set of P2 points (example: random points)
num_points = 20;
P2_set = rand(3, num_points) * 10;  % 3x20 matrix, each column is a P2

% Initialize storage for T2 vectors
T2_set = zeros(3, num_points);

% Compute T2 for each P2
for i = 1:num_points
    P2 = P2_set(:, i);
    C = circle_center(P1', P2, T1');  % Assuming column vectors in circle_center
    T2_set(:, i) = tangent_at_P2(P1', P2, T1', C);  % Store each T2
end

% Plot
figure; hold on;

% Plot each P2 and its corresponding T2 using quiver3
for i = 1:num_points
    P2 = P2_set(:, i);
    T2 = T2_set(:, i);
    quiver3(P2(1), P2(2), P2(3), T2(1), T2(2), T2(3), 1, 'r');  % '1' scales the arrow, 'r' specifies red color
end

title('Tangent Vectors T2 at Different P2 Points');
xlabel('x'); ylabel('y'); zlabel('z');
axis equal;
grid on;

% If you want to save the figure for your paper:
print('TangentVectors','-dpng','-r300');  % Save as PNG with 300 dpi

function visualize_tangents(start_point, end_points, start_tangent)
    % Number of end points
    num_end_points = size(end_points, 1);

    % Plot the start and end points
    figure; hold on;
    plot3(start_point(1), start_point(2), start_point(3), 'ro');
    plot3(end_points(:, 1), end_points(:, 2), end_points(:, 3), 'bo');
    
    % Plot the tangent at the start point
    quiver3(start_point(1), start_point(2), start_point(3), ...
            start_tangent(1), start_tangent(2), start_tangent(3), 'r');

    for i = 1:num_end_points
        C = circle_center(start_point,end_points(i,:),start_tangent);
        
        end_tangent = tangent_at_P2(start_point,end_points(i,:),start_tangent,C);
        
        % Plot the tangent at the end point
        quiver3(end_points(i, 1), end_points(i, 2), end_points(i, 3), ...
                end_tangent(1), end_tangent(2), end_tangent(3), 'b');
    end

    % Set plot properties
    grid on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Visualization of Tangent Vectors at End Points');
end

function C = circle_center(P1, P2, T1)

    % Calculate normal to the plane
    N = cross(T1,P2 - P1);

    % Find the midpoint
    M = 0.5 * (P1 + P2);

    % Bisector of the segment joining P1 and P2
    B = cross(N, P2 - P1);

    % Find the angle between T1 and B
    theta = acos(dot(T1, B) / (norm(T1) * norm(B)))-pi/2;

    % Radius of the circle
    r = norm(P2 - P1) / (2 * sin(theta));

    % Center of the circle
    C = M + r * (B / norm(B))*cos(theta);

end

function T2 = tangent_at_P2(P1, P2, T1, C)
    % Calculate normal to the plane
    N = cross(P2 - P1, T1);

    % Radius vector at P2
    R = P2 - C;

    % Tangent vector at P2
    T2 = -cross(N, R);
    
    % Normalize the tangent vector
    T2 = T2 / norm(T2);
end
