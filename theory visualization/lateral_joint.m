clear all;
clc;
clf;

p3 = [0 1 0 1];
p1 = [0.5 1 0 1];
p2 = [-0.5 1 0 1];

syms q2 q3 x0 z0;
H34 = [[1 0 0; 0 cos(q2) -sin(q2); 0 sin(q2) cos(q2)] [0;0;0];[0 0 0 1]];

H45 = [[cos(pi/2) -sin(pi/2) 0; sin(pi/2) cos(pi/2) 0; 0 0 1]*...
    [1 0 0; 0 cos(q3) -sin(q3); 0 sin(q3) cos(q3)] [x0; 0; z0];[0 0 0 1]];
H = H34*H45;



%% plot1 translate z0
p1_w = H*p1';
p2_w = H*p2';
p3_w = H*p3';
pz_w = H(1:3,1:3)*[0;0;1];

start_p1 = [-1 0.5 0];
start_p2 = [-1 -0.5 0];

figure(1);
p1_wsub = subs(p1_w,[x0, z0],[0, 0.2]);
fsurf(p1_wsub(1),p1_wsub(2),p1_wsub(3),[-pi/6,pi/6,0,pi/2]); hold on;
p2_wsub = subs(p2_w,[x0, z0],[0, 0.2]);
fsurf(p2_wsub(1),p2_wsub(2),p2_wsub(3),[-pi/6,pi/6,0,pi/2]); hold on;


theta1 = -pi/6:0.1:pi/6;
theta2 = 0:0.1:pi/2;

for i = theta1
    for j = theta2
        p1_end = double(subs(p1_wsub,[q2,q3],[i,j]));
        p2_end = double(subs(p2_wsub,[q2,q3],[i,j]));
        visualize_tangents(start_p1,p1_end(1:3)',[0,0,1]);
        % Plot the tangent at the end point
        pz_wsub = subs(pz_w,[x0, z0,q2,q3],[0, 0.2 i j]);
        quiver3(p1_end(1), p1_end(2), p1_end(3), ...
                pz_wsub(1), pz_wsub(2), pz_wsub(3), 'r');
%         visualize_tangents(start_p2,p2_end(1:3)',[0,0,1]);
    end
end

xlabel('x', 'FontSize', 16, 'FontWeight','bold', 'Interpreter', 'latex', 'rotation',15);
ylabel('y', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex', 'rotation',-15);
zlabel('z', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');

hold off;

%% plot 2 translate x0
% figure(2);
% for i = 1:1
%     p1_w = H*(p1-[0 i/10 0 0])';
%     p2_w = H*(p2-[0 i/10 0 0])';
%     p3_w = H*p3';
%     p1_wsub = subs(p1_w,[z0, x0],[0.2, -i/10]);
%     fsurf(p1_wsub(1),p1_wsub(2),p1_wsub(3),[-pi/6,pi/6,0,pi/2]); hold on;
%     p2_wsub = subs(p2_w,[z0, x0],[0.2, -i/10]);
%     fsurf(p2_wsub(1),p2_wsub(2),p2_wsub(3),[-pi/6,pi/6,0,pi/2]); hold on;
% end
% xlabel('x', 'FontSize', 16, 'FontWeight','bold', 'Interpreter', 'latex', 'rotation',15);
% ylabel('y', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex', 'rotation',-15);
% zlabel('z', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');
% 


%% functions
function visualize_tangents(start_point, end_points, start_tangent)
    % Number of end points
    num_end_points = size(end_points, 1);

    % Plot the start and end points
    plot3(start_point(1), start_point(2), start_point(3), 'ro');
    plot3(end_points(:, 1), end_points(:, 2), end_points(:, 3), 'bo');
    
    % Plot the tangent at the start point
    quiver3(start_point(1), start_point(2), start_point(3), ...
            start_tangent(1), start_tangent(2), start_tangent(3), 'g');

    for i = 1:num_end_points
        C = circle_center(start_point,end_points(i,:),start_tangent);
        
        end_tangent = tangent_at_P2(start_point,end_points(i,:),start_tangent,C);
        
        % Plot the tangent at the end point
        quiver3(end_points(i, 1), end_points(i, 2), end_points(i, 3), ...
                end_tangent(1), end_tangent(2), end_tangent(3), 'b');
    end

    % Set plot properties
    axis equal;
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