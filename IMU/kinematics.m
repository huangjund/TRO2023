clear all;
clc;

syms q1 q2 q3 q4 q5 q6 q7

H23 = [[cos(-pi/6) 0 -sin(-pi/6); 0 1 0; sin(-pi/6) 0 cos(-pi/6)] [-59.07/2+q1; 0; 60.16];...
    0 0 0 1];
H34 = [[1 0 0; 0 cos(q2) -sin(q2); 0 sin(q2) cos(q2)] [0; 0; 0];...
    0 0 0 1];
H45 = [[cos(pi/2) -sin(pi/2) 0; sin(pi/2) cos(pi/2) 0; 0 0 1]*...
    [1 0 0; 0 cos(q3) -sin(q3); 0 sin(q3) cos(q3)] [0; 0; 9.56];...
    0 0 0 1];
H56 = [[1 0 0; 0 cos(q4) -sin(q4); 0 sin(q4) cos(q4)] [0;0;52.06];...
    0 0 0 1];

H27 = [[cos(pi/6) 0 -sin(pi/6); 0 1 0; sin(pi/6) 0 cos(pi/6)] [59.07/2-q1; 0; 60.16];...
    0 0 0 1];
H78 = [[1 0 0; 0 cos(q5) -sin(q5); 0 sin(q5) cos(q5)] [0; 0; 0];...
    0 0 0 1];
H89 = [[cos(pi/2) -sin(pi/2) 0; sin(pi/2) cos(pi/2) 0; 0 0 1]*...
    [1 0 0; 0 cos(q3) -sin(q6); 0 sin(q6) cos(q6)] [0; 0; 9.56];...
    0 0 0 1];
H910 = [[1 0 0; 0 cos(q7) -sin(q7); 0 sin(q7) cos(q7)] [0;0;52.06];...
    0 0 0 1];

H25 = simplify(H23*H34*H45);
H29 = simplify(H27*H78*H89);

se3mat1 = MatrixLog6(simplify(TransInv(H29)*H25)); %% SE3 -> se3 R4x4
V1 = se3ToVec(se3mat1); %% R6

%[theta1,w_hat1]=R2so3(H29(1:3,1:3)'*H25(1:3,1:3));


function [theta, w_hat] = R2so3(R)
    if trace(R) == -1
        theta = pi;
        w_hat = 1/sqrt(2*(1+R(3,3)))*[R(1,3);R(2,3);R(3,3)+1]; %% unit vector
        so3
    else
        theta = acos((trace(R)-1)/2);
        W_hat = (R-R')/(2*sin(theta));
        w_hat = [W_hat(3,2);W_hat(1,3);W_hat(2,1)];
    end
end