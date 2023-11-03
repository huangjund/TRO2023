clear all;
clc;
%clf;

syms x2 q1 z2 x3 z3 theta3 q2 z4 q3 q5 q6

H23 = [eye(3) [-x2+q1;0;z2];[0 0 0 1]];
H34 = [[cos(theta3) 0 -sin(theta3);0 1 0; sin(theta3) 0 cos(theta3)]*...
    [1 0 0; 0 cos(q2) -sin(q2); 0 sin(q2) cos(q2)], [x3;0;z3];[0 0 0 1]];
H45 = [[0 -1 0; 1 0 0; 0 0 1]*...
    [1 0 0; 0 cos(q3) -sin(q3); 0 sin(q3) cos(q3)], [0;0;z4];[0 0 0 1]];

H27 = [eye(3) [x2-q1;0;z2];[0 0 0 1]];
H78 = [[cos(-theta3) 0 -sin(-theta3);0 1 0; sin(-theta3) 0 cos(-theta3)]*...
    [1 0 0; 0 cos(q5) -sin(q5); 0 sin(q5) cos(q5)], [-x3;0;z3];[0 0 0 1]];
H89 = [[0 -1 0; 1 0 0; 0 0 1]*...
    [1 0 0; 0 cos(q6) -sin(q6); 0 sin(q6) cos(q6)], [0;0;z4];[0 0 0 1]];

H25 = simplify(H23*H34*H45)
H29 = simplify(H27*H78*H89)

syms R23 x23 R34 x34 R45 x45 R27 x27 R78 x78 R89 x89
syms invR23 invR34 invR45 invR27 invR78 invR89

symH25 = simplify([R23 x23;0 1]*[R34 x34;0 1]*[R45 x45; 0 1])
symH29 = simplify([R27 x27;0 1]*[R78 x78;0 1]*[R89 x89; 0 1])


