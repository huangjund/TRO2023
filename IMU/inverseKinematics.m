clear all;
clc;

syms q1 q2 q3 q4 q5 q6 q7 real
rotx = @(t) [1 0 0; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)] ;
rotz = @(t) [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1] ;
roty = @(t) [cos(t) 0 sin(t) ; 0 1 0 ; -sin(t) 0  cos(t)] ;

xrange = linspace(0,pi/2,40);
yrange = linspace(0,pi/2,40);
[X, Y] = meshgrid(xrange, yrange);

R25_fake = rotx(q2)*rotz(pi/2)*rotx(q3);
R25_omg = simplify(so3ToVec(MatrixLog3(R25_fake)));
[R25_omghat, R25_theta] = AxisAng3(R25_omg);

for x = 1:length(xrange)
    for y = 1:length(yrange)
        R25_omg_array((x-1)*length(yrange)+y,:) = simplify(vpa(subs(R25_omg,{q2,q3},{xrange(x),yrange(y)})));
    end
end
for i = 1:length(R25_omg_array)
    [R25_omghat_array(i,:), R25_theta(i)] = AxisAng3(R25_omg_array(i,:));
end

%fplot3(double(subs(R25_omg(1), {q2, q3}, {X, Y})),double(subs(R25_omg(2), {q2, q3}, {X, Y})),double(subs(R25_omg(3), {q2, q3}, {X, Y})))