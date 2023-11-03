clear all;
clc;
%clf;

syms H12;

se3mat = MatrixLog6(H12); %% T -> [S]theta
V = se3ToVec(se3mat); %% [S]theta -> S*thetad
%% [S]theta -> S*theta
[S, theta] = AxisAng(SE3vec6); %% S*theta -> S, theta