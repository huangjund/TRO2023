clc;
clear all;
%close all;

syms d t alpha F delta_y miu

E = 37; %MPa 7.7
miu = 0.38;    %
t = 0.2*10^-3;
h0 = 4.8*10^-3;
D0 = 21.5*10^-3;
d = 13.9*10^-3;
N = 6;
l0 = sqrt((D0-d)^2/4+h0^2/4);

max_alpha = (2*l0+d)/d;

[x,y,force] = readMatrix();

% %% computation method in paper
% alpha = 1+sqrt(h0^2+(D0-d)^2-(y/N+h0)^2)/d;
% 
% ky1 = pi*E*t^3/(3*N*(1-miu^2)*d^2*(log(alpha)-(alpha-1)+(alpha-1)^2/2));
% F1 = ky1*y;
% 
% figure();
% fplot(F1,[-15,15]*10^-3);
% %% computation method in 1
% 
% ky2 = pi*E*t*d*(alpha-1)/(N*(alpha^4-1+log(alpha)-2*miu*alpha^2+2*miu)*2*l0);
% F2 = ky2*y;
% 
% figure();
% fplot(F2,[-15,15]*10^-3);

%% plot experimental data
% 绘制三维图形
figure();
plot3(repmat(x, size(y)), repmat(y, size(x)), -force, 'o-');
grid on;

% 设置图形标题和轴标签
title('Actuator Experiment');
xlabel('Initial Volumn/mm^3');
ylabel('Displacement/mm');
zlabel('Force/N');

%% computation method in 2
%syms r d D l0 E dtheta F
syms r dtheta F distance_y
I = t^3/(12*(1-miu^2))*r*dtheta;
D = d+sqrt(h0^2+(D0-d)^2-(distance_y/N+h0)^2);
f = F*dtheta/2/pi;
cosphi = (D-d)/(2*l0);
x0 = (r-d/2)/cosphi;
M = f*(r-d/2)+f*8/d^2*(cosphi/3*x0^3+d/4*x0^2);
dM = r-d/2+8/d^2*(cosphi/3*x0^3+d/4*x0^2);
G = simplify(2*M*dM/E/I);
%result = simplify(int(G,r,d/2,D/2));
delta_y = -2*(125*F*(12*miu^2 - 12)*((27*d^10)/2 - 72*D*d^9 + (315*D^2*d^8)/2 - 180*D^3*d^7 + (225*D^4*d^6)/2 - 36*D^5*d^5 + (9*D^6*d^4)/2 + (32*D^6*l0^4)/3 + (92*d^6*l0^4)/5 - 32*d^8*l0^2 - 9*d^10*log(d/D) + 88*D*d^7*l0^2 - (192*D^5*d*l0^4)/5 + 36*D*d^9*log(d/D) - 48*D^2*d^4*l0^4 - 44*D^2*d^6*l0^2 + (64*D^3*d^3*l0^4)/3 - 88*D^3*d^5*l0^2 + 36*D^4*d^2*l0^4 + 128*D^4*d^4*l0^2 - 64*D^5*d^3*l0^2 + 12*D^6*d^2*l0^2 - 54*D^2*d^8*log(d/D) + 36*D^3*d^7*log(d/D) - 9*D^4*d^6*log(d/D) - 16*d^6*l0^4*log(d/D) + 24*d^8*l0^2*log(d/D) - 48*D*d^7*l0^2*log(d/D) + 24*D^2*d^6*l0^2*log(d/D)))/(72*E*d^4*pi*(D - d)^4);

ky3 = F/(delta_y*N);
F3 = ky3*distance_y;

figure();
fplot(F3,[-15,15]*10^-3); hold on;
plot(y(17:67)*10^-3,-force(17:67,9)');

%% function
function [x,y,z] = readMatrix()

% 读取Excel数据
filename = 'data/actuator_data.xlsx';
sheet = 'Average';

% 读取x轴数据（第一行）
y = xlsread(filename, sheet, 'A2:A74');

% 读取y轴数据（第一列）
x = xlsread(filename, sheet, 'B1:R1');

% 读取z轴数据（第二列开始，除去第一行）
z = xlsread(filename, sheet, 'B2:R74');


end