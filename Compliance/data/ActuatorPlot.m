% 清除之前的变量和命令窗口
clear;
clc;

% 读取Excel数据
filename = 'actuator_data.xlsx';
sheet = 'Average';

% 读取x轴数据（第一行）
y = xlsread(filename, sheet, 'A2:A74');

% 读取y轴数据（第一列）
x = xlsread(filename, sheet, 'B1:R1');

% 读取z轴数据（第二列开始，除去第一行）
z = xlsread(filename, sheet, 'B2:R74');



% 绘制三维图形
plot3(repmat(x, size(y)), repmat(y, size(x)), z, 'o-');
grid on;
% surf(x,y,z);


% 设置图形标题和轴标签
title('Actuator Experiment');
xlabel('Pump Distance/mm');
ylabel('Screw Distance/mm');
zlabel('Force/N');