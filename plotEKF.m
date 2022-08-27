clear; clc; close all;

T = csvread('EKF.csv');

time = T(:,1);
x_true = T(:,2:5);
x_dr = T(:,6:9);
x_Est = T(:,10:13);
z = T(:,14:15);


%% plot

figure(100);
plot(x_true(:,1), x_true(:,2), '-r', 'LineWidth', 2)
axis equal
hold on
plot(x_Est(:,1), x_Est(:,2), '--b')
hold on
plot(x_dr(:,1), x_dr(:,2), '-.k')
grid on;