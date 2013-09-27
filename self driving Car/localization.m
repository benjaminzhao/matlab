%%
% ___________________________________________________
% Author : benjamin_zhao
% date :   2012-08-20
% 
% localization functions
% ___________________________________________________
%%
clear;
clc;
%%
% p1D = [0.2, 0.2, 0.2, 0.2, 0.2];
% %p = [0, 1, 0, 0, 0];
% world1D = {'green', 'red', 'red', 'green', 'green'};
% measurements1D = {'red','red'};
% motions1D = {[0,1],[0,1]};
% 
% pHit = 0.6; % probability of hit
% pMiss = 0.2; % probability of miss
% pExact = 0.8; % probabilty of move to exact location
% pOvershoot = 0.5*(1-pExact); % probability of move overshoot
% pUndershoot = 0.5*(1-pExact); % probability of move undershoot
% 
% for i = 1:length(measurements1D)
%     p1D = sense(p1D, world1D, measurements1D{i}, pHit, pMiss);
%     p1D = move(p1D, motions1D{i}, pExact, pOvershoot, pUndershoot);
% end

%%
world2D = { 'red', 'green', 'green', 'red',   'red';
            'red', 'red',   'green', 'red',   'red';
            'red', 'red',   'green', 'green', 'red';
            'red', 'red',   'red',   'red',   'red',};
[m,n] = size(world2D);
p2D = ones(m,n)/sum(sum(ones(m,n)));
measurements2D = {'green',  'green', 'green', 'green', 'green'};
motions2D = {'no', 'right', 'down', 'down', 'right'};

sensor_right = 0.7;
sensor_wrong = 1.0 - sensor_right;
p_move = 0.8;
p_stay = 1.0 - p_move;
p_over = 0;

for i = 1:length(measurements2D)
    p2D = move(p2D, motions2D{i}, p_move, p_over, p_stay);
    p2D = sense(p2D, world2D, measurements2D{i}, sensor_right, sensor_wrong);
    %p2D = move(p2D, motions2D{i}, p_move, p_over, p_stay);
end



%%
% figure;
% plot(p1D);
mesh(p2D);


