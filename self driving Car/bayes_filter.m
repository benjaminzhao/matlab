%%
% ___________________________________________________
% Author : benjamin_zhao
% date :   2013-05-20
% beyes filter algorithm
% example of a robot measures a door wether open
% input:
%   p: probablity distribution -- data table
%   bel_predict: belief table in prediction @ step(t)
%   bel: belief table @ step(t)
%   world X(t): state data -- cell array
%   measurement: Z(t) one measurement data at step(t) -- cell
%   control: U(t) control at step(t) -- numbel
%   
% output:
%   new probability distribution
% ___________________________________________________
%%
% function belx = bayes_filter(belx, Ut, Zt)
% 
function V = bayes_filter()
U = {'NULL', 'nothing', 'push', 'nothing', 'nothing'};    % U=control sequence
Z = {'open','open', 'open', 'open', 'open'};         % Z=measure sequence
%%X = ;   %% state sequence

X_name = {'open', 'close'}; % state type
Ut = {'nothing', 'push'};   % control type
Xt_1 = X_name;   % state type @ step(t-1)
Xt = X_name;     % state type @ step(t)

P = zeros(2, 2, 2); % prediction probabilty table with control
P(find_i(Xt, 'open'), find_i(Ut, 'nothing'), find_i(Xt_1, 'open')) = 1;
P(find_i(Xt, 'close'), find_i(Ut, 'nothing'), find_i(Xt_1, 'open')) = 0;
P(find_i(Xt, 'open'), find_i(Ut, 'nothing'), find_i(Xt_1, 'close')) = 0;
P(find_i(Xt, 'close'), find_i(Ut, 'nothing'), find_i(Xt_1, 'close')) = 1;
P(find_i(Xt, 'open'), find_i(Ut, 'push'), find_i(Xt_1, 'open')) = 1;
P(find_i(Xt, 'close'), find_i(Ut, 'push'), find_i(Xt_1, 'open')) = 0;
P(find_i(Xt, 'open'), find_i(Ut, 'push'), find_i(Xt_1, 'close')) = 0.8;
P(find_i(Xt, 'close'), find_i(Ut, 'push'), find_i(Xt_1, 'close')) = 0.2;

Z_name = {'open', 'close'}; % measure type
P_Z = zeros(2, 2);  % probabilty table with measurement
P_Z(find_i(Z_name, 'open'), find_i(X_name, 'open')) = 0.6;
P_Z(find_i(Z_name, 'close'), find_i(X_name, 'open')) = 0.4;
P_Z(find_i(Z_name, 'open'), find_i(X_name, 'close')) = 0.2;
P_Z(find_i(Z_name, 'close'), find_i(X_name, 'close')) = 0.8;

bel_name = {'open', 'close'};   % belief type
bel_predict = zeros(length(Xt), length(U));      % belief table in prediction @ step(t)
bel_temp = zeros(length(Xt), length(U));         % buffer for normalization
bel = zeros(length(Xt), length(U));              % belief table @ step(t)
bel(find_i(bel_name, 'open'), 1) = 0.5; %initial belief
bel(find_i(bel_name, 'close'), 1) = 0.5;%initial belief


for t = 2:length(U) % step(1)is initial state, calc from step(2)
%%   step(t) bel_prediction calc
%     bel_predict(find_i(bel_name, 'open'), t) = P(find_i(Xt, 'open'), find_i(Ut, U{t}), find_i(Xt_1, 'open')) * bel(find_i(bel_name, 'open'), t-1)...
%                                              + P(find_i(Xt, 'open'), find_i(Ut, U{t}), find_i(Xt_1, 'close')) * bel(find_i(bel_name, 'close'), t-1);
%     bel_predict(find_i(bel_name, 'close'), t) = P(find_i(Xt, 'close'), find_i(Ut, U{t}), find_i(Xt_1, 'open')) * bel(find_i(bel_name, 'open'), t-1)...
%                                               + P(find_i(Xt, 'close'), find_i(Ut, U{t}), find_i(Xt_1, 'close')) * bel(find_i(bel_name, 'close'), t-1);
%     code as below1:
%     for type = 1:length(Xt)
%         bel_predict(find_i(bel_name, 'open'), t) = P(find_i(Xt, 'open'), find_i(Ut, U{t}), find_i(Xt_1, Xt{type})) * bel(find_i(bel_name, Xt{type}), t-1)...
%                                                  + bel_predict(find_i(bel_name, 'open'), t);
%     end
%     for type = 1:length(Xt)
%         bel_predict(find_i(bel_name, 'close'), t) = P(find_i(Xt, 'close'), find_i(Ut, U{t}), find_i(Xt_1, Xt{type})) * bel(find_i(bel_name, Xt{type}), t-1) + bel_predict(find_i(bel_name, 'close'), t);
%     end
%   code as below2:
    for type1 = 1:length(Xt)
        for type2 = 1:length(Xt)
            bel_predict(find_i(bel_name, Xt{type1}), t) = P(find_i(Xt, Xt{type1}), find_i(Ut, U{t}), find_i(Xt_1, Xt{type2})) * bel(find_i(bel_name, Xt{type2}), t-1)...
                                                        + bel_predict(find_i(bel_name, Xt{type1}), t);
        end
    end
%%  step(t) bel calc
%    bel_temp(find_i(bel_name, 'open'), t) = P_Z(find_i(Z_name, Z{t}), find_i(X_name, 'open')) * bel_predict(find_i(bel_name, 'open'), t);
%    bel_temp(find_i(bel_name, 'close'), t) = P_Z(find_i(Z_name, Z{t}), find_i(X_name, 'close')) * bel_predict(find_i(bel_name, 'close'), t);
%    sum = bel_temp(find_i(bel_name, 'open'), t) + bel_temp(find_i(bel_name, 'close'), t);
%    bel(find_i(bel_name, 'open'), t) = bel_temp(find_i(bel_name, 'open'), t) / sum;
%    bel(find_i(bel_name, 'close'), t) = bel_temp(find_i(bel_name, 'close'), t) / sum;
    sum = 0;
    for type3 = 1:length(Xt)
        bel_temp(find_i(bel_name, Xt{type3}), t) = P_Z(find_i(Z_name, Z{t}), find_i(X_name, Xt{type3})) * bel_predict(find_i(bel_name, Xt{type3}), t);
        sum = sum + bel_temp(find_i(bel_name, Xt{type3}), t);
    end
    for type3 = 1:length(Xt)
        bel(find_i(bel_name, Xt{type3}), t) = bel_temp(find_i(bel_name, Xt{type3}), t) / sum;
    end

end
plot(1:length(U), bel(1,1:length(U)), 1:length(U), bel(2,1:length(U)));
V=bel;
