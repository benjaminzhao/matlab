%%
% ___________________________________________________
% Author : benjamin_zhao
% date :   2012-08-20
% 2D move functions
% find the probability distribution of the world
% with any measurements sequence
% input:
%   pm: current probablity distribution -- data array
%   motion: world description sequence -- cell
%   pExact: probabilty of move to exact location -- number
%   pOvershoot: probability of move overshoot -- number
%   pUndershoot: probability of move undershoot -- number
% output:
%   new probability distribution
% ___________________________________________________
%%
function r = move(pm, motionSTR, pExact, pOvershoot, pUndershoot)

[row, col] = size(pm);
r = zeros(row, col);

%% move translation
move = {[0,0],[0,1],[0,-1],[1,0],[-1,0]};
move_name = {'no','right','left','down','up'};
move_label = {'x','>','<','v','^'};
for i = 1:length(move);
    if isequal(motionSTR, move_name{i})
        action = i;
        break
    end
end
motion = move{action};
motion_print = move_label(action);
disp(motion_print);

%%
for i = 1:row
    for j = 1:col
        tempROW = mod(i-motion(1,1),row);
        tempCOL = mod(j-motion(1,2),col);
        if tempROW == 0
            tempROW = row;
        end
        if tempCOL == 0
            tempCOL = col;
        end
        r(i,j) = r(i,j) + pExact * pm(tempROW, tempCOL) + pUndershoot * pm(i, j);
    end
end

%%
%for i = 1:row
%    for j = 1:col
%         tempROW = mod(i-motion(1,1),row);
%         tempCOL = mod(j-motion(1,2),col);
%         if tempROW == 0
%             tempROW = row;
%         end
%         if tempCOL == 0
%             tempCOL = col;
%         end
%         r(i,j) = r(i,j) + pExact * pm(tempROW, tempCOL);
% 
%         tempROW = mod(i-(motion(1)+1),row);
%         tempCOL = mod(j-(motion(2)+1),col);
%         if tempROW == 0
%             tempROW = row;
%         end
%         if tempCOL == 0
%             tempCOL = col;
%         end
%         r(i,j) = r(i,j) + pOvershoot * pm(tempROW, tempCOL);
% 
%         tempROW = mod(i-(motion(1)-1),row);
%         tempCOL = mod(j-(motion(2)-1),col);
%         if tempROW == 0
%             tempROW = row;
%         end
%         if tempCOL == 0
%             tempCOL = col;
%         end
%         r(i,j) = r(i,j) + pUndershoot * pm(tempROW, tempCOL);
%    end
%end

