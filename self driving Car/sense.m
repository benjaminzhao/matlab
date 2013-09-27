%%
% ___________________________________________________
% Author : benjamin_zhao
% date :   2012-08-20
% sense functions
% find the probability distribution of the world
% with one measurement
% input:
%   p:     current probablity distribution -- data array
%   world: world description sequence -- cell array
%   measurement: one measurement -- cell
%   hit: hit probabilty -- number
%   miss: miss probability -- number
% output:
%   new probability distribution
% ___________________________________________________
%%
function q = sense(p, world, measurement, hit, miss)

[m,n] = size(p); % get p size = m*n
q = zeros(m,n);  %allocate q as new p for speed in the loop

[row, col] = size(world); % get world size = m*n
%% check
if (m ~= row)||(n ~= col)
    %errargt(mfilename,'p & world is unequal in size','msg');
    error('probability & world are unequal in size');
    return;
end

%% sense one measurement
for i = 1:row
    for j = 1:col
        ishit = isequal(world{i,j}, measurement);
        q(i,j) = p(i,j) * (ishit * hit + (1 - ishit) * miss);
    end
end
qsum = sum(sum(q));
q = q / qsum; % normalize


