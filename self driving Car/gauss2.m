%%
% ___________________________________________________
% Author : benjamin_zhao
% date :   2012-08-20
% The parameters for gauss fi^2=sig2 and u 
% listed in order in the vector [sig2 u].
% ___________________________________________________
%%
function f = gauss2(x, sig2, u)

f = 1/sqrt(2 * pi * sig2) * gaussmf(x,[sqrt(sig2) u]);