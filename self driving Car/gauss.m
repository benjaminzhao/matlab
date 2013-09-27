%%
% ___________________________________________________
% Author : benjamin_zhao
% date :   2012-08-20
% The parameters for gauss fi=sig and u 
% listed in order in the vector [sig u].
% ___________________________________________________
%%
function f = gauss(x, sig, u)

f = 1/sqrt(2 * pi * sig^2) * gaussmf(x,[sig u]);

