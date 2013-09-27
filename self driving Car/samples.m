%%
% ___________________________________________________
% Author : benjamin_zhao
% date :   2013-05-22
% produce sample data array
% 
% ___________________________________________________
function p = samples(Ze, sig, n)
p = normrnd(Ze, sig, 1,n);

% clc
% clear;
% Ze = 10;
% n = 10;
% sig = 0.5;
% Z = normrnd(Ze, sig, n,1);
% figure, hist(Z, 8:0.1:12);
% p = mvnpdf(Z,Ze,sig);
% figure, hist(p);

