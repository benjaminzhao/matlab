%%
% ___________________________________________________
% Author : benjamin_zhao
% date :   2013-05-22
% The beam range finder model algorithm
% 
% ___________________________________________________
%%
function q = beam_model()
clc;
clear;
u1 = 0;
u2 = 20;
sig1 = 2;
sig2 = 3.5;
K = 100;
Z1 = samples(u1, sig1, K);
Z2 = samples(u2, sig2, K);
Z = [Z1, Z2];
%x = mean(Z);

% init value
alpha = [0.2, 0.8];
mu = [max(Z),min(Z)];
d = 10*rand();
sigma = [d, d];
e = 0.1;
maxIteration = 1000;
loglikelihood = 0;
oldlikelihood = -Inf;

P1t = zeros(1,length(Z));
P2t = zeros(1,length(Z));
W1t = zeros(1,length(Z));
W2t = zeros(1,length(Z));


for iter = 1:maxIteration
    %E step
    for i = 1 : length(Z)
        %P1t(i) = gauss2(Z(i), sigma(1), mu(1)) * alpha(1);
        %P2t(i) = gauss2(Z(i), sigma(2), mu(2)) * alpha(2);
        aaa = mvnpdf(Z(i), mu(1), sigma(1));
        bbb = mvnpdf(Z(i), mu(2), sigma(2));
        ccc = gauss2(Z(i), sigma(1), mu(1));
        ddd = gauss2(Z(i), sigma(2), mu(2));
        P1t(i) = aaa * alpha(1);
        P2t(i) = bbb * alpha(2);
        Pt = P1t(i) + P2t(i);
        loglikelihood = loglikelihood + log(Pt);
        W1t(i) = P1t(i) / Pt;
        W2t(i) = P2t(i) / Pt;
    end

    if abs(loglikelihood-oldlikelihood) < e
        break;
    else
        oldlikelihood = loglikelihood;
    end
    
    %M step
    nW1t = 0;
    nW2t = 0;
    for i = 1 : length(Z)
        nW1t = nW1t + W1t(i);
        nW2t = nW2t + W2t(i);
    end
    alpha(1) = nW1t / length(Z);
    alpha(2) = nW2t / length(Z);
    
    nmu1 = 0;
    nmu2 = 0;
    for i = 1 : length(Z)
        nmu1 = nmu1 + Z(i)*W1t(i);
        nmu2 = nmu2 + Z(i)*W2t(i);
    end
    mu(1) = nmu1 / (length(Z) * alpha(1));
    mu(2) = nmu2 / (length(Z) * alpha(2));
    
    nsigma1 = 0;
    nsigma2 = 0;
    for i = 1 : length(Z)
        nsigma1 = nsigma1 + alpha(1) * (Z(i) - mu(1))^2;
        nsigma2 = nsigma2 + alpha(2) * (Z(i) - mu(2))^2;
    end
    sigma(1) = nsigma1 / (length(Z) * alpha(1));
    sigma(2) = nsigma2 / (length(Z) * alpha(2));
end   
 q = [mu,sigma,alpha, iter];   
    


% q = [x,sig2];
% for k = 1:K
%      p = Phit(Z(k, u,));
% %       + Zshort * Pshort(z, x, m)...
% %       + Zmax * Pmax(z, x, m)...
% %       + Zrand * Prand(z, x, m);
%     q = q * p;
% end



%%
% function p = Phit(z, sig2, u)
% sig2 = 0.5;
% p = gauss2(z, sig2, u);

%%
% function p = Pshort(z, x, m)


%%
% function p = Pmax(z, x, m)


%%
% function p = Prand(z, x, m)


