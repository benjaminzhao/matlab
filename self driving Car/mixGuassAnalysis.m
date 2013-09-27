

%用两个高斯函数生成数据，这两个高斯密度函数分别为N(0,2)和N(20,3.5)，分别生成160和240个数据，也就是说其在混合模型中的权重分别为0.4 和 0.6.
%混合高斯的EM算法的matlab实现如下：
function [mu, sigma, phi] = mixGuassAnalysis( k, maxIteration, epsilon)
%UNTITLED Analysis the k guass distribution by the input matrix m
%sampleMatrix   the matrix of sample, in which each row represents a sample.
%k  the number of guass distriubtion
%maxIteration   the max times of iteration,default 1000.
%epsilon    the epsilon of loglikelihood,default 100.


%数据准备阶段：
%generate data from normal distribution.
mu1= 0;
sigma1= 2;
R1 = normrnd(mu1,sigma1,160,1);

mu2 = 20;
sigma2 = 3.5;
R2 = normrnd(mu2,sigma2,240,1);

%merge
R = [R1;R2];
%shuffel
r = randperm(size(R,1));
R=R(r,:);

sampleMatrix = R;

%check parameters
if nargin < 4
    epsilon = 0.1;
    if nargin < 3
        maxIteration = 1000;
        if nargin < 2
            error('MATLAB:mixGuassAnalysis:NotEnoughInputs',...
                'Not enough input arguments.  See mixGuassAnalysis.');
        end
    end
end

if k==1
    mu = mean(sampleMatrix);
    sigma = var(sampleMatrix);
    phi = 1;
    return;
end

%400,1
[sampleNum, dimensionality] = size(sampleMatrix);

%init k guass distribution 2,1
mu = zeros(k,dimensionality);
for i=1:1:dimensionality
    colVector = sampleMatrix(:,i);
    maxV = max(colVector);
    minV = min(colVector);
    mu(1,i) = minV;
    mu(k,i) = maxV;
    for j=2:1:k-1
        mu(j:i) =  minv+(j-1)*(maxV-minV)/(k-1);
    end
end

%2,1,1
sigma = zeros(k,dimensionality,dimensionality);
for i=1:1:k
    d = rand();
    sigma(i,:) = 10*d*eye(dimensionality);%init sigma
end

phi = zeros(1,k);
for i=1:1:k
    phi(1,i) = 1.0/k;
end

%the weight of sample i is generated by guass distribution j
%400,2
weight = zeros(sampleNum,k);

oldlikelihood = -inf;

%1 to 1000
for iter=1:maxIteration
    loglikelihood = 0;
    %E-step
    for i=1:1:sampleNum
        for j = 1:1:k
            weight(i,j)=mvnpdf(sampleMatrix(i,:),mu(j,:),reshape(sigma(j,:),dimensionality,dimensionality))*phi(j);
        end
        
        sum = 0;
        for j = 1:1:k
            sum = sum+weight(i,j);
        end
        
        loglikelihood = loglikelihood + log(sum);
        
        % normalize
        for j = 1:1:k
            weight(i,j)=weight(i,j)/sum;
        end
    end
     
    if abs(loglikelihood-oldlikelihood)<epsilon
        break;
    else
        oldlikelihood = loglikelihood;
    end

    %M-step
    %update phi
    for i=1:1:k
        sum = 0;
        for j=1:1:sampleNum
            sum = sum+weight(j,i);
        end
        phi(i) = sum/sampleNum;
    end
    
    %update mu
    for i=1:1:k
        sum = zeros(1,dimensionality);
        for j=1:1:sampleNum
            sum =  sum+weight(j,i)*sampleMatrix(j,:);
        end
        
        mu(i,:) =  sum/(phi(i)*sampleNum);
    end
    
    %update sigma
    for i=1:1:k
        sum = zeros(dimensionality,dimensionality);
        for j=1:1:sampleNum
            sum = sum+ weight(j,i)*(sampleMatrix(j,:)-mu(i,:))'*(sampleMatrix(j,:)-mu(i,:));
        end
        sigma(i,:) = sum/(phi(i)*sampleNum);
    end
    
end
sigma = sqrt(sigma);

