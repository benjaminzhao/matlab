%%clear the way
% author : liuweimin
% ictcas
close all;
clear;
clc;

%% settings

M=3;        % number of Gaussian
N=1000;     % total number of data samples

th=1e-3;    % convergent threshold
Nit=2000;    % maximal iteration
Nrep=10;     % number of repetation to find global maximal

K=2;        % demention of output signal

pi=3.141592653589793; % in case it is overwriten by smae name variable

cond_num =1000;  % prevent the singular covariance matrix in simulation data

plot_flag=1;
print_flag=1;

%% paramethers for random signal genrator

% random parameters for M Gaussian signals
 mu_real = randn(K,M); % mean
 
cov_real =zeros(K,K,M);  % covariance matrix
covd_real=zeros(K,K,M);  % covariance matrix decomposition
for cm=1:M
    while 1
        covd_real(:,:,cm)=randn(K,K);
        cov_real(:,:,cm)=covd_real(:,:,cm)*covd_real(:,:,cm)';
        if cond(cov_real(:,:,cm))>cond_num
            continue;
        else
            break;
        end
    end
end

% probablilty of a channel being selected
a_real = abs(randn(M,1));
a_real = a_real/sum(a_real);    % normlize

if print_flag==1
      a_real;%类别概率的真值
     mu_real;%均值的真值
    cov_real;%协方差的真值
end

%% generate random sample of Gaussian vectors

%m=randdist(1,N,[1:M],a_real);   % selector

rand_num_a=rand();%产生随机数
rand_num_b=rand();%产生随机数
while rand_num_a>=rand_num_b
    rand_num_a=rand();%产生随机数
    rand_num_b=rand();%产生随机数
end 

x=randn(K,N);
for c=1:round(rand_num_a*N)
    %sel=m(c);
    x(:,c)=covd_real(:,:,1)*x(:,c)+mu_real(1);
end
for c=round(rand_num_a*N)+1:round(rand_num_b*N)
    %sel=m(c);
    x(:,c)=covd_real(:,:,2)*x(:,c)+mu_real(2);
end
for c=round(rand_num_b*N)+1:(N)
    %sel=m(c);
    x(:,c)=covd_real(:,:,3)*x(:,c)+mu_real(3);
end

%% EM Algorothm

% loop
f_best=-inf;
for crep=1:Nrep
    c=1;    
    
    % initial values of parameters for EM
    a=abs(randn(M,1));  % randomly generated 
    a=a/sum(a); % normlize, such that sum(a_EM)=1
    mu=randn(K,M);
    cov =zeros(K,K,M);  % covariance matrix
    covd=zeros(K,K,M);  % covariance matrix decomposition
    for cm=1:M
        while 1
            covd(:,:,cm)=randn(K,K);
            cov(:,:,cm)=covd(:,:,cm)*covd(:,:,cm)';
            if cond(cov(:,:,cm))>cond_num
                continue;
            else
                break;
            end
        end
    end

    % iteration to find local maxima
    break_flag=0;
    while 1
          a_old=  a;
         mu_old= mu;
        cov_old=cov;
        
        fprintf(1,'calculating probability pmx...\n');
        pause(0);
        % pmx(m,x|param)
        pmx=zeros(M,N);
        for cm=1:M
            cov_cm=cov(:,:,cm);
            if cond(cov_cm) > cond_num
                break_flag=1;
            end
            inv_cov_cm=inv(cov_cm);
            mu_cm=mu(:,cm);
            for cn=1:N
                %p_cm=exp(-0.5*(x(:,cn)-mu_cm)'*inv_cov_cm*(x(:,cn)-mu_cm));
                p_cm=a(cm,:)*exp(-0.5*(x(:,cn)-mu_cm)'*inv_cov_cm*(x(:,cn)-mu_cm));%这里加上a
                pmx(cm,cn)=p_cm;
            end
            pmx(cm,:)=pmx(cm,:)/sqrt(det(cov_cm));
        end
        pmx=pmx*(2*pi)^(-K/2);

        fprintf(1,'calculating conditional probability, p...\n');
        pause(0);
        % conditional probability p(m|x,param) for estimated parameters
        p=pmx./kron(ones(M,1),sum(pmx));
    
        fprintf(1,'updating parametres\n');
        pause(0);
        a = 1/N*sum(p')';    
        mu = 1/N*x*p'*diag(1./a);
        for cm=1:M
            a_cm=a(cm);
            mu_cm=mu(:,cm);             
            tmp=x-kron(ones(1,N),mu_cm);
            cov(:,:,cm)=1/N*(kron(ones(K,1),p(cm,:)).*tmp)*tmp'*diag(1./a_cm);
           
        end
         
        t=max([norm(a_old(:)-a(:))/norm(a_old(:));
               norm(mu_old(:)-mu(:))/norm(mu_old(:));
               norm(cov_old(:)-cov(:))/norm(cov_old(:))]);
        if print_flag==1
            fprintf('c=%04d: t=%f\n',c,t);
            c=c+1;
        end
        
        if t<th
            break;
        end
    
        if c>Nit
            disp('reach maximal iteration')
            break;
        end
        
        if break_flag==1
            disp('***** break on singular covariance matrix *****');
            break;
        end
    end

    f=sum(log(sum(pmx.*kron(ones(1,N),a))));
    if f>f_best
          a_best=a;
         mu_best=mu;
        cov_best=cov;
          f_best=f;
    end
end
!echo 真实的rand_num_a
rand_num_a
!echo 真实的rand_num_b
rand_num_b-rand_num_a
!echo 真实的rand_num_c
1-rand_num_b
!echo 迭代出来的a_best
a_best
!echo 真实的均值矩阵
mu_real
!echo 迭代出来的均值矩阵
mu_best
for cs=1:M
  !echo 真实的协方差矩阵
  cov_real(:,:,cs)
  !echo 迭代出来的协方差矩阵
  cov_best(:,:,cs)
end

%% plot all
% for 2D (K=2) only
x1_vect=-1:0.02:1;
x2_vect=-1:0.02:1;
px=zeros(length(x1_vect), length(x2_vect));
for c1=1:length(x1_vect)
    for c2=1:length(x2_vect)
        for cm=1:3
            cov_real_cm=cov_real(:,:,cm);
            mu_real_cm=mu_real(:,cm);
            a_real_cm=a_real(cm);
            x_cm=[x1_vect(c1);
                  x2_vect(c2)];
            pm=a_real_cm*(2*pi)^(-0.5*K)*det(cov_real_cm)^(-0.5)*exp(-0.5*x_cm'*inv(cov_real_cm)*x_cm);
            px(c1,c2)=px(c1,c2)+pm;
        end
    end
end

px_hat=zeros(length(x1_vect), length(x2_vect));
for c1=1:length(x1_vect)
    for c2=1:length(x2_vect)
        for cm=1:3
            cov_cm=cov(:,:,cm);
            mu_cm=mu(:,cm);
            a_cm=a(cm);
            x_cm=[x1_vect(c1);
                  x2_vect(c2)];
            pm=a_cm*(2*pi)^(-0.5*K)*det(cov_cm)^(-0.5)*exp(-0.5*x_cm'*inv(cov_cm)*x_cm);
            px_hat(c1,c2)=px_hat(c1,c2)+pm;
        end
    end
end

figure(1); clf; hold on; 

hold on
title('理论图');
xlabel('x');
ylabel('y');
mesh(x1_vect, x2_vect, px);
figure(2); clf; hold on; 
mesh(x1_vect, x2_vect, px_hat);
title('统计图');
xlabel('x');
ylabel('y');