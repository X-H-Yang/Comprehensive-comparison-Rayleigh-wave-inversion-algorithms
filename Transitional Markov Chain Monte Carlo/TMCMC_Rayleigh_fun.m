function [x,ln_S,log_like,x_stage,p_stage,time_index] = TMCMC_Rayleigh_fun(log_like_fun,model_true_dispersion, ...
    N_earth,f,x_low,x_up,N,modes_num_vec,Vp_true,den_true,index_vec_all)
% Acknowledgement: the Transitional Markov Chain Monte Carlo (TMCMC) codes 
%                  available on the website 
%   (https://www.researchgate.net/publication/288345281_TMCMC_matlab_codes, 
%                  Ching and Chen 2007) were also applied for the 
%                  accomplishment of this code package.
%
% Ching and Chen (2007): Transitional Markov chain Monte Carlo method for 
%                        Bayesian model updating, model class selection, 
%                        and model averaging. Journal of engineering 
%                        mechanics, 133(7), 816-832.

%
%    log_like_fun: filename for the m-file that computes the log
%                  likelihood.  The m-file must have the format log_like = f(x,y),
%                  where x contains the uncertain variable and y contains
%                  other variables necessary to compute log_like
%    y: other variables necessary to compute the log likelihood
%    x_low & x_up: lower and upper bounds for the uncertain variables x
%    N: number of TMCMC samples (per stage)
%    x: resulting TMCMC samples
%    ln_S: resulting log evidence 
%
%    Note: this TMCMC function assumes that the prior PDF for x is an 
%          uniform PDF defined by x_low & x_up
%

%% beta setting
b = 0.2;
%% initialing
p = 0; ln_S = 0;
%% drawing prior samples
[log_like,x] = getSamples_log_like(model_true_dispersion,...
    modes_num_vec,f,Vp_true,den_true,x_low,x_up,N,index_vec_all);

time_index = [];

x_stage{1} = x;
stage = 1;
nnn = length(x_low);

%% compute the weights
% adaptively choose p
p_index = 0;
while p<1,
    limit_value = 2;
    if p > 0.3
        limit_value = 0.5;
    end
    if p > 0.5
        limit_value = 0.3;
    end
    low_p = p; up_p = 2; old_p = p;
    while up_p-low_p > 1e-6
        current_p = (low_p + up_p)/2;
        for j = 1:1:numel(log_like)
            if log_like(j) == Inf
                log_like(j) = 400000;
            end
        end
        temp = exp((current_p-p)*(log_like-max(log_like)+2));
        cov_temp = std(temp)/mean(temp);
        if cov_temp > limit_value
            up_p = current_p;
        else
            low_p = current_p;
        end
    end
    p = current_p;
    p_index = p_index + 1;
    if p < 1, disp(['p',num2str(p_index),' =  ',num2str(p)]); end
    p_stage(stage) = p;
    if p > 1, break; end % breakout if approaching the final stage
    weight = temp/sum(temp); % weights are normalized
    ln_S = ln_S+log(mean(temp))+(p-old_p)*max(log_like);
    %% compute the covariance matrix
    mu_x = x*weight; sigma = zeros(nnn,nnn);
    for i=1:N,
        sigma = sigma + b^2*weight(i)*(x(:,i)-mu_x)*(x(:,i)-mu_x)';
    end
    sqrt_s = sqrtm(sigma);
    %% do MCMC
    sam_ind = deterministicR((1:N),weight); current_x = x; current_log_like = log_like;
    for i=1:N,
        now_ind = sam_ind(i);
        x_c = current_x(:,now_ind) + sqrt_s*randn(nnn,1);
        theta = x_c';
%         log_like_c = feval(log_like_fun,theta,N_earth,f,model_true_dispersion);
        isSameLast = 0;
        [isSameLast,log_like_c] = feval(log_like_fun,theta,N_earth,f,model_true_dispersion,modes_num_vec,index_vec_all,den_true,Vp_true);
%         i
        if(isSameLast)
            x_c = x(:,i);
            log_like_c = log_like(i);
        end
        r = exp(p*(log_like_c - current_log_like(now_ind)));
        if r > rand & sum(x_c < x_up)==nnn & sum(x_c > x_low)==nnn,
            x(:,i) = x_c; current_x(:,now_ind) = x_c; current_log_like(now_ind) = log_like_c; log_like(i) = log_like_c;
        else
            x(:,i) = current_x(:,now_ind); log_like(i) = current_log_like(now_ind);
        end
    end
    stage = stage+1;
    x_stage{stage} = x;
    
    temp_time = toc;
    time_index = [time_index;temp_time];
end
p_stage(stage) = 1;
%% finalize the final stage
for j = 1:1:numel(log_like)
    if log_like(j) == Inf
        log_like(j) = 4000;
    end
end
temp = exp((1-old_p)*(log_like-max(log_like)));
weight = temp/sum(temp);
ln_S = ln_S+log(mean(temp))+(1-old_p)*max(log_like);
mu_x = x*weight; sigma = zeros(nnn,nnn);
for i=1:N,
    sigma = sigma + b^2*weight(i)*(x(:,i)-mu_x)*(x(:,i)-mu_x)';
end
sqrt_s = sqrtm(sigma);
sam_ind = deterministicR((1:N),weight); current_x = x; current_log_like = log_like;
for i=1:N,
    now_ind = sam_ind(i);
    x_c = current_x(:,now_ind) + sqrt_s*randn(nnn,1);
    theta = x_c';
%     log_like_c = feval(log_like_fun,theta,N_earth,f,model_true_dispersion);
    isSameLast = 0;
    [isSameLast,log_like_c] = feval(log_like_fun,theta,N_earth,f,model_true_dispersion,modes_num_vec,index_vec_all,den_true,Vp_true);
    if(isSameLast)
        x_c = x(:,i);
        log_like_c = log_like(i);
    end
    r = exp(log_like_c - current_log_like(now_ind));
    if r > rand & sum(x_c < x_up)==nnn & sum(x_c > x_low)==nnn,
        x(:,i) = x_c; current_x(:,now_ind) = x_c; current_log_like(now_ind) = log_like_c; log_like(i) = log_like_c;
    else
        x(:,i) = current_x(:,now_ind); log_like(i) = current_log_like(now_ind);
    end
end
stage = stage+1;
x_stage{stage} = x;

temp_time = toc;
time_index = [time_index;temp_time];
