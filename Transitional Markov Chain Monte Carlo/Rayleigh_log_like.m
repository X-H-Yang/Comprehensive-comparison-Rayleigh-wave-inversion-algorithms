function [isSameLast,ln_likelihood] = Rayleigh_log_like(theta,N_earth,f,model_true_dispersion,modes_num_vec,index_vec_all,den_true,Vp_true)

[isSameLast,sigma2] = Jg_TMCMC(theta,N_earth,f,model_true_dispersion,modes_num_vec,index_vec_all,den_true,Vp_true);
sigma = sqrt(sigma2);
dispersion_num = max(size(model_true_dispersion,1),size(model_true_dispersion,2));
ln_likelihood = -dispersion_num*log(2*pi)/2 - dispersion_num ...
        *log(sigma) - dispersion_num/2;

end

