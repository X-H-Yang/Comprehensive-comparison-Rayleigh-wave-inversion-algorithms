function y = BFGSobjfun(theta_BFGS, layers_num, model_dispersion, Vp, ...
    den, f,modes_num_vec,index_vec_all,Vs_profile_lower,Vs_profile_upper)

[y, num] = Jg(theta_BFGS, layers_num, model_dispersion, Vp, den, f,...
    modes_num_vec,index_vec_all,Vs_profile_lower,Vs_profile_upper);

end

