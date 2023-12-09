function y = PSOobjfun(theta_SA, layers_num, model_dispersion, Vp, den, f,modes_num_vec,index_vec_all)

[y, num] = Jg(theta_SA, layers_num, model_dispersion, Vp, den, f,modes_num_vec,index_vec_all);

end

