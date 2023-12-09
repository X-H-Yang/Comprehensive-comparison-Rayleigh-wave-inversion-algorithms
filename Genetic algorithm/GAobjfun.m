function y = GAobjfun(theta_GA, layers_num, model_dispersion, Vp, den, f,modes_num_vec,index_vec_all)

thetaNum = size(theta_GA,1);
y = zeros(thetaNum,1);
for i =1:1:thetaNum
    [y(i), num] = Jg(theta_GA(i,:), layers_num, model_dispersion, Vp, ...
        den, f,modes_num_vec,index_vec_all);
end

end

