function gradValue = my_gradFun(theta_BFGS,layers_num, ...
    model_dispersion_true, Vp, den, f,modes_num_vec,index_vec_all)
% Acknowledgement: The forward modeling program used to generate 
%                  theoretical Rayleigh wave dispersion curves in this 
%                  code package was obtained from the  website 
%                  (https://github.com/eespr/MuLTI) provided by 
%                  Killingbeck et al. (2018)
%
% Killingbeck et al. (2018): Killingbeck, S. F., Livermore, P. W., 
%                            Booth, A. D., & West, L. J. (2018). Multimodal 
%                            layered transdimensional inversion of seismic 
%                            dispersion curves with depth constraints. 
%                            Geochemistry, Geophysics, Geosystems, 19(12), 
%                            4957-4971.

theta_BFGS = theta_BFGS(:)';

h = [theta_BFGS(1:layers_num-1) 0];
Vs = theta_BFGS(layers_num:end);

pre_dispersion_BFGS = [];
try
    out = gpdc(h,Vp,Vs,den,'fV',f);
    out2 = rdivide(1, out(:, 2:end));
    
    for jj = 1:1:length(modes_num_vec)
        temp = modes_num_vec(jj);
        point_temp = index_vec_all{temp};
        star_point = point_temp(1);
        end_point = point_temp(2);
        pre_dispersion_BFGS = [pre_dispersion_BFGS out2(...
            star_point:end_point,temp)'];
    end
    pre_dispersion_BFGS(isnan(pre_dispersion_BFGS)) = 0;
    
catch
    pre_dispersion_BFGS = model_dispersion_true;
    
end







m = length(model_dispersion_true);
n = 2*layers_num-1;

d_factor = 0.0001;

Jaco_matrix = zeros(m,n);

for i = 1:n
    theta_new = theta_BFGS;
    d_theta = d_factor*theta_BFGS(i);
    
    theta_new(i) = theta_new(i)+d_theta;
    
    h = [theta_new(1:layers_num-1) 0];
    Vs = theta_new(layers_num:end);
    
    pre_dispersion = [];
    
    try
        out = gpdc(h,Vp,Vs,den,'fV',f);
        out2 = rdivide(1, out(:, 2:end));
        
        for jj = 1:1:length(modes_num_vec)
            temp = modes_num_vec(jj);
            point_temp = index_vec_all{temp};
            star_point = point_temp(1);
            end_point = point_temp(2);
            pre_dispersion = [pre_dispersion out2(...
                star_point:end_point,temp)'];
        end
        pre_dispersion(isnan(pre_dispersion)) = 0;
        
    catch
        pre_dispersion = pre_dispersion_BFGS;
        
    end
    
    temp_2 = (pre_dispersion - pre_dispersion_BFGS)/d_theta;
    
    for j = 1:m
        Jaco_matrix(j,i) = temp_2(j);
    end
    
end



gradValue = -2*Jaco_matrix'*(pre_dispersion_BFGS-...
    model_dispersion_true)'/length(model_dispersion_true);

% gradValue = gradValue(:)';
gradValue = -gradValue(:); %%% negtive sign is the key point

end

