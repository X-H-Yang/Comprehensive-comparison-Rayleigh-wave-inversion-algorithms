function [y, dispersion_num] = Jg(theta, layers_num, model_dispersion, ...
    Vp, den, f,modes_num_vec,index_vec_all)

% Acknowledgement: The forward modeling program used to generate 
%                  theoretical Rayleigh wave dispersion curves in this 
%                  study was obtained from the  website 
%                  (https://github.com/eespr/MuLTI) provided by 
%                  Killingbeck et al. (2018)
%
% Killingbeck et al. (2018): Killingbeck, S. F., Livermore, P. W., 
%                            Booth, A. D., & West, L. J. (2018). Multimodal 
%                            layered transdimensional inversion of seismic 
%                            dispersion curves with depth constraints. 
%                            Geochemistry, Geophysics, Geosystems, 19(12), 
%                            4957-4971.

dispersion_num = max(size(model_dispersion,1),size(model_dispersion,2));
N = layers_num-1;
model_parameter = [N theta Vp den];

pre_dispersion = [];

h = [theta(1:layers_num-1) 0];
Vs = theta(layers_num:end);

try
    out = gpdc(h,Vp,Vs,den,'fV',f);
    out2 = rdivide(1, out(:, 2:end));
    
    for jj = 1:1:length(modes_num_vec)
        temp = modes_num_vec(jj);
        point_temp = index_vec_all{temp};
        star_point = point_temp(1);
        end_point = point_temp(2);
        pre_dispersion = [pre_dispersion out2(star_point:end_point,temp)'];
    end
    pre_dispersion(isnan(pre_dispersion)) = 0;
    
    y = sum((pre_dispersion - model_dispersion).^2)/dispersion_num;
catch
    y = Inf;
end

time_recorded = toc;
global time_global_index Obj_global_index

time_global_index = [time_global_index;time_recorded];
Obj_global_index = [Obj_global_index;y];

end