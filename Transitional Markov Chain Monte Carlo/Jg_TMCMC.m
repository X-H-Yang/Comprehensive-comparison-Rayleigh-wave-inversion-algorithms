function [isSameLast,y] = Jg_TMCMC(theta,N_earth,f,model_true_dispersion,modes_num_vec,index_vec_all,den_true,Vp_true)
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

model_true_dispersion = model_true_dispersion(:)';

model_dispersion_R = [];

dispersion_num = max(size(model_true_dispersion,1),size(model_true_dispersion,2));
h = theta(1:N_earth);
h2 = [h 0];
Vp = Vp_true; % primary wave velocity (all layers)
den = den_true; % density (all layers)
Vs = theta(end-N_earth:end);

try
    out = gpdc(h2,Vp,Vs,den,'fV',f);
    out2 = rdivide(1, out(:, 2:end));
    
    for jj = 1:1:length(modes_num_vec)
        temp = modes_num_vec(jj);
        point_temp = index_vec_all{temp};
        star_point = point_temp(1);
        end_point = point_temp(2);
        model_dispersion_R = [model_dispersion_R out2(star_point:end_point,temp)'];
    end
    model_dispersion_R(isnan(model_dispersion_R)) = 0;
    
    y = sum((model_dispersion_R - model_true_dispersion).^2)/dispersion_num;

    isSameLast=0;
catch
    isSameLast = 1;
    y = 100000;
end

end

