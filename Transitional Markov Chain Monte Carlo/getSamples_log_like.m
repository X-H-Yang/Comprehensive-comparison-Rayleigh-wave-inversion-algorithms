function [log_like_all,Vs_profile_all] = getSamples_log_like(dispersions_R_true,...
    modes_num_vec,f_vector,Vp_true,den_true,Vs_profile_lower,Vs_profile_upper,samples_N,index_vec_all)
% function purpose: generate samples of Rayleigh wave phase velocities
%                   of multi-modes : generate training dataset
%
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

log_like_all = zeros(samples_N,1);
N_true = length(den_true)-1;
dispersions_R_true = dispersions_R_true(:)';
Vs_profile_all = [];
dispersion_num = length(dispersions_R_true);

%% generate samples

x_low = Vs_profile_lower;
x_up = Vs_profile_upper;

num_unknown = length(x_low);
Y_raw = zeros(samples_N,num_unknown);
X_raw = zeros(samples_N,length(dispersions_R_true));

i = 1;

while i <= samples_N
    
    for j = 1:num_unknown
        Y_raw(i,j) = x_low(j) + (x_up(j)-x_low(j))*rand(1,1);
    end
    
    
    h0 = Y_raw(i,1:N_true);
    h = [Y_raw(i,1:N_true) 0];
    Vp = Vp_true;
    Vs = Y_raw(i,N_true+1:2*N_true+1);
    den = den_true;
    
    model_dispersion_R = [];
    
    try
        out = gpdc(h,Vp,Vs,den,'fV',f_vector);
        out2 = rdivide(1, out(:, 2:end));
        
        for jj = 1:1:length(modes_num_vec)
            temp = modes_num_vec(jj);
            point_temp = index_vec_all{temp};
            star_point = point_temp(1);
            end_point = point_temp(2);
            model_dispersion_R = [model_dispersion_R out2(star_point:end_point,temp)'];
        end
        model_dispersion_R(isnan(model_dispersion_R)) = 0;
        X_raw(i,:) = model_dispersion_R;
        
        pro_temp = [h0 Vs];
        Vs_profile_all = [Vs_profile_all;pro_temp];
        
        sigma2 = sum((model_dispersion_R - dispersions_R_true).^2)/dispersion_num;
        sigma = sqrt(sigma2);
        ln_likelihood = -dispersion_num*log(2*pi)/2 - dispersion_num ...
                *log(sigma) - dispersion_num/2;
        log_like_all(i) = ln_likelihood;
        
        i = i+1;
    catch
%         X_raw(i,:) = X_raw(i-1,:);
%         Y_raw(i,:) = Y_raw(i-1,:);
    end
    
end
Vs_profile_all = Vs_profile_all';
end

