function [samples_x,samples_y] = getSamples_RW_fieldData_sub_2( samples_N,...
    modes_num_vec,index_vec_all,f_vector,Vp_true,den_true,...
    Vs_profile_lower,Vs_profile_upper,dispersions_R_true)
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

assert(size(modes_num_vec,1)==1 || size(modes_num_vec,2)==1,...
    'Error! modes_num_vec should be a vector!')
assert(max(modes_num_vec)<=5,...
    'Error! the maximum of modes_num_vec should be no larger than 5!')

%% Row vector
Vs_profile_lower = Vs_profile_lower(:)';
Vs_profile_upper = Vs_profile_upper(:)';

%% generate samples
% lower and upper boundary of samples
x_low = Vs_profile_lower;
x_up = Vs_profile_upper;
% samples of 1st generation
num_unknown = length(x_low);
Y_raw = zeros(samples_N,num_unknown);
X_raw = zeros(samples_N,length(dispersions_R_true));

N_true = length(Vp_true) - 1;

i = 1;
while i <= samples_N
    for j = 1:num_unknown
        Y_raw(i,j) = x_low(j) + (x_up(j)-x_low(j))*rand(1,1);
    end

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
        
        i = i+1;
    catch

    end
    
end
samples_x = X_raw;
samples_y = Y_raw;


end

